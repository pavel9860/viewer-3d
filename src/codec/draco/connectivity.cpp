#include "codec/draco/connectivity.hpp"

#include <algorithm>
#include <limits>
#include <unordered_map>

namespace cadcore::codec::draco {

// ---- attribute corner table (declared in corner_table.hpp) ---------------------------------------------

void AttributeCornerTable::add_seam_edge(std::int32_t c) {
    edge_on_seam[std::size_t(c)] = true;
    vertex_on_seam[std::size_t(pos->Vertex(Next(c)))] = true;
    vertex_on_seam[std::size_t(pos->Vertex(Previous(c)))] = true;
    const std::int32_t opp = pos->Opposite(c);
    if (opp != kInvalid) {
        edge_on_seam[std::size_t(opp)] = true;
        vertex_on_seam[std::size_t(pos->Vertex(Next(opp)))] = true;
        vertex_on_seam[std::size_t(pos->Vertex(Previous(opp)))] = true;
    }
}

bool AttributeCornerTable::recompute_vertices() {
    left_most.clear();
    std::int32_t num_new = 0;
    for (std::int32_t v = 0; v < pos->num_vertices(); ++v) {
        const std::int32_t c = pos->LeftMostCorner(v);
        if (c == kInvalid) continue;  // isolated vertex
        std::int32_t first_vert = num_new++;
        std::int32_t first_c = c;
        if (vertex_on_seam[std::size_t(v)]) {
            // Walk left to the first seam so the vertex split starts at a seam edge.
            std::int32_t act = SwingLeft(first_c);
            while (act != kInvalid) {
                first_c = act;
                act = SwingLeft(act);
                if (act == c) return false;
            }
        }
        vertex[std::size_t(first_c)] = first_vert;
        left_most.push_back(first_c);
        std::int32_t act = pos->SwingRight(first_c);
        while (act != kInvalid && act != first_c) {
            if (IsCornerOppositeToSeamEdge(Next(act))) {
                first_vert = num_new++;
                left_most.push_back(act);
            }
            vertex[std::size_t(act)] = first_vert;
            act = pos->SwingRight(act);
        }
    }
    return true;
}

// ---- sequential -------------------------------------------------------------------------------------------

void decode_sequential_connectivity(Stream& s, Connectivity& out) {
    std::uint32_t num_faces, num_points;
    if (s.ver < version(2, 2)) {
        num_faces = s.u32();
        num_points = s.u32();
    } else {
        num_faces = s.varint32();
        num_points = s.varint32();
    }
    if (std::uint64_t(num_faces) > 0xffffffffu / 3 || num_faces > s.remaining() / 3) bad("sequential face count exceeds the data");
    const int method = s.u8();
    out.num_points = num_points;
    out.faces.resize(std::size_t(num_faces) * 3);
    if (method == 0) {  // compressed: rANS-coded signed deltas
        std::vector<std::uint32_t> sym(out.faces.size());
        decode_symbols(s, std::uint32_t(sym.size()), 1, sym.data());
        std::int64_t last = 0;
        for (std::size_t i = 0; i < sym.size(); ++i) {
            const std::int64_t diff = sym[i] >> 1;
            last += (sym[i] & 1) ? -diff : diff;
            if (last < 0 || last >= num_points) bad("sequential index out of range");
            out.faces[i] = std::uint32_t(last);
        }
    } else {
        for (auto& f : out.faces) {
            std::uint32_t v;
            if (num_points < 256) v = s.u8();
            else if (num_points < (1u << 16)) v = s.u16();
            else if (num_points < (1u << 21) && s.ver >= version(2, 2)) v = s.varint32();
            else v = s.u32();
            if (v >= num_points) bad("sequential index out of range");
            f = v;
        }
    }
}

// ---- EdgeBreaker ---------------------------------------------------------------------------------------------

namespace {

constexpr int TOPOLOGY_C = 0, TOPOLOGY_S = 1, TOPOLOGY_L = 3, TOPOLOGY_R = 5, TOPOLOGY_E = 7;
constexpr int kSymbolToTopology[5] = {TOPOLOGY_C, TOPOLOGY_S, TOPOLOGY_L, TOPOLOGY_R, TOPOLOGY_E};
constexpr int STANDARD = 0, PREDICTIVE = 1, VALENCE = 2;
constexpr int RIGHT_FACE_EDGE = 1;

struct SplitEvent {
    std::uint32_t split_symbol_id = 0, source_symbol_id = 0, source_edge = 0;
};

std::uint32_t count_field(Stream& s) { return s.ver < version(2, 0) ? s.u32() : s.varint32(); }

// Topology split events (and the legacy hole events, which carry no decodable information).
std::vector<SplitEvent> decode_events(Stream& s, std::uint32_t num_faces) {
    std::vector<SplitEvent> ev;
    const std::uint32_t n = count_field(s);
    if (n > 0) {
        if (n > num_faces) bad("more topology splits than faces");
        if (s.ver < version(1, 2)) {
            for (std::uint32_t i = 0; i < n; ++i) {
                SplitEvent e;
                e.split_symbol_id = s.u32();
                e.source_symbol_id = s.u32();
                e.source_edge = s.u8() & 1;
                ev.push_back(e);
            }
        } else {
            std::uint32_t last = 0;
            for (std::uint32_t i = 0; i < n; ++i) {
                SplitEvent e;
                e.source_symbol_id = s.varint32() + last;
                const std::uint32_t delta = s.varint32();
                if (delta > e.source_symbol_id) bad("topology split before its source");
                e.split_symbol_id = e.source_symbol_id - delta;
                last = e.source_symbol_id;
                ev.push_back(e);
            }
            s.start_bits();
            for (auto& e : ev) e.source_edge = s.bits(s.ver < version(2, 2) ? 2 : 1) & 1;
            s.end_bits();
        }
    }
    std::uint32_t holes = 0;
    if (s.ver < version(2, 0)) holes = s.u32();
    else if (s.ver < version(2, 1)) holes = s.varint32();
    for (std::uint32_t i = 0; i < holes; ++i) {
        if (s.ver < version(1, 2)) s.i32();
        else s.varint32();
    }
    return ev;
}

// Symbol source for the three traversal flavours.
class Traversal {
public:
    Traversal(int type, const CornerTable& ct) : type_(type), ct_(ct) {}

    // Parses everything between the split events and the attribute data.
    void start(Stream& s, std::uint32_t num_vertices, std::uint32_t num_faces, std::uint32_t num_attribute_data) {
        const bool legacy = s.ver < version(2, 2);
        auto sized_bits = [&](BitSpan& span) {
            const std::uint64_t n = legacy ? s.u64() : s.varint64();
            if (n > s.remaining()) bad("EdgeBreaker bit buffer exceeds the stream");
            span = BitSpan(s.bytes(std::size_t(n)), std::size_t(n));
        };
        if (type_ != VALENCE || legacy) sized_bits(symbols_);
        if (legacy) sized_bits(start_face_bits_);
        else start_faces_.start(s);
        seams_.resize(num_attribute_data);
        for (auto& d : seams_) d.start(s);
        if (type_ == PREDICTIVE) {
            const std::int32_t splits = s.i32();
            if (splits < 0 || std::uint32_t(splits) >= num_vertices) bad("invalid predictive split count");
            valences_.assign(num_vertices, 0);
            prediction_.start(s);
        } else if (type_ == VALENCE) {
            if (legacy) {
                const std::uint32_t splits = count_field(s);
                if (splits >= num_vertices) bad("invalid valence split count");
                if (s.i8() != 0) unsupported("valence EdgeBreaker mode other than 2..7");
            }
            valences_.assign(num_vertices, 0);
            context_symbols_.resize(6);
            context_counters_.assign(6, 0);
            for (int i = 0; i < 6; ++i) {
                const std::uint32_t n = s.varint32();
                if (n > num_faces) bad("valence context larger than the face count");
                if (n > 0) {
                    context_symbols_[std::size_t(i)].resize(n);
                    decode_symbols(s, n, 1, context_symbols_[std::size_t(i)].data());
                    context_counters_[std::size_t(i)] = std::int32_t(n);
                }
            }
        }
        legacy_ = legacy;
    }

    int symbol() {
        if (type_ == VALENCE) {
            if (active_context_ != -1) {
                const std::int32_t k = --context_counters_[std::size_t(active_context_)];
                if (k < 0) bad("valence context exhausted");
                const std::uint32_t id = context_symbols_[std::size_t(active_context_)][std::size_t(k)];
                if (id > 4) bad("invalid valence symbol");
                last_ = kSymbolToTopology[id];
            } else {
                last_ = legacy_ ? standard_symbol() : TOPOLOGY_E;
            }
            return last_;
        }
        if (type_ == PREDICTIVE && predicted_ != -1 && prediction_.next()) return last_ = predicted_;
        return last_ = standard_symbol();
    }

    void new_active_corner(std::int32_t corner) {
        if (type_ == STANDARD) return;
        const std::int32_t next = ct_.Next(corner), prev = ct_.Previous(corner);
        auto add = [&](std::int32_t c, int k) {
            const std::int32_t v = ct_.Vertex(c);
            if (v < 0 || std::size_t(v) >= valences_.size()) bad("valence vertex out of range");
            valences_[std::size_t(v)] += k;
        };
        switch (last_) {
            case TOPOLOGY_C:
            case TOPOLOGY_S: add(next, 1); add(prev, 1); break;
            case TOPOLOGY_R: add(corner, 1); add(next, 1); add(prev, 2); break;
            case TOPOLOGY_L: add(corner, 1); add(next, 2); add(prev, 1); break;
            case TOPOLOGY_E: add(corner, 2); add(next, 2); add(prev, 2); break;
            default: break;
        }
        if (type_ == VALENCE) {
            const int val = valences_[std::size_t(ct_.Vertex(next))];
            active_context_ = std::clamp(val, 2, 7) - 2;
        } else if (last_ == TOPOLOGY_C || last_ == TOPOLOGY_R) {
            predicted_ = valences_[std::size_t(ct_.Vertex(next))] < 6 ? TOPOLOGY_R : TOPOLOGY_C;
        } else {
            predicted_ = -1;
        }
    }

    void merge_vertices(std::int32_t dest, std::int32_t src) {
        if (type_ == STANDARD) return;
        valences_[std::size_t(dest)] += valences_[std::size_t(src)];
    }

    bool start_face_configuration() { return legacy_ ? start_face_bits_.bits(1) != 0 : start_faces_.next(); }
    bool attribute_seam(std::size_t i) { return seams_[i].next(); }

private:
    int standard_symbol() {
        int sym = int(symbols_.bits(1));
        if (sym != TOPOLOGY_C) sym |= int(symbols_.bits(2)) << 1;
        return sym;
    }

    int type_;
    const CornerTable& ct_;
    bool legacy_ = false;
    BitSpan symbols_, start_face_bits_;
    BitDecoder start_faces_, prediction_;
    std::vector<BitDecoder> seams_;
    std::vector<int> valences_;
    std::vector<std::vector<std::uint32_t>> context_symbols_;
    std::vector<std::int32_t> context_counters_;
    int last_ = -1, active_context_ = -1, predicted_ = -1;
};

class EdgeBreaker {
public:
    EdgeBreaker(CornerTable& ct, Traversal& tr, std::vector<SplitEvent> splits, std::uint32_t max_vertices, bool remove_invalid)
        : ct_(ct), tr_(tr), splits_(std::move(splits)), max_vertices_(max_vertices), remove_invalid_(remove_invalid) {
        is_vert_hole.assign(max_vertices, true);
    }

    std::vector<bool> is_vert_hole;

    // Returns the number of connectivity vertices.
    std::int32_t run(std::int32_t num_symbols) {
        std::vector<std::int32_t> stack;
        std::unordered_map<std::int32_t, std::int32_t> split_active;
        std::vector<std::int32_t> invalid_vertices;
        std::int32_t num_faces = 0;
        auto need = [&](bool ok, const char* what) { if (!ok) bad(std::string("EdgeBreaker: ") + what); };
        for (std::int32_t symbol_id = 0; symbol_id < num_symbols; ++symbol_id) {
            const std::int32_t face = num_faces++;
            need(face < ct_.num_faces(), "more symbols than faces");
            bool check_split = false;
            const int symbol = tr_.symbol();
            const std::int32_t corner = 3 * face;
            if (symbol == TOPOLOGY_C) {
                need(!stack.empty(), "C without an active corner");
                const std::int32_t corner_a = stack.back();
                const std::int32_t vertex_x = ct_.Vertex(ct_.Next(corner_a));
                const std::int32_t corner_b = ct_.Next(ct_.LeftMostCorner(vertex_x));
                need(corner_a != corner_b && corner_b != kInvalid, "degenerate C");
                need(ct_.Opposite(corner_a) == kInvalid && ct_.Opposite(corner_b) == kInvalid, "C on an interior edge");
                set_opposite(corner_a, corner + 1);
                set_opposite(corner_b, corner + 2);
                const std::int32_t vert_a_prev = ct_.Vertex(ct_.Previous(corner_a));
                const std::int32_t vert_b_next = ct_.Vertex(ct_.Next(corner_b));
                need(vertex_x != vert_a_prev && vertex_x != vert_b_next, "degenerate C face");
                ct_.vertex[std::size_t(corner)] = vertex_x;
                ct_.vertex[std::size_t(corner + 1)] = vert_b_next;
                ct_.vertex[std::size_t(corner + 2)] = vert_a_prev;
                ct_.left_most[std::size_t(vert_a_prev)] = corner + 2;
                is_vert_hole[std::size_t(vertex_x)] = false;
                stack.back() = corner;
            } else if (symbol == TOPOLOGY_R || symbol == TOPOLOGY_L) {
                need(!stack.empty(), "R/L without an active corner");
                const std::int32_t corner_a = stack.back();
                need(ct_.Opposite(corner_a) == kInvalid, "R/L on an interior edge");
                std::int32_t opp, corner_l, corner_r;
                if (symbol == TOPOLOGY_R) { opp = corner + 2; corner_l = corner + 1; corner_r = corner; }
                else { opp = corner + 1; corner_l = corner; corner_r = corner + 2; }
                set_opposite(opp, corner_a);
                const std::int32_t nv = new_vertex();
                ct_.vertex[std::size_t(opp)] = nv;
                ct_.left_most[std::size_t(nv)] = opp;
                const std::int32_t vertex_r = ct_.Vertex(ct_.Previous(corner_a));
                ct_.vertex[std::size_t(corner_r)] = vertex_r;
                ct_.left_most[std::size_t(vertex_r)] = corner_r;
                ct_.vertex[std::size_t(corner_l)] = ct_.Vertex(ct_.Next(corner_a));
                stack.back() = corner;
                check_split = true;
            } else if (symbol == TOPOLOGY_S) {
                need(!stack.empty(), "S without an active corner");
                const std::int32_t corner_b = stack.back();
                stack.pop_back();
                if (auto it = split_active.find(symbol_id); it != split_active.end()) stack.push_back(it->second);
                need(!stack.empty(), "S without a second active corner");
                const std::int32_t corner_a = stack.back();
                need(corner_a != corner_b, "degenerate S");
                need(ct_.Opposite(corner_a) == kInvalid && ct_.Opposite(corner_b) == kInvalid, "S on an interior edge");
                set_opposite(corner_a, corner + 2);
                set_opposite(corner_b, corner + 1);
                const std::int32_t vertex_p = ct_.Vertex(ct_.Previous(corner_a));
                ct_.vertex[std::size_t(corner)] = vertex_p;
                ct_.vertex[std::size_t(corner + 1)] = ct_.Vertex(ct_.Next(corner_a));
                const std::int32_t vert_b_prev = ct_.Vertex(ct_.Previous(corner_b));
                ct_.vertex[std::size_t(corner + 2)] = vert_b_prev;
                ct_.left_most[std::size_t(vert_b_prev)] = corner + 2;
                std::int32_t corner_n = ct_.Next(corner_b);
                const std::int32_t vertex_n = ct_.Vertex(corner_n);
                tr_.merge_vertices(vertex_p, vertex_n);
                ct_.left_most[std::size_t(vertex_p)] = ct_.LeftMostCorner(vertex_n);
                const std::int32_t first = corner_n;
                while (corner_n != kInvalid) {
                    ct_.vertex[std::size_t(corner_n)] = vertex_p;
                    corner_n = ct_.SwingLeft(corner_n);
                    need(corner_n != first, "vertex merge loops");
                }
                ct_.left_most[std::size_t(vertex_n)] = kInvalid;
                if (remove_invalid_) invalid_vertices.push_back(vertex_n);
                stack.back() = corner;
            } else if (symbol == TOPOLOGY_E) {
                const std::int32_t v0 = new_vertex(), v1 = new_vertex(), v2 = new_vertex();
                ct_.vertex[std::size_t(corner)] = v0;
                ct_.vertex[std::size_t(corner + 1)] = v1;
                ct_.vertex[std::size_t(corner + 2)] = v2;
                ct_.left_most[std::size_t(v0)] = corner;
                ct_.left_most[std::size_t(v1)] = corner + 1;
                ct_.left_most[std::size_t(v2)] = corner + 2;
                stack.push_back(corner);
                check_split = true;
            } else {
                bad("invalid EdgeBreaker symbol " + std::to_string(symbol));
            }
            tr_.new_active_corner(stack.back());
            if (check_split) {
                const std::int32_t encoder_symbol_id = num_symbols - symbol_id - 1;
                while (!splits_.empty()) {
                    const SplitEvent& e = splits_.back();
                    if (e.source_symbol_id > std::uint32_t(encoder_symbol_id)) bad("topology split event missed");
                    if (e.source_symbol_id != std::uint32_t(encoder_symbol_id)) break;
                    const std::int32_t act_top = stack.back();
                    const std::int32_t new_active = e.source_edge == RIGHT_FACE_EDGE ? ct_.Next(act_top) : ct_.Previous(act_top);
                    split_active[num_symbols - std::int32_t(e.split_symbol_id) - 1] = new_active;
                    splits_.pop_back();
                }
            }
        }
        // Remaining active corners: holes or interior (closing) faces.
        while (!stack.empty()) {
            const std::int32_t corner = stack.back();
            stack.pop_back();
            if (!tr_.start_face_configuration()) continue;
            need(num_faces < ct_.num_faces(), "more faces than declared");
            const std::int32_t corner_a = corner;
            const std::int32_t vert_n = ct_.Vertex(ct_.Next(corner_a));
            const std::int32_t corner_b = ct_.Next(ct_.LeftMostCorner(vert_n));
            const std::int32_t vert_x = ct_.Vertex(ct_.Next(corner_b));
            const std::int32_t corner_c = ct_.Next(ct_.LeftMostCorner(vert_x));
            need(corner_b != kInvalid && corner_c != kInvalid, "open interior face");
            need(corner != corner_b && corner != corner_c && corner_b != corner_c, "degenerate interior face");
            need(ct_.Opposite(corner) == kInvalid && ct_.Opposite(corner_b) == kInvalid && ct_.Opposite(corner_c) == kInvalid,
                 "interior face on an interior edge");
            const std::int32_t vert_p = ct_.Vertex(ct_.Next(corner_c));
            const std::int32_t nc = 3 * num_faces++;
            set_opposite(nc, corner);
            set_opposite(nc + 1, corner_b);
            set_opposite(nc + 2, corner_c);
            ct_.vertex[std::size_t(nc)] = vert_x;
            ct_.vertex[std::size_t(nc + 1)] = vert_p;
            ct_.vertex[std::size_t(nc + 2)] = vert_n;
            for (int k = 0; k < 3; ++k) is_vert_hole[std::size_t(ct_.vertex[std::size_t(nc + k)])] = false;
        }
        need(num_faces == ct_.num_faces(), "face count differs from the header");
        // Compact vertices freed by S merges (only when no attribute seams refer to them).
        std::int32_t num_vertices = ct_.num_vertices();
        for (const std::int32_t invalid : invalid_vertices) {
            std::int32_t src = num_vertices - 1;
            while (src >= 0 && ct_.LeftMostCorner(src) == kInvalid) src = --num_vertices - 1;
            if (src < invalid) continue;
            std::vector<std::int32_t> corners;
            for_each_vertex_corner(ct_, ct_.LeftMostCorner(src), [&](std::int32_t c) { corners.push_back(c); });
            for (const std::int32_t c : corners) {
                need(ct_.Vertex(c) == src, "inconsistent vertex during compaction");
                ct_.vertex[std::size_t(c)] = invalid;
            }
            ct_.left_most[std::size_t(invalid)] = ct_.LeftMostCorner(src);
            ct_.left_most[std::size_t(src)] = kInvalid;
            is_vert_hole[std::size_t(invalid)] = is_vert_hole[std::size_t(src)];
            is_vert_hole[std::size_t(src)] = false;
            --num_vertices;
        }
        return num_vertices;
    }

private:
    void set_opposite(std::int32_t a, std::int32_t b) {
        ct_.opposite[std::size_t(a)] = b;
        ct_.opposite[std::size_t(b)] = a;
    }
    std::int32_t new_vertex() {
        if (std::uint32_t(ct_.num_vertices()) >= max_vertices_) bad("EdgeBreaker: more vertices than declared");
        return ct_.add_vertex();
    }

    CornerTable& ct_;
    Traversal& tr_;
    std::vector<SplitEvent> splits_;
    std::uint32_t max_vertices_;
    bool remove_invalid_;
};

}  // namespace

void decode_edgebreaker_connectivity(Stream& s, Connectivity& out) {
    out.edgebreaker = true;
    const int traversal_type = s.u8();
    if (traversal_type != STANDARD && traversal_type != PREDICTIVE && traversal_type != VALENCE)
        bad("unknown EdgeBreaker traversal type " + std::to_string(traversal_type));
    if (s.ver < version(2, 2)) count_field(s);  // num_new_vertices (not used for decoding)
    const std::uint32_t num_encoded_vertices = count_field(s);
    const std::uint32_t num_faces = count_field(s);
    if (num_faces > std::uint32_t(std::numeric_limits<std::int32_t>::max()) / 3) bad("too many faces");
    if (std::uint64_t(num_encoded_vertices) > std::uint64_t(num_faces) * 3) bad("more vertices than corners");
    const std::uint64_t nv = num_encoded_vertices;
    if (nv * (nv > 0 ? nv - 1 : 0) / 2 < 3ull * num_faces / 2) bad("too few vertices for the face count");
    const std::uint32_t num_attribute_data = s.u8();
    const std::uint32_t num_encoded_symbols = count_field(s);
    if (num_faces < num_encoded_symbols) bad("more symbols than faces");
    if (num_faces > num_encoded_symbols + num_encoded_symbols / 3) bad("too many faces for the symbol count");
    const std::uint32_t num_encoded_split_symbols = count_field(s);
    if (num_encoded_split_symbols > num_encoded_symbols) bad("more split symbols than symbols");
    const std::uint32_t max_vertices = num_encoded_vertices + num_encoded_split_symbols;

    out.table.reset(num_faces);
    out.attribute_data.resize(num_attribute_data);

    std::vector<SplitEvent> splits;
    std::size_t legacy_event_bytes = 0;
    if (s.ver < version(2, 2)) {
        const std::uint32_t connectivity_size = count_field(s);
        if (connectivity_size == 0 || connectivity_size > s.remaining()) bad("invalid encoded connectivity size");
        Stream ev(s.head() + connectivity_size, s.remaining() - connectivity_size);
        ev.ver = s.ver;
        splits = decode_events(ev, num_faces);
        legacy_event_bytes = ev.pos();
    } else {
        splits = decode_events(s, num_faces);
    }

    Traversal tr(traversal_type, out.table);
    tr.start(s, max_vertices, num_faces, num_attribute_data);
    EdgeBreaker eb(out.table, tr, std::move(splits), max_vertices, num_attribute_data == 0);
    const std::int32_t num_connectivity_verts = eb.run(std::int32_t(num_encoded_symbols));
    if (s.ver < version(2, 2)) s.skip(legacy_event_bytes);

    // Attribute seams, read face by face.
    std::vector<std::vector<std::int32_t>> seam_corners(num_attribute_data);
    if (num_attribute_data > 0) {
        const CornerTable& ct = out.table;
        for (std::int32_t c0 = 0; c0 < ct.num_corners(); c0 += 3) {
            const std::int32_t corners[3] = {c0, ct.Next(c0), ct.Previous(c0)};
            for (const std::int32_t c : corners) {
                const std::int32_t opp = ct.Opposite(c);
                if (opp == kInvalid) {
                    for (auto& sc : seam_corners) sc.push_back(c);
                    continue;
                }
                if (s.ver >= version(2, 1) && opp / 3 < c0 / 3) continue;  // decoded with the other face
                for (std::size_t i = 0; i < num_attribute_data; ++i)
                    if (tr.attribute_seam(i)) seam_corners[i].push_back(c);
            }
        }
    }
    for (std::size_t i = 0; i < num_attribute_data; ++i) {
        AttributeCornerTable& t = out.attribute_data[i].table;
        t.init(out.table);
        for (const std::int32_t c : seam_corners[i]) t.add_seam_edge(c);
        if (!t.recompute_vertices()) bad("attribute seams form a cycle");
    }

    // Points: one per (vertex, attribute-value combination); corners map to points.
    const CornerTable& ct = out.table;
    out.faces.assign(std::size_t(ct.num_corners()), 0);
    if (num_attribute_data == 0) {
        for (std::int32_t c = 0; c < ct.num_corners(); ++c) out.faces[std::size_t(c)] = std::uint32_t(ct.Vertex(c));
        out.num_points = std::uint32_t(num_connectivity_verts);
        return;
    }
    std::vector<std::int32_t> corner_to_point(std::size_t(ct.num_corners()), kInvalid);
    std::uint32_t num_points = 0;
    for (std::int32_t v = 0; v < ct.num_vertices(); ++v) {
        std::int32_t c = ct.LeftMostCorner(v);
        if (c == kInvalid) continue;
        std::int32_t dedup_first = c;
        if (!eb.is_vert_hole[std::size_t(v)]) {
            for (auto& ad : out.attribute_data) {
                if (!ad.table.IsCornerOnSeam(c)) continue;
                const std::int32_t vert_id = ad.table.Vertex(c);
                std::int32_t act = ct.SwingRight(c);
                bool seam_found = false;
                while (act != c) {
                    if (act == kInvalid) bad("open vertex fan on an interior vertex");
                    if (ad.table.Vertex(act) != vert_id) {
                        dedup_first = act;
                        seam_found = true;
                        break;
                    }
                    act = ct.SwingRight(act);
                }
                if (seam_found) break;
            }
        }
        c = dedup_first;
        corner_to_point[std::size_t(c)] = std::int32_t(num_points++);
        std::int32_t prev = c;
        c = ct.SwingRight(c);
        while (c != kInvalid && c != dedup_first) {
            bool seam = false;
            for (auto& ad : out.attribute_data)
                if (ad.table.Vertex(c) != ad.table.Vertex(prev)) { seam = true; break; }
            corner_to_point[std::size_t(c)] = seam ? std::int32_t(num_points++) : corner_to_point[std::size_t(prev)];
            prev = c;
            c = ct.SwingRight(c);
        }
    }
    for (std::size_t c = 0; c < corner_to_point.size(); ++c) {
        if (corner_to_point[c] < 0) bad("corner without a point");
        out.faces[c] = std::uint32_t(corner_to_point[c]);
    }
    out.num_points = num_points;
}

}  // namespace cadcore::codec::draco
