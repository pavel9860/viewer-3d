// meshopt bitstreams, written from KHR_meshopt_compression Appendix A/B (which includes
// the EXT_meshopt_compression v0 format).
#include "codec/meshopt.hpp"

#include <cmath>
#include <cstring>

#include "cadcore/error.hpp"

namespace cadcore::codec::meshopt {
namespace {

[[noreturn]] void bad(const std::string& w) { fail(ErrorCode::Malformed, "meshopt: " + w); }

// ---- mode 0: attributes -------------------------------------------------------------------------

// Decodes the 16 deltas of one group stored with `bits` per value (0, 1, 2, 4 or 8) and sentinel escapes.
const std::uint8_t* decode_group(const std::uint8_t* p, const std::uint8_t* end, int bits, std::uint8_t out[16]) {
    if (bits == 0) {
        std::memset(out, 0, 16);
        return p;
    }
    if (bits == 8) {
        if (end - p < 16) bad("attribute stream ends inside a group");
        std::memcpy(out, p, 16);
        return p + 16;
    }
    const int packed = 16 * bits / 8;
    if (end - p < packed) bad("attribute stream ends inside a group");
    const unsigned sentinel = (1u << bits) - 1;
    for (int i = 0; i < 16; ++i) {
        unsigned v;
        if (bits == 1) v = (p[i / 8] >> (i % 8)) & 1;                     // LSB-first
        else if (bits == 2) v = (p[i / 4] >> (6 - 2 * (i % 4))) & 3;      // MSB-first
        else v = (p[i / 2] >> (i % 2 ? 0 : 4)) & 15;                       // high nibble first
        out[i] = std::uint8_t(v);
    }
    p += packed;
    for (int i = 0; i < 16; ++i)
        if (out[i] == sentinel) {
            if (p >= end) bad("attribute stream ends inside escaped deltas");
            out[i] = *p++;
        }
    return p;
}

std::uint8_t unzigzag8(std::uint8_t v) { return std::uint8_t((v & 1) ? ~(v >> 1) : (v >> 1)); }
std::uint16_t unzigzag16(std::uint16_t v) { return std::uint16_t((v & 1) ? ~(v >> 1) : (v >> 1)); }
std::uint32_t rotl(std::uint32_t v, unsigned r) { r &= 31; return r ? (v << r) | (v >> (32 - r)) : v; }

std::vector<std::uint8_t> decode_attributes(ByteView src, std::size_t count, std::size_t stride) {
    if (stride == 0 || stride > 256 || stride % 4) bad("attribute byteStride must be a multiple of 4 in [4, 256]");
    if (src.size < 1) bad("empty attribute stream");
    const int version = src.data[0] == 0xa0 ? 0 : src.data[0] == 0xa1 ? 1 : -1;
    if (version < 0) bad("attribute stream header is not 0xa0/0xa1");
    const std::size_t tail_size = stride + (version == 1 ? stride / 4 : 0);
    const std::size_t min_tail = version == 1 ? 24 : 32;
    const std::size_t padding = tail_size < min_tail ? min_tail - tail_size : 0;
    if (src.size < 1 + tail_size + padding) bad("attribute stream shorter than its tail");
    const std::uint8_t* tail = src.data + src.size - tail_size;
    const std::uint8_t* end_of_blocks = tail - padding;
    std::uint8_t channel_mode[64] = {};
    if (version == 1) {
        for (std::size_t c = 0; c < stride / 4; ++c) {
            channel_mode[c] = tail[stride + c];
            const int m = channel_mode[c] & 15;
            if (m > 2 || (m < 2 && (channel_mode[c] >> 4))) bad("invalid channel mode");
        }
    }
    std::vector<std::uint8_t> out(count * stride);
    std::vector<std::uint8_t> prev(tail, tail + stride);  // baseline element
    const std::size_t max_block = std::min<std::size_t>((8192 / stride) & ~std::size_t(15), 256);
    if (max_block == 0) bad("byteStride too large");
    const std::uint8_t* p = src.data + 1;
    std::vector<std::uint8_t> deltas(stride * 256);
    for (std::size_t first = 0; first < count; first += max_block) {
        const std::size_t n = std::min(max_block, count - first);
        const std::size_t groups = (n + 15) / 16;
        std::uint8_t control[64] = {};
        if (version == 1) {
            if (std::size_t(end_of_blocks - p) < stride / 4) bad("attribute stream ends in a control header");
            for (std::size_t k = 0; k < stride; ++k) control[k] = std::uint8_t((p[k / 4] >> (2 * (k % 4))) & 3);
            p += stride / 4;
        }
        for (std::size_t k = 0; k < stride; ++k) {
            std::uint8_t* d = &deltas[k * 256];
            if (version == 1 && control[k] == 2) {
                std::memset(d, 0, groups * 16);
                continue;
            }
            if (version == 1 && control[k] == 3) {  // literal: one raw delta byte per element of the block
                if (std::size_t(end_of_blocks - p) < n) bad("attribute stream ends in a literal block");
                std::memcpy(d, p, n);
                std::memset(d + n, 0, groups * 16 - n);
                p += n;
                continue;
            }
            const std::size_t hdr = (groups + 3) / 4;
            if (std::size_t(end_of_blocks - p) < hdr) bad("attribute stream ends in header bits");
            const std::uint8_t* header = p;
            p += hdr;
            for (std::size_t g = 0; g < groups; ++g) {
                const int mode = (header[g / 4] >> (2 * (g % 4))) & 3;
                static constexpr int v0_bits[4] = {0, 2, 4, 8}, c0_bits[4] = {0, 1, 2, 4}, c1_bits[4] = {1, 2, 4, 8};
                const int bits = version == 0 ? v0_bits[mode] : control[k] == 0 ? c0_bits[mode] : c1_bits[mode];
                p = decode_group(p, end_of_blocks, bits, d + g * 16);
            }
        }
        // Reconstruct elements from deltas, channel by channel (4 bytes per channel).
        for (std::size_t i = 0; i < n; ++i) {
            std::uint8_t* dst = &out[(first + i) * stride];
            for (std::size_t c = 0; c < stride / 4; ++c) {
                const int mode = version == 1 ? channel_mode[c] & 15 : 0;
                const std::size_t b = c * 4;
                if (mode == 0) {
                    for (std::size_t k = b; k < b + 4; ++k) prev[k] = std::uint8_t(prev[k] + unzigzag8(deltas[k * 256 + i]));
                } else if (mode == 1) {
                    for (std::size_t k = b; k < b + 4; k += 2) {
                        const std::uint16_t dv = std::uint16_t(deltas[k * 256 + i] | deltas[(k + 1) * 256 + i] << 8);
                        const std::uint16_t pv = std::uint16_t(prev[k] | prev[k + 1] << 8);
                        const std::uint16_t v = std::uint16_t(pv + unzigzag16(dv));
                        prev[k] = std::uint8_t(v);
                        prev[k + 1] = std::uint8_t(v >> 8);
                    }
                } else {
                    std::uint32_t dv = 0, pv = 0;
                    for (int j = 0; j < 4; ++j) {
                        dv |= std::uint32_t(deltas[(b + std::size_t(j)) * 256 + i]) << (8 * j);
                        pv |= std::uint32_t(prev[b + std::size_t(j)]) << (8 * j);
                    }
                    const std::uint32_t v = pv ^ rotl(dv, 32 - (channel_mode[c] >> 4));
                    for (int j = 0; j < 4; ++j) prev[b + std::size_t(j)] = std::uint8_t(v >> (8 * j));
                }
            }
            std::memcpy(dst, prev.data(), stride);
        }
    }
    if (p != end_of_blocks) bad("unprocessed bytes before the attribute tail");
    return out;
}

// ---- mode 1: triangles / mode 2: indices --------------------------------------------------------

std::uint32_t leb128(const std::uint8_t*& p, const std::uint8_t* end) {
    std::uint32_t v = 0;
    for (int i = 0; i < 5; ++i) {
        if (p >= end) bad("index stream ends inside a varint");
        const std::uint8_t b = *p++;
        v |= std::uint32_t(b & 0x7f) << (7 * i);
        if (!(b & 0x80)) return v;
    }
    bad("varint longer than 5 bytes");
}

void store_index(std::vector<std::uint8_t>& out, std::size_t i, std::uint32_t v, std::size_t stride) {
    if (stride == 2) { const std::uint16_t s = std::uint16_t(v); std::memcpy(&out[i * 2], &s, 2); }
    else std::memcpy(&out[i * 4], &v, 4);
}

std::vector<std::uint8_t> decode_triangles(ByteView src, std::size_t count, std::size_t stride) {
    if (stride != 2 && stride != 4) bad("index byteStride must be 2 or 4");
    if (count % 3) bad("triangle index count is not a multiple of 3");
    const std::size_t tris = count / 3;
    if (src.size < 1 + tris + 16 || src.data[0] != 0xe1) bad("triangle stream header is not 0xe1 or the stream is too short");
    const std::uint8_t* code = src.data + 1;
    const std::uint8_t* data = code + tris;
    const std::uint8_t* const data_end = src.data + src.size - 16;
    const std::uint8_t* codeaux = data_end;
    std::uint32_t edge[16][2] = {}, vtx[16] = {};
    unsigned edge_off = 0, vtx_off = 0, edges_written = 0, vtx_written = 0;
    std::uint32_t next = 0, last = 0;
    auto push_edge = [&](std::uint32_t a, std::uint32_t b) { edge[edge_off & 15][0] = a; edge[edge_off & 15][1] = b; ++edge_off; ++edges_written; };
    auto push_vtx = [&](std::uint32_t v) { vtx[vtx_off & 15] = v; ++vtx_off; ++vtx_written; };
    auto get_edge = [&](unsigned x, std::uint32_t& a, std::uint32_t& b) {
        if (x >= edges_written) bad("triangle stream reads an unwritten edge FIFO entry");
        const auto* e = edge[(edge_off - 1 - x) & 15];
        a = e[0];
        b = e[1];
    };
    auto get_vtx = [&](unsigned y) -> std::uint32_t {
        if (y >= vtx_written) bad("triangle stream reads an unwritten vertex FIFO entry");
        return vtx[(vtx_off - 1 - y) & 15];
    };
    auto decode_index = [&]() {
        const std::uint32_t v = leb128(data, data_end);
        last += (v & 1) ? ~(v >> 1) : (v >> 1);
        return last;
    };
    std::vector<std::uint8_t> out(count * stride);
    for (std::size_t t = 0; t < tris; ++t) {
        const std::uint8_t c = code[t];
        std::uint32_t a, b, cc;
        const unsigned X = c >> 4, Y = c & 15;
        if (X < 15) {
            get_edge(X, a, b);
            if (Y == 0) { cc = next++; push_vtx(cc); }
            else if (Y < 13) { cc = get_vtx(Y); }
            else if (Y == 13 || Y == 14) { cc = Y == 13 ? last - 1 : last + 1; last = cc; push_vtx(cc); }
            else { cc = decode_index(); push_vtx(cc); }
            push_edge(cc, b);
            push_edge(a, cc);
        } else {
            unsigned Z, W;
            if (Y < 14) {
                const std::uint8_t aux = codeaux[Y];
                Z = aux >> 4;
                W = aux & 15;
                a = next++;
                if (Z == 0) b = next++; else b = get_vtx(Z - 1);
                if (W == 0) cc = next++; else cc = get_vtx(W - 1);
                push_vtx(a);
                if (Z == 0) push_vtx(b);
                if (W == 0) push_vtx(cc);
            } else {
                if (data >= data_end) bad("triangle stream ends early");
                const std::uint8_t zw = *data++;
                if (zw == 0) next = 0;
                Z = zw >> 4;
                W = zw & 15;
                a = Y == 14 ? next++ : decode_index();
                if (Z == 0) b = next++; else if (Z < 15) b = get_vtx(Z - 1); else b = decode_index();
                if (W == 0) cc = next++; else if (W < 15) cc = get_vtx(W - 1); else cc = decode_index();
                push_vtx(a);
                if (Z == 0 || Z == 15) push_vtx(b);
                if (W == 0 || W == 15) push_vtx(cc);
            }
            push_edge(b, a);
            push_edge(cc, b);
            push_edge(a, cc);
        }
        store_index(out, t * 3, a, stride);
        store_index(out, t * 3 + 1, b, stride);
        store_index(out, t * 3 + 2, cc, stride);
    }
    if (data != data_end) bad("unprocessed bytes before the triangle tail");
    return out;
}

std::vector<std::uint8_t> decode_indices(ByteView src, std::size_t count, std::size_t stride) {
    if (stride != 2 && stride != 4) bad("index byteStride must be 2 or 4");
    if (src.size < 1 + 4 || src.data[0] != 0xd1) bad("index stream header is not 0xd1 or the stream is too short");
    const std::uint8_t* p = src.data + 1;
    const std::uint8_t* const end = src.data + src.size - 4;
    std::uint32_t last[2] = {0, 0};
    std::vector<std::uint8_t> out(count * stride);
    for (std::size_t i = 0; i < count; ++i) {
        const std::uint32_t v = leb128(p, end);
        const std::uint32_t base = v & 1;
        last[base] += (v & 2) ? ~(v >> 2) : (v >> 2);
        store_index(out, i, last[base], stride);
    }
    if (p != end) bad("unprocessed bytes before the index tail");
    return out;
}

// ---- filters ------------------------------------------------------------------------------------

double round_away(double v) { return v < 0 ? -std::floor(-v + 0.5) : std::floor(v + 0.5); }

template <class T>
void octahedral(std::uint8_t* data, std::size_t count) {
    const double maxv = double((1 << (8 * sizeof(T) - 1)) - 1);
    for (std::size_t i = 0; i < count; ++i) {
        T in[4];
        std::memcpy(in, data + i * 4 * sizeof(T), sizeof in);
        const float one = float(in[2]);
        float x = in[0] / one, y = in[1] / one;
        const float z0 = 1.0f - std::fabs(x) - std::fabs(y);
        const float t = std::min(z0, 0.0f);
        x -= std::copysign(t, x);
        y -= std::copysign(t, y);
        const float z = z0;
        const float len = std::sqrt(x * x + y * y + z * z);
        T o[4] = {T(round_away(x / len * maxv)), T(round_away(y / len * maxv)), T(round_away(z / len * maxv)), in[3]};
        std::memcpy(data + i * 4 * sizeof(T), o, sizeof o);
    }
}

void quaternion(std::uint8_t* data, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        std::int16_t in[4];
        std::memcpy(in, data + i * 8, 8);
        const float range = 1.0f / std::sqrt(2.0f);
        const float one = float(in[3] | 3);
        const float x = in[0] / one * range, y = in[1] / one * range, z = in[2] / one * range;
        const float w = std::sqrt(std::max(0.0f, 1.0f - x * x - y * y - z * z));
        const int mc = in[3] & 3;
        std::int16_t o[4];
        o[(mc + 1) % 4] = std::int16_t(round_away(x * 32767.0));
        o[(mc + 2) % 4] = std::int16_t(round_away(y * 32767.0));
        o[(mc + 3) % 4] = std::int16_t(round_away(z * 32767.0));
        o[(mc + 0) % 4] = std::int16_t(round_away(w * 32767.0));
        std::memcpy(data + i * 8, o, 8);
    }
}

void exponential(std::uint8_t* data, std::size_t words) {
    for (std::size_t i = 0; i < words; ++i) {
        std::int32_t v;
        std::memcpy(&v, data + i * 4, 4);
        const std::int32_t e = v >> 24;
        const std::int32_t m = std::int32_t(std::uint32_t(v) << 8) >> 8;
        const float f = float(std::ldexp(double(m), e));
        std::memcpy(data + i * 4, &f, 4);
    }
}

template <class T, class S>
void color(std::uint8_t* data, std::size_t count) {
    const double umax = double(T(~T(0)));
    for (std::size_t i = 0; i < count; ++i) {
        T in[4];
        std::memcpy(in, data + i * 4 * sizeof(T), sizeof in);
        int msb = -1;
        for (int b = int(8 * sizeof(T)) - 1; b >= 0; --b)
            if (in[3] & (T(1) << b)) { msb = b; break; }
        const int as = (1 << (msb + 1)) - 1;
        const int y = in[0], co = S(in[1]), cg = S(in[2]);
        const int r = y + co - cg, g = y + cg, b = y - co - cg;
        int a = in[3] & (as >> 1);
        a = (a << 1) | (a & 1);
        const double ss = as > 0 ? umax / double(as) : 0;
        T o[4] = {T(round_away(r * ss)), T(round_away(g * ss)), T(round_away(b * ss)), T(round_away(a * ss))};
        std::memcpy(data + i * 4 * sizeof(T), o, sizeof o);
    }
}

}  // namespace

std::vector<std::uint8_t> decode(ByteView src, std::size_t count, std::size_t stride, Mode mode, Filter filter) {
    if (count > (std::size_t(1) << 28)) bad("implausible element count");
    std::vector<std::uint8_t> out;
    switch (mode) {
        case Mode::Attributes: out = decode_attributes(src, count, stride); break;
        case Mode::Triangles: out = decode_triangles(src, count, stride); break;
        case Mode::Indices: out = decode_indices(src, count, stride); break;
    }
    if (filter != Filter::None && mode != Mode::Attributes) bad("filters apply only to attribute streams");
    switch (filter) {
        case Filter::None: break;
        case Filter::Octahedral:
            if (stride == 4) octahedral<std::int8_t>(out.data(), count);
            else if (stride == 8) octahedral<std::int16_t>(out.data(), count);
            else bad("octahedral filter requires byteStride 4 or 8");
            break;
        case Filter::Quaternion:
            if (stride != 8) bad("quaternion filter requires byteStride 8");
            quaternion(out.data(), count);
            break;
        case Filter::Exponential:
            exponential(out.data(), out.size() / 4);
            break;
        case Filter::Color:
            if (stride == 4) color<std::uint8_t, std::int8_t>(out.data(), count);
            else if (stride == 8) color<std::uint16_t, std::int16_t>(out.data(), count);
            else bad("color filter requires byteStride 4 or 8");
            break;
    }
    return out;
}

}  // namespace cadcore::codec::meshopt
