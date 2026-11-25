// glTF 2.0 decoder, written from the Khronos glTF 2.0 specification.
//
// Containers: .gltf (JSON) and .glb (12-byte header "glTF" v2, JSON chunk 0x4E4F534A,
// optional BIN chunk 0x004E4942). Buffers come from the GLB BIN chunk, base64 data
// URIs, or external files through DecodeOptions::load_resource.
//
// Geometry output: every primitive becomes one Mesh (positions, indices, optional
// COLOR_0, material base-colour factor as Mesh::color). Triangle strips and fans
// are expanded; points become Topology::Points meshes; lines, line strips and
// loops become polylines. The node hierarchy of the default scene (or scene 0)
// is kept with local transforms (matrix or T*R*S).
//
// Extensions decoded: KHR_mesh_quantization (normalized / integer attributes),
// KHR_draco_mesh_compression, EXT_meshopt_compression, KHR_meshopt_compression,
// EXT_mesh_gpu_instancing. Appearance-only extensions are reported as notes;
// unknown required extensions as warnings.
//
// Validation (spec requirements that affect geometry): JSON types of every member
// read; ids in range; nodes form a strict forest (one parent, no cycles); buffer
// views inside buffers; accessors inside buffer views; index values below the
// vertex count; "extensions" members are objects of objects; "name" members are
// strings. Components are read with memcpy, so alignment never causes UB.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <string>

#include "codec/draco.hpp"
#include "core/parallel.hpp"
#include "codec/json.hpp"
#include "codec/meshopt.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

using json::Type;
using json::Value;
using MeshRefs = std::pair<std::vector<std::uint32_t>, std::vector<std::uint32_t>>;  // meshes, polylines

[[noreturn]] void bad(const std::string& what) { fail(ErrorCode::Malformed, "glTF: " + what); }

int sniff(ByteView d, std::string_view ext) {
    if (d.starts_with("glTF")) return 2;
    std::size_t i = d.starts_with("\xEF\xBB\xBF") ? 3 : 0;
    while (i < d.size && (d[i] == ' ' || d[i] == '\t' || d[i] == '\r' || d[i] == '\n')) ++i;
    const std::size_t window = std::min<std::size_t>(d.size, 1 << 16);
    if (i < d.size && d[i] == '{' && d.find("\"asset\"", i, window) != ByteView::npos && d.find("\"version\"", i, window) != ByteView::npos) return 2;
    return ext == "gltf" || ext == "glb" ? 1 : 0;
}

// ---- typed JSON access -------------------------------------------------------------------------

// `where` is passed by value (string_view) so no temporary binds to a reference parameter.
const Value& need(const Value& o, std::string_view k, Type t, std::string_view where) {
    const Value* v = o.get(k);
    if (!v) bad(std::string(where) + " lacks required \"" + std::string(k) + "\"");
    if (!v->is(t)) bad(std::string(where) + "." + std::string(k) + " has the wrong JSON type");
    return *v;
}
const Value* opt(const Value& o, std::string_view k, Type t, std::string_view where) {
    const Value* v = o.get(k);
    if (v && !v->is(t)) bad(std::string(where) + "." + std::string(k) + " has the wrong JSON type");
    return v;
}
std::size_t as_index(const Value& v, std::size_t limit, const std::string& where) {
    if (!v.is(Type::Number) || !v.is_integer || v.integer < 0) bad(where + " is not a non-negative integer");
    if (std::uint64_t(v.integer) >= limit) bad(where + " = " + std::to_string(v.integer) + " is out of range");
    return std::size_t(v.integer);
}
std::size_t uint_member(const Value& o, std::string_view k, const std::string& where, bool required, std::size_t dflt = 0) {
    const Value* v = o.get(k);
    if (!v) {
        if (required) bad(where + " lacks required \"" + std::string(k) + "\"");
        return dflt;
    }
    if (!v->is(Type::Number) || !v->is_integer || v->integer < 0) bad(where + "." + std::string(k) + " is not a non-negative integer");
    return std::size_t(v->integer);
}
const std::vector<Value>& array_of(const Value& root, std::string_view k) {
    static const std::vector<Value> empty;
    const Value* v = opt(root, k, Type::Array, "root");
    return v ? v->array : empty;
}

// Generic schema rules (glTFProperty): "extensions" is an object whose members are objects,
// and a core object's "name" is a string. Extension payloads follow their own schemas, so
// the "name" rule is not applied inside them; "extras" is free-form and skipped.
void check_schema(const Value& v, const std::string& path, int depth, bool in_extension) {
    if (depth > 64) return;
    if (v.is(Type::Object)) {
        for (const auto& [k, c] : v.object) {
            if (k == "extras") continue;
            if (k == "extensions") {
                if (!c.is(Type::Object)) bad(path + ".extensions is not an object");
                for (const auto& [ek, ev] : c.object) {
                    if (!ev.is(Type::Object)) bad(path + ".extensions." + ek + " is not an object");
                    check_schema(ev, path + ".extensions." + ek, depth + 1, true);
                }
                continue;
            }
            if (!in_extension && k == "name" && !c.is(Type::String)) bad(path + ".name is not a string");
            check_schema(c, path + "." + k, depth + 1, in_extension);
        }
    } else if (v.is(Type::Array)) {
        for (std::size_t i = 0; i < v.array.size(); ++i) check_schema(v.array[i], path + "[" + std::to_string(i) + "]", depth + 1, in_extension);
    }
}

std::size_t component_size(std::int64_t ct) {
    switch (ct) {
        case 5120: case 5121: return 1;
        case 5122: case 5123: return 2;
        case 5125: case 5126: return 4;
        default: return 0;
    }
}
int component_count(const std::string& t) {
    if (t == "SCALAR") return 1;
    if (t == "VEC2") return 2;
    if (t == "VEC3") return 3;
    if (t == "VEC4") return 4;
    if (t == "MAT2") return 4;
    if (t == "MAT3") return 9;
    if (t == "MAT4") return 16;
    return 0;
}

// Every core cross-reference must name an existing object, and every accessor must use a
// glTF 2.0 component type and accessor type, whether or not the geometry uses it.
void validate_references(const Value& root) {
    auto count = [&](std::string_view k) { return array_of(root, k).size(); };
    auto ref = [&](const Value& o, std::string_view k, std::string_view target, const std::string& w) {
        if (const Value* v = o.get(k)) as_index(*v, count(target), w + "." + std::string(k));
    };
    auto each = [&](std::string_view k, auto&& fn) {
        const auto& arr = array_of(root, k);
        for (std::size_t i = 0; i < arr.size(); ++i) {
            const std::string w = std::string(k) + "[" + std::to_string(i) + "]";
            if (!arr[i].is(Type::Object)) bad(w + " is not an object");
            fn(arr[i], w);
        }
    };
    auto texture_info = [&](const Value* t, const std::string& w) {
        if (!t) return;
        if (!t->is(Type::Object)) bad(w + " is not an object");
        ref(*t, "index", "textures", w);
    };
    each("accessors", [&](const Value& a, const std::string& w) {
        ref(a, "bufferView", "bufferViews", w);
        const Value& ct = need(a, "componentType", Type::Number, w);
        if (!ct.is_integer || !component_size(ct.integer)) bad(w + ".componentType " + std::to_string(std::int64_t(ct.number)) + " is not allowed in glTF 2.0");
        if (!component_count(need(a, "type", Type::String, w).string)) bad(w + ".type is not a glTF accessor type");
        uint_member(a, "count", w, true);
        if (const Value* sp = opt(a, "sparse", Type::Object, w)) {
            ref(need(*sp, "indices", Type::Object, w + ".sparse"), "bufferView", "bufferViews", w + ".sparse.indices");
            ref(need(*sp, "values", Type::Object, w + ".sparse"), "bufferView", "bufferViews", w + ".sparse.values");
        }
    });
    each("bufferViews", [&](const Value& v, const std::string& w) { ref(v, "buffer", "buffers", w); });
    each("images", [&](const Value& v, const std::string& w) { ref(v, "bufferView", "bufferViews", w); });
    each("textures", [&](const Value& v, const std::string& w) {
        ref(v, "sampler", "samplers", w);
        ref(v, "source", "images", w);
    });
    each("materials", [&](const Value& m, const std::string& w) {
        if (const Value* pbr = opt(m, "pbrMetallicRoughness", Type::Object, w)) {
            texture_info(pbr->get("baseColorTexture"), w + ".baseColorTexture");
            texture_info(pbr->get("metallicRoughnessTexture"), w + ".metallicRoughnessTexture");
        }
        texture_info(m.get("normalTexture"), w + ".normalTexture");
        texture_info(m.get("occlusionTexture"), w + ".occlusionTexture");
        texture_info(m.get("emissiveTexture"), w + ".emissiveTexture");
    });
    each("meshes", [&](const Value& m, const std::string& w) {
        for (const Value& p : need(m, "primitives", Type::Array, w).array) {
            if (!p.is(Type::Object)) bad(w + " primitive is not an object");
            ref(p, "indices", "accessors", w);
            ref(p, "material", "materials", w);
            for (const auto& [k, v] : need(p, "attributes", Type::Object, w).object) as_index(v, count("accessors"), w + ".attributes." + k);
            if (const Value* targets = opt(p, "targets", Type::Array, w))
                for (const Value& t : targets->array)
                    for (const auto& [k, v] : t.object) as_index(v, count("accessors"), w + ".targets." + k);
        }
    });
    each("nodes", [&](const Value& n, const std::string& w) {
        ref(n, "mesh", "meshes", w);
        ref(n, "camera", "cameras", w);
        ref(n, "skin", "skins", w);
    });
    each("skins", [&](const Value& sk, const std::string& w) {
        ref(sk, "inverseBindMatrices", "accessors", w);
        ref(sk, "skeleton", "nodes", w);
        for (const Value& j : need(sk, "joints", Type::Array, w).array) as_index(j, count("nodes"), w + ".joints");
    });
    each("animations", [&](const Value& an, const std::string& w) {
        const std::size_t nsamplers = need(an, "samplers", Type::Array, w).array.size();
        for (const Value& smp : an.get("samplers")->array) {
            if (!smp.is(Type::Object)) bad(w + " sampler is not an object");
            ref(smp, "input", "accessors", w);
            ref(smp, "output", "accessors", w);
        }
        for (const Value& ch : need(an, "channels", Type::Array, w).array) {
            if (!ch.is(Type::Object)) bad(w + " channel is not an object");
            as_index(need(ch, "sampler", Type::Number, w), nsamplers, w + ".channel.sampler");
            if (const Value* t = opt(ch, "target", Type::Object, w)) ref(*t, "node", "nodes", w + ".target");
        }
    });
    each("scenes", [&](const Value& sc, const std::string& w) {
        if (const Value* ns = opt(sc, "nodes", Type::Array, w))
            for (const Value& n : ns->array) as_index(n, count("nodes"), w + ".nodes");
    });
    if (const Value* sv = root.get("scene")) as_index(*sv, count("scenes"), "root.scene");
}

// ---- data URIs ----------------------------------------------------------------------------------

std::vector<std::uint8_t> base64(std::string_view s) {
    auto val = [](char c) -> int {
        if (c >= 'A' && c <= 'Z') return c - 'A';
        if (c >= 'a' && c <= 'z') return c - 'a' + 26;
        if (c >= '0' && c <= '9') return c - '0' + 52;
        if (c == '+' || c == '-') return 62;
        if (c == '/' || c == '_') return 63;
        return -1;
    };
    std::vector<std::uint8_t> out;
    out.reserve(s.size() * 3 / 4);
    std::uint32_t acc = 0;
    int bits = 0;
    for (char c : s) {
        if (c == '=') break;
        if (c == '\r' || c == '\n' || c == ' ') continue;
        const int v = val(c);
        if (v < 0) bad("invalid base64 in data URI");
        acc = acc << 6 | std::uint32_t(v);
        bits += 6;
        if (bits >= 8) { bits -= 8; out.push_back(std::uint8_t(acc >> bits)); }
    }
    return out;
}

std::string percent_decode(std::string_view s) {
    auto h = [](char c) { return c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1; };
    std::string out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '%' && i + 2 < s.size() && h(s[i + 1]) >= 0 && h(s[i + 2]) >= 0) {
            out += char(h(s[i + 1]) * 16 + h(s[i + 2]));
            i += 2;
        } else {
            out += s[i];
        }
    }
    return out;
}

// ---- buffers, views, accessors ------------------------------------------------------------------

struct BufferView {
    std::size_t buffer = 0, offset = 0, length = 0, stride = 0;
    std::vector<std::uint8_t> decoded;  // meshopt-decompressed contents
    bool has_decoded = false;
};

class Document {
public:
    Document(const Value& root, std::vector<std::vector<std::uint8_t>> buffers) : root_(root), buffers_(std::move(buffers)) {
        const auto& views = array_of(root, "bufferViews");
        views_.resize(views.size());
        for (std::size_t i = 0; i < views.size(); ++i) {
            const std::string w = "bufferViews[" + std::to_string(i) + "]";
            if (!views[i].is(Type::Object)) bad(w + " is not an object");
            BufferView& bv = views_[i];
            bv.buffer = as_index(need(views[i], "buffer", Type::Number, w), buffers_.size(), w + ".buffer");
            bv.offset = uint_member(views[i], "byteOffset", w, false);
            bv.length = uint_member(views[i], "byteLength", w, true);
            bv.stride = uint_member(views[i], "byteStride", w, false);
            if (bv.offset > buffers_[bv.buffer].size() || bv.length > buffers_[bv.buffer].size() - bv.offset)
                bad(w + " extends past the end of buffer " + std::to_string(bv.buffer));
        }
    }

    const Value& root() const { return root_; }
    std::vector<BufferView>& views() { return views_; }

    // Draco streams decoded ahead (in parallel) by buffer view; a failure is kept and raised
    // when its primitive is converted, so errors appear in document order as without threads.
    struct Draco { codec::draco::Mesh mesh; std::exception_ptr error; };
    std::map<std::size_t, Draco> draco;
    const codec::draco::Mesh& draco_mesh(std::size_t view) const {
        if (auto it = draco.find(view); it != draco.end()) {
            if (it->second.error) std::rethrow_exception(it->second.error);
            return it->second.mesh;
        }
        fail(ErrorCode::Malformed, "glTF: Draco stream of bufferView " + std::to_string(view) + " was not decoded");
    }
    const std::vector<std::vector<std::uint8_t>>& buffers() const { return buffers_; }

    ByteView view_bytes(std::size_t i) const {
        const BufferView& bv = views_.at(i);
        if (bv.has_decoded) return ByteView(bv.decoded.data(), bv.decoded.size());
        return ByteView(buffers_[bv.buffer].data() + bv.offset, bv.length);
    }

    // Accessor values as doubles (exact for every glTF component type), sparse substitution applied.
    std::vector<double> values(std::size_t ai, int* ncomp_out) const {
        const auto& accs = array_of(root_, "accessors");
        const std::string w = "accessors[" + std::to_string(ai) + "]";
        const Value& a = accs.at(ai);
        if (!a.is(Type::Object)) bad(w + " is not an object");
        const std::int64_t ct = need(a, "componentType", Type::Number, w).integer;
        const std::size_t cs = component_size(ct);
        if (!cs) bad(w + ".componentType " + std::to_string(ct) + " is not allowed in glTF 2.0");
        const int nc = component_count(need(a, "type", Type::String, w).string);
        if (!nc) bad(w + ".type is not a glTF accessor type");
        const std::size_t count = uint_member(a, "count", w, true);
        if (count > (std::size_t(1) << 31)) bad(w + ".count is implausible");
        const Value* norm = opt(a, "normalized", Type::Bool, w);
        const bool normalized = norm && norm->boolean;
        std::vector<double> out(count * std::size_t(nc), 0.0);
        auto convert = [&](const std::uint8_t* p) -> double {
            switch (ct) {
                case 5120: { std::int8_t v; std::memcpy(&v, p, 1); return normalized ? std::max(v / 127.0, -1.0) : double(v); }
                case 5121: return normalized ? p[0] / 255.0 : double(p[0]);
                case 5122: { std::int16_t v; std::memcpy(&v, p, 2); return normalized ? std::max(v / 32767.0, -1.0) : double(v); }
                case 5123: { std::uint16_t v; std::memcpy(&v, p, 2); return normalized ? v / 65535.0 : double(v); }
                case 5125: { std::uint32_t v; std::memcpy(&v, p, 4); return double(v); }
                default: { float v; std::memcpy(&v, p, 4); return double(v); }
            }
        };
        if (const Value* bvv = a.get("bufferView")) {
            const std::size_t vi = as_index(*bvv, views_.size(), w + ".bufferView");
            const ByteView data = view_bytes(vi);
            const std::size_t elem = cs * std::size_t(nc);
            const std::size_t stride = views_[vi].stride ? views_[vi].stride : elem;
            const std::size_t off = uint_member(a, "byteOffset", w, false);
            if (count && (off > data.size || stride * (count - 1) + elem > data.size - off)) bad(w + " reads past the end of its buffer view");
            for (std::size_t i = 0; i < count; ++i)
                for (int k = 0; k < nc; ++k) out[i * std::size_t(nc) + std::size_t(k)] = convert(data.data + off + i * stride + std::size_t(k) * cs);
        }
        if (const Value* sp = opt(a, "sparse", Type::Object, w)) {
            const std::string ws = w + ".sparse";
            const std::size_t scount = uint_member(*sp, "count", ws, true);
            const Value& idx = need(*sp, "indices", Type::Object, ws);
            const Value& vals = need(*sp, "values", Type::Object, ws);
            const std::size_t iv = as_index(need(idx, "bufferView", Type::Number, ws + ".indices"), views_.size(), ws + ".indices.bufferView");
            const std::size_t vv = as_index(need(vals, "bufferView", Type::Number, ws + ".values"), views_.size(), ws + ".values.bufferView");
            const std::int64_t ict = need(idx, "componentType", Type::Number, ws + ".indices").integer;
            const std::size_t ics = (ict == 5121 || ict == 5123 || ict == 5125) ? component_size(ict) : 0;
            if (!ics) bad(ws + ".indices.componentType is not an unsigned integer type");
            const ByteView id = view_bytes(iv), vd = view_bytes(vv);
            const std::size_t ioff = uint_member(idx, "byteOffset", ws, false), voff = uint_member(vals, "byteOffset", ws, false);
            if (ioff > id.size || scount * ics > id.size - ioff || voff > vd.size || scount * cs * std::size_t(nc) > vd.size - voff) bad(ws + " reads past its buffer views");
            for (std::size_t s = 0; s < scount; ++s) {
                std::uint32_t target = 0;
                const std::uint8_t* ip = id.data + ioff + s * ics;
                if (ics == 1) target = ip[0];
                else if (ics == 2) { std::uint16_t t; std::memcpy(&t, ip, 2); target = t; }
                else std::memcpy(&target, ip, 4);
                if (target >= count) bad(ws + " index out of range");
                for (int k = 0; k < nc; ++k) out[target * std::size_t(nc) + std::size_t(k)] = convert(vd.data + voff + (s * std::size_t(nc) + std::size_t(k)) * cs);
            }
        }
        *ncomp_out = nc;
        return out;
    }

    std::vector<float> floats(std::size_t ai, int* nc) const {
        const std::vector<double> d = values(ai, nc);
        return std::vector<float>(d.begin(), d.end());
    }

    // floats() decoded once per accessor: primitives commonly share one vertex buffer.
    const std::vector<float>& cached_floats(std::size_t ai, int* nc) const {
        auto it = float_cache_.find(ai);
        if (it == float_cache_.end()) {
            int n = 0;
            std::vector<float> v = floats(ai, &n);
            it = float_cache_.emplace(ai, std::make_pair(std::move(v), n)).first;
        }
        *nc = it->second.second;
        return it->second.first;
    }

    std::vector<std::uint32_t> indices(std::size_t ai, std::size_t vertex_count) const {
        const std::string w = "accessors[" + std::to_string(ai) + "]";
        const std::int64_t ct = need(array_of(root_, "accessors").at(ai), "componentType", Type::Number, w).integer;
        if (ct != 5121 && ct != 5123 && ct != 5125) bad(w + " is used as indices but is not an unsigned integer type");
        int nc = 0;
        const std::vector<double> d = values(ai, &nc);
        if (nc != 1) bad(w + " is used as indices but is not SCALAR");
        std::vector<std::uint32_t> out(d.size());
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (d[i] >= double(vertex_count)) bad("index value " + std::to_string(std::uint64_t(d[i])) + " >= vertex count " + std::to_string(vertex_count) + " (" + w + ")");
            out[i] = std::uint32_t(d[i]);
        }
        return out;
    }

private:
    mutable std::map<std::size_t, std::pair<std::vector<float>, int>> float_cache_;
    const Value& root_;
    std::vector<std::vector<std::uint8_t>> buffers_;
    std::vector<BufferView> views_;
};

std::uint8_t unit_to_u8(double v) { return std::uint8_t(v <= 0 ? 0 : v >= 1 ? 255 : v * 255.0 + 0.5); }

// ---- extensions ---------------------------------------------------------------------------------

bool geometry_supported(std::string_view e) {
    return e == "KHR_mesh_quantization" || e == "KHR_draco_mesh_compression" || e == "EXT_meshopt_compression" || e == "KHR_meshopt_compression" ||
           e == "EXT_mesh_gpu_instancing";
}
bool appearance_only(std::string_view e) {
    return e.rfind("KHR_materials_", 0) == 0 || e.rfind("KHR_texture_", 0) == 0 || e.rfind("EXT_texture_", 0) == 0 || e == "KHR_lights_punctual" ||
           e == "KHR_technique_webgl" || e == "KHR_xmp_json_ld" || e == "KHR_xmp" || e == "KHR_animation_pointer" || e == "EXT_lights_image_based" ||
           e == "MSFT_texture_dds" || e == "KHR_interactivity" || e == "KHR_node_visibility";
}

void classify_extensions(const Value& root, Scene& scene) {
    std::set<std::string> required;
    for (const Value& v : array_of(root, "extensionsRequired")) {
        if (!v.is(Type::String)) bad("extensionsRequired entries must be strings");
        required.insert(v.string);
    }
    for (const Value& v : array_of(root, "extensionsUsed")) {
        if (!v.is(Type::String)) bad("extensionsUsed entries must be strings");
        const std::string& e = v.string;
        if (geometry_supported(e)) continue;
        if (appearance_only(e)) scene.notes.push_back("appearance extension " + e + " not rendered (geometry unaffected)");
        else if (required.count(e)) scene.warnings.push_back("required extension " + e + " is not supported; geometry may be incomplete");
        else scene.notes.push_back("optional extension " + e + " ignored");
    }
}

void decode_meshopt_views(Document& doc) {
    const auto& views = array_of(doc.root(), "bufferViews");
    for (std::size_t i = 0; i < views.size(); ++i) {
        const Value* ext = views[i].get("extensions");
        if (!ext) continue;
        const Value* mo = ext->get("EXT_meshopt_compression");
        if (!mo) mo = ext->get("KHR_meshopt_compression");
        if (!mo) continue;
        const std::string w = "bufferViews[" + std::to_string(i) + "].meshopt";
        const std::size_t buf = as_index(need(*mo, "buffer", Type::Number, w), doc.buffers().size(), w + ".buffer");
        const std::size_t off = uint_member(*mo, "byteOffset", w, false), len = uint_member(*mo, "byteLength", w, true);
        const std::size_t stride = uint_member(*mo, "byteStride", w, true), count = uint_member(*mo, "count", w, true);
        const auto& src = doc.buffers()[buf];
        if (off > src.size() || len > src.size() - off) bad(w + " source range outside buffer");
        using codec::meshopt::Filter;
        using codec::meshopt::Mode;
        const std::string& ms = need(*mo, "mode", Type::String, w).string;
        Mode mode = Mode::Attributes;
        if (ms == "TRIANGLES") mode = Mode::Triangles;
        else if (ms == "INDICES") mode = Mode::Indices;
        else if (ms != "ATTRIBUTES") bad(w + ".mode '" + ms + "' is unknown");
        Filter filter = Filter::None;
        if (const Value* fv = opt(*mo, "filter", Type::String, w)) {
            const std::string& fs = fv->string;
            if (fs == "OCTAHEDRAL") filter = Filter::Octahedral;
            else if (fs == "QUATERNION") filter = Filter::Quaternion;
            else if (fs == "EXPONENTIAL") filter = Filter::Exponential;
            else if (fs == "COLOR") filter = Filter::Color;
            else if (fs != "NONE") bad(w + ".filter '" + fs + "' is unknown");
        }
        BufferView& bv = doc.views()[i];
        bv.decoded = codec::meshopt::decode(ByteView(src.data() + off, len), count, stride, mode, filter);
        bv.has_decoded = true;
        if (!bv.stride) bv.stride = stride;
    }
}

// ---- primitives ---------------------------------------------------------------------------------

struct Geometry {
    std::vector<float> positions, colors;  // colors: color_components per vertex, may be empty
    int color_components = 0;
    std::vector<std::uint32_t> indices;
    bool indexed = false;
};

Geometry draco_geometry(const Document& doc, const Value& ext, const std::string& w) {
    const std::size_t vi = as_index(need(ext, "bufferView", Type::Number, w), array_of(doc.root(), "bufferViews").size(), w + ".bufferView");
    const codec::draco::Mesh& dm = doc.draco_mesh(vi);
    const Value& attrs = need(ext, "attributes", Type::Object, w);
    const Value* pos = attrs.get("POSITION");
    if (!pos || !pos->is_integer) bad(w + " lacks a POSITION attribute id");
    const codec::draco::Attribute* dpos = dm.find(std::uint32_t(pos->integer));
    if (!dpos || dpos->components != 3) bad(w + " POSITION id is not a 3-component attribute of the Draco stream");
    Geometry g;
    g.positions = dpos->values;
    if (const Value* col = attrs.get("COLOR_0"); col && col->is_integer) {
        const codec::draco::Attribute* c = dm.find(std::uint32_t(col->integer));
        if (!c) bad(w + " COLOR_0 id not present in the Draco stream");
        g.color_components = c->components;
        g.colors = c->values;
    }
    g.indices = dm.faces;
    g.indexed = true;
    return g;
}

// Triangle primitives of one glTF mesh that share their vertex accessors become one Mesh
// with a face group per primitive (the vertex data is copied once, not once per primitive).
struct Batch {
    std::size_t mesh = std::size_t(-1);       // scene mesh being extended
    std::size_t pos = 0, col = std::size_t(-1);
};

void convert_primitive(const Document& doc, const Value& prim, const std::string& name, const std::string& w, Scene& scene, MeshRefs& refs, Batch& batch) {
    const Value& attrs = need(prim, "attributes", Type::Object, w);
    const std::size_t nacc = array_of(doc.root(), "accessors").size();
    const std::size_t mode = uint_member(prim, "mode", w, false, 4);
    if (mode > 6) bad(w + ".mode " + std::to_string(mode) + " is not a glTF primitive mode");
    Geometry g;
    std::size_t pos_ai = std::size_t(-1), col_ai = std::size_t(-1);
    const std::vector<float>* shared_pos = nullptr;  // accessor data owned by the cache
    const std::vector<float>* shared_col = nullptr;
    const Value* exts = opt(prim, "extensions", Type::Object, w);
    if (const Value* draco = exts ? exts->get("KHR_draco_mesh_compression") : nullptr) {
        g = draco_geometry(doc, *draco, w + ".KHR_draco_mesh_compression");
    } else {
        const Value* pos = attrs.get("POSITION");
        if (!pos) bad(w + " has no POSITION attribute");
        int nc = 0;
        pos_ai = as_index(*pos, nacc, w + ".attributes.POSITION");
        shared_pos = &doc.cached_floats(pos_ai, &nc);
        if (nc != 3) bad(w + " POSITION is not VEC3");
        if (const Value* col = attrs.get("COLOR_0")) {
            col_ai = as_index(*col, nacc, w + ".attributes.COLOR_0");
            shared_col = &doc.cached_floats(col_ai, &g.color_components);
            if (g.color_components != 3 && g.color_components != 4) bad(w + " COLOR_0 is not VEC3/VEC4");
        }
        if (const Value* ind = prim.get("indices")) {
            g.indices = doc.indices(as_index(*ind, nacc, w + ".indices"), shared_pos->size() / 3);
            g.indexed = true;
        }
    }
    const std::vector<float>& P = shared_pos ? *shared_pos : g.positions;
    const std::vector<float>& C = shared_col ? *shared_col : g.colors;
    const std::size_t nv = P.size() / 3;
    if (!g.indexed) {
        g.indices.resize(nv);
        for (std::uint32_t i = 0; i < nv; ++i) g.indices[i] = i;
    }
    for (auto i : g.indices)
        if (i >= nv) bad(w + " index out of range");
    Color base{255, 255, 255, 255};
    if (const Value* mat = prim.get("material")) {
        const auto& mats = array_of(doc.root(), "materials");
        const Value& m = mats[as_index(*mat, mats.size(), w + ".material")];
        if (!m.is(Type::Object)) bad("material is not an object");
        if (const Value* pbr = opt(m, "pbrMetallicRoughness", Type::Object, "material")) {
            if (const Value* f = opt(*pbr, "baseColorFactor", Type::Array, "pbrMetallicRoughness")) {
                if (f->array.size() != 4) bad("baseColorFactor must have 4 numbers");
                double c[4];
                for (int k = 0; k < 4; ++k) {
                    if (!f->array[std::size_t(k)].is(Type::Number)) bad("baseColorFactor must have 4 numbers");
                    c[k] = f->array[std::size_t(k)].number;
                }
                base = {unit_to_u8(c[0]), unit_to_u8(c[1]), unit_to_u8(c[2]), unit_to_u8(c[3])};
            }
        }
    }
    const auto& idx = g.indices;
    if (mode >= 1 && mode <= 3) {  // LINES, LINE_LOOP, LINE_STRIP
        auto start_line = [&]() -> Polyline& {
            refs.second.push_back(std::uint32_t(scene.polylines.size()));
            scene.polylines.emplace_back().color = base;
            return scene.polylines.back();
        };
        auto put = [&](Polyline& pl, std::uint32_t i) { pl.points.insert(pl.points.end(), &P[i * 3], &P[i * 3] + 3); };
        if (mode == 1) {
            for (std::size_t i = 0; i + 1 < idx.size(); i += 2) {
                Polyline& pl = start_line();
                put(pl, idx[i]);
                put(pl, idx[i + 1]);
            }
        } else if (!idx.empty()) {
            Polyline& pl = start_line();
            for (auto i : idx) put(pl, i);
            if (mode == 2) put(pl, idx[0]);
        }
        return;
    }
    // Triangles from the batched vertex buffer.
    std::vector<std::uint32_t> tri;
    if (mode >= 4) {
        if (mode == 4) {
            tri.assign(idx.begin(), idx.end() - std::ptrdiff_t(idx.size() % 3));
        } else if (mode == 5) {  // TRIANGLE_STRIP
            for (std::size_t i = 0; i + 2 < idx.size(); ++i) {
                if (i % 2 == 0) tri.insert(tri.end(), {idx[i], idx[i + 1], idx[i + 2]});
                else tri.insert(tri.end(), {idx[i + 1], idx[i], idx[i + 2]});
            }
        } else {  // TRIANGLE_FAN
            for (std::size_t i = 1; i + 1 < idx.size(); ++i) tri.insert(tri.end(), {idx[0], idx[i], idx[i + 1]});
        }
        if (shared_pos && batch.mesh != std::size_t(-1) && batch.pos == pos_ai && batch.col == col_ai) {
            Mesh& bm = scene.meshes[batch.mesh];
            if (bm.groups.empty()) bm.groups.push_back({0, std::uint32_t(bm.triangle_count()), bm.color});
            const auto first = std::uint32_t(bm.triangle_count());
            bm.indices.insert(bm.indices.end(), tri.begin(), tri.end());
            bm.groups.push_back({first, std::uint32_t(bm.triangle_count()) - first, base});
            return;
        }
    }
    Mesh m;
    m.name = name;
    m.color = base;
    std::vector<std::uint8_t> vcol;
    if (!C.empty()) {
        const std::size_t nc = std::size_t(g.color_components);
        if (C.size() / nc != nv) bad(w + " COLOR_0 count differs from POSITION count");
        vcol.resize(nv * 4);
        for (std::size_t i = 0; i < nv; ++i)
            for (std::size_t k = 0; k < 4; ++k) vcol[i * 4 + k] = k < nc ? unit_to_u8(C[i * nc + k]) : 255;
    }
    if (mode == 0) {  // POINTS
        m.topology = Topology::Points;
        for (auto i : idx) {
            m.positions.insert(m.positions.end(), &P[i * 3], &P[i * 3] + 3);
            if (!vcol.empty()) m.colors.insert(m.colors.end(), &vcol[i * 4], &vcol[i * 4] + 4);
        }
    } else {
        m.positions = shared_pos ? *shared_pos : std::move(g.positions);
        m.colors = std::move(vcol);
        m.indices = std::move(tri);
    }
    if (mode >= 4 && shared_pos) {
        batch.mesh = scene.meshes.size();
        batch.pos = pos_ai;
        batch.col = col_ai;
    } else {
        batch.mesh = std::size_t(-1);
    }
    refs.first.push_back(std::uint32_t(scene.meshes.size()));
    scene.meshes.push_back(std::move(m));
}

// Keeps only the vertices a mesh references when it uses a small part of a shared buffer
// (many meshes over one accessor would otherwise copy the whole buffer each).
void compact(Mesh& m) {
    if (m.topology != Topology::Triangles || m.indices.size() * 4 >= m.vertex_count()) return;
    std::unordered_map<std::uint32_t, std::uint32_t> remap;
    remap.reserve(m.indices.size());
    std::vector<float> pos;
    std::vector<std::uint8_t> col;
    for (auto& i : m.indices) {
        auto [it, fresh] = remap.emplace(i, std::uint32_t(remap.size()));
        if (fresh) {
            pos.insert(pos.end(), &m.positions[std::size_t(i) * 3], &m.positions[std::size_t(i) * 3] + 3);
            if (!m.colors.empty()) col.insert(col.end(), &m.colors[std::size_t(i) * 4], &m.colors[std::size_t(i) * 4] + 4);
        }
        i = it->second;
    }
    m.positions = std::move(pos);
    m.colors = std::move(col);
}

// ---- transforms ---------------------------------------------------------------------------------

Matrix4 trs(const double t[3], const double q[4], const double s[3]) {
    const double x = q[0], y = q[1], z = q[2], w = q[3];
    return {(1 - 2 * (y * y + z * z)) * s[0], 2 * (x * y + z * w) * s[0], 2 * (x * z - y * w) * s[0], 0,
            2 * (x * y - z * w) * s[1], (1 - 2 * (x * x + z * z)) * s[1], 2 * (y * z + x * w) * s[1], 0,
            2 * (x * z + y * w) * s[2], 2 * (y * z - x * w) * s[2], (1 - 2 * (x * x + y * y)) * s[2], 0,
            t[0], t[1], t[2], 1};
}

void read_numbers(const Value* v, double* out, std::size_t n, const std::string& w) {
    if (!v) return;
    if (!v->is(Type::Array) || v->array.size() != n) bad(w + " must be an array of " + std::to_string(n) + " numbers");
    for (std::size_t i = 0; i < n; ++i) {
        if (!v->array[i].is(Type::Number)) bad(w + " must contain numbers");
        out[i] = v->array[i].number;
    }
}

Matrix4 node_transform(const Value& n, const std::string& w) {
    if (const Value* m = n.get("matrix")) {
        Matrix4 out;
        read_numbers(m, out.data(), 16, w + ".matrix");
        return out;
    }
    double t[3] = {0, 0, 0}, q[4] = {0, 0, 0, 1}, s[3] = {1, 1, 1};
    read_numbers(n.get("translation"), t, 3, w + ".translation");
    read_numbers(n.get("rotation"), q, 4, w + ".rotation");
    read_numbers(n.get("scale"), s, 3, w + ".scale");
    return trs(t, q, s);
}

// ---- top level ----------------------------------------------------------------------------------

std::vector<std::vector<std::uint8_t>> load_buffers(const Value& root, std::vector<std::uint8_t>& bin, bool have_bin, const DecodeOptions& o) {
    std::vector<std::vector<std::uint8_t>> buffers;
    const auto& bufs = array_of(root, "buffers");
    for (std::size_t i = 0; i < bufs.size(); ++i) {
        const std::string w = "buffers[" + std::to_string(i) + "]";
        if (!bufs[i].is(Type::Object)) bad(w + " is not an object");
        const std::size_t len = uint_member(bufs[i], "byteLength", w, true);
        std::vector<std::uint8_t> data;
        if (const Value* uri = opt(bufs[i], "uri", Type::String, w)) {
            const std::string& u = uri->string;
            if (u.rfind("data:", 0) == 0) {
                const auto comma = u.find(',');
                if (comma == std::string::npos || u.substr(0, comma).find(";base64") == std::string::npos) bad(w + " data URI is not base64");
                data = base64(std::string_view(u).substr(comma + 1));
            } else {
                if (!o.load_resource) fail(ErrorCode::MissingResource, "external buffer '" + u + "' (no resource loader)");
                auto bytes = o.load_resource(percent_decode(u));
                if (!bytes) fail(ErrorCode::MissingResource, "external buffer '" + u + "' not found");
                data = std::move(*bytes);
            }
        } else if (i == 0 && have_bin) {
            data = std::move(bin);
        } else {
            // EXT/KHR_meshopt_compression fallback buffer: only referenced by compressed views, never loaded.
            const Value* ex = bufs[i].get("extensions");
            const Value* mo = ex ? (ex->get("EXT_meshopt_compression") ? ex->get("EXT_meshopt_compression") : ex->get("KHR_meshopt_compression")) : nullptr;
            const Value* fb = mo ? mo->get("fallback") : nullptr;
            if (!(fb && fb->is(Type::Bool) && fb->boolean)) bad(w + " has no uri and is not the GLB binary chunk");
            data.assign(len, 0);
        }
        if (data.size() < len) bad(w + " holds " + std::to_string(data.size()) + " bytes but byteLength is " + std::to_string(len));
        data.resize(len);
        buffers.push_back(std::move(data));
    }
    return buffers;
}

Scene decode(ByteView d, const DecodeOptions& o) {
    Scene scene;
    std::string_view json_text;
    std::vector<std::uint8_t> bin;
    bool have_bin = false;
    if (d.starts_with("glTF")) {
        scene.format = "glb";
        Reader r(d, 4);
        const std::uint32_t version = r.u32le(), length = r.u32le();
        if (version != 2) fail(ErrorCode::Unsupported, "GLB version " + std::to_string(version) + " (only glTF 2.0)");
        if (length > d.size) fail(ErrorCode::Truncated, "GLB header length exceeds the file size");
        const std::uint32_t jlen = r.u32le(), jtype = r.u32le();
        if (jtype != 0x4E4F534A) bad("first GLB chunk is not JSON");
        json_text = std::string_view(reinterpret_cast<const char*>(r.bytes(jlen)), jlen);
        if (r.pos() + 8 <= length) {
            const std::uint32_t blen = r.u32le(), btype = r.u32le();
            if (btype == 0x004E4942) {
                const std::uint8_t* p = r.bytes(blen);
                bin.assign(p, p + blen);
                have_bin = true;
            }
        }
    } else {
        scene.format = "gltf";
        json_text = d.text();
        if (json_text.substr(0, 3) == "\xEF\xBB\xBF") json_text.remove_prefix(3);
    }
    const Value root = json::parse(json_text);
    if (!root.is(Type::Object)) bad("top level is not an object");
    check_schema(root, "root", 0, false);
    validate_references(root);
    const Value& asset = need(root, "asset", Type::Object, "root");
    const std::string& ver = need(asset, "version", Type::String, "asset").string;
    if (ver.rfind("2.", 0) != 0) fail(ErrorCode::Unsupported, "glTF version " + ver + " (only 2.x)");
    classify_extensions(root, scene);
    Document doc(root, load_buffers(root, bin, have_bin, o));
    decode_meshopt_views(doc);
    // every KHR_draco_mesh_compression stream (by buffer view), decoded in parallel
    {
        const std::size_t nviews = array_of(root, "bufferViews").size();
        for (const Value& m : array_of(root, "meshes")) {
            const Value* prims = m.is(Type::Object) ? m.get("primitives") : nullptr;
            if (!prims || !prims->is(Type::Array)) continue;
            for (const Value& p : prims->array) {
                const Value* exts = p.is(Type::Object) ? p.get("extensions") : nullptr;
                const Value* dr = exts && exts->is(Type::Object) ? exts->get("KHR_draco_mesh_compression") : nullptr;
                const Value* bv = dr && dr->is(Type::Object) ? dr->get("bufferView") : nullptr;
                if (bv && bv->is_integer && bv->integer >= 0 && std::size_t(bv->integer) < nviews) doc.draco.emplace(std::size_t(bv->integer), Document::Draco{});
            }
        }
        std::vector<std::function<void()>> tasks;
        for (auto& [view, slot] : doc.draco)
            tasks.push_back([&doc, view = view, &slot = slot] {
                try {
                    slot.mesh = codec::draco::decode(doc.view_bytes(view));
                } catch (...) {
                    slot.error = std::current_exception();
                }
            });
        run_all(tasks, o.threads);
    }

    const auto& meshes = array_of(root, "meshes");
    std::map<std::size_t, MeshRefs> converted;
    auto mesh_of = [&](std::size_t mi) -> const MeshRefs& {
        if (auto it = converted.find(mi); it != converted.end()) return it->second;
        MeshRefs& entry = converted[mi];
        const std::string w = "meshes[" + std::to_string(mi) + "]";
        const Value& m = meshes[mi];
        if (!m.is(Type::Object)) bad(w + " is not an object");
        const Value* nm = m.get("name");
        const std::string name = nm ? nm->string : "mesh";
        const Value& prims = need(m, "primitives", Type::Array, w);
        Batch batch;
        for (std::size_t p = 0; p < prims.array.size(); ++p) {
            const std::string wp = w + ".primitives[" + std::to_string(p) + "]";
            if (!prims.array[p].is(Type::Object)) bad(wp + " is not an object");
            convert_primitive(doc, prims.array[p], name, wp, scene, entry, batch);
        }
        for (const std::uint32_t mi2 : entry.first) compact(scene.meshes[mi2]);
        return entry;
    };

    const auto& nodes = array_of(root, "nodes");
    std::vector<int> parent(nodes.size(), -1);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const std::string w = "nodes[" + std::to_string(i) + "]";
        if (!nodes[i].is(Type::Object)) bad(w + " is not an object");
        if (const Value* ch = opt(nodes[i], "children", Type::Array, w))
            for (const Value& c : ch->array) {
                const std::size_t ci = as_index(c, nodes.size(), w + ".children");
                if (parent[ci] != -1 || ci == i) bad("node " + std::to_string(ci) + " has more than one parent (glTF requires a strict tree)");
                parent[ci] = int(i);
            }
    }
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::size_t steps = 0;
        for (int p = parent[i]; p != -1; p = parent[std::size_t(p)])
            if (++steps > nodes.size()) bad("node hierarchy contains a cycle");
    }

    scene.nodes.emplace_back().name = "root";
    const std::size_t nacc = array_of(root, "accessors").size();
    auto add_node = [&](auto&& self, std::size_t ni, std::uint32_t parent_id) -> void {
        const Value& n = nodes[ni];
        const std::string w = "nodes[" + std::to_string(ni) + "]";
        Node out;
        if (const Value* nm = n.get("name")) out.name = nm->string;
        out.transform = node_transform(n, w);
        const Value* ex = n.get("extensions");
        const Value* inst = ex ? ex->get("EXT_mesh_gpu_instancing") : nullptr;
        const Value* mv = n.get("mesh");
        const std::size_t mi = mv ? as_index(*mv, meshes.size(), w + ".mesh") : 0;
        if (mv && !inst) {
            const MeshRefs& e = mesh_of(mi);
            out.meshes = e.first;
            out.polylines = e.second;
        }
        const auto id = std::uint32_t(scene.nodes.size());
        scene.nodes.push_back(std::move(out));
        scene.nodes[parent_id].children.push_back(id);
        if (mv && inst) {  // EXT_mesh_gpu_instancing: one child node per instance (T * R * S)
            const Value& attrs = need(*inst, "attributes", Type::Object, w + ".EXT_mesh_gpu_instancing");
            std::vector<double> T, R, S;
            int c = 0;
            if (const Value* a = attrs.get("TRANSLATION")) T = doc.values(as_index(*a, nacc, w + ".TRANSLATION"), &c);
            if (const Value* a = attrs.get("ROTATION")) R = doc.values(as_index(*a, nacc, w + ".ROTATION"), &c);
            if (const Value* a = attrs.get("SCALE")) S = doc.values(as_index(*a, nacc, w + ".SCALE"), &c);
            const std::size_t count = !T.empty() ? T.size() / 3 : !R.empty() ? R.size() / 4 : S.size() / 3;
            if ((!T.empty() && T.size() / 3 != count) || (!R.empty() && R.size() / 4 != count) || (!S.empty() && S.size() / 3 != count))
                bad(w + " instancing attributes have different counts");
            const MeshRefs& e = mesh_of(mi);
            for (std::size_t k = 0; k < count; ++k) {
                double t[3] = {0, 0, 0}, q[4] = {0, 0, 0, 1}, s[3] = {1, 1, 1};
                if (!T.empty()) std::copy_n(&T[k * 3], 3, t);
                if (!R.empty()) std::copy_n(&R[k * 4], 4, q);
                if (!S.empty()) std::copy_n(&S[k * 3], 3, s);
                Node in;
                in.name = scene.nodes[id].name + "#" + std::to_string(k);
                in.transform = trs(t, q, s);
                in.meshes = e.first;
                in.polylines = e.second;
                const auto iid = std::uint32_t(scene.nodes.size());
                scene.nodes.push_back(std::move(in));
                scene.nodes[id].children.push_back(iid);
            }
        }
        if (const Value* ch = n.get("children"))
            for (const Value& c : ch->array) self(self, std::size_t(c.integer), id);
    };
    const auto& scenes = array_of(root, "scenes");
    auto add_scene = [&](std::size_t si) {
        const Value& sc = scenes[si];
        if (!sc.is(Type::Object)) bad("scenes[" + std::to_string(si) + "] is not an object");
        if (const Value* roots = opt(sc, "nodes", Type::Array, "scene"))
            for (const Value& r : roots->array) {
                const std::size_t ni = as_index(r, nodes.size(), "scene.nodes");
                if (parent[ni] != -1) bad("scene root node " + std::to_string(ni) + " has a parent");
                add_node(add_node, ni, 0);
            }
    };
    if (const Value* sv = root.get("scene")) add_scene(as_index(*sv, scenes.size(), "root.scene"));
    else if (!scenes.empty()) add_scene(0);
    else
        for (std::size_t i = 0; i < nodes.size(); ++i)
            if (parent[i] == -1) add_node(add_node, i, 0);
    if (scene.nodes.size() == 1)  // no nodes at all: show every mesh untransformed
        for (std::size_t mi = 0; mi < meshes.size(); ++mi) {
            const MeshRefs& e = mesh_of(mi);
            scene.nodes[0].meshes.insert(scene.nodes[0].meshes.end(), e.first.begin(), e.first.end());
            scene.nodes[0].polylines.insert(scene.nodes[0].polylines.end(), e.second.begin(), e.second.end());
        }
    return scene;
}

}  // namespace

const Decoder& gltf() {
    static const SimpleDecoder d({"gltf", "glTF 2.0 / GLB", {"gltf", "glb"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
