#include "codec/draco/kd_tree.hpp"

#include <limits>
#include <memory>

namespace cadcore::codec::draco {
namespace {

// Bit sources used by the compression levels: raw 32-bit words read MSB first ("direct"),
// rABS ("RAnsBitDecoder") and 32 rABS coders, one per bit position ("folded").
class BitSource {
public:
    virtual ~BitSource() = default;
    virtual void start(Stream& s) = 0;
    virtual bool bit() = 0;
    virtual std::uint32_t bits(int n) = 0;
};

class DirectBits final : public BitSource {
public:
    void start(Stream& s) override {
        const std::uint32_t size = s.u32();
        if (size == 0 || (size & 3) || size > s.remaining()) bad("invalid direct bit buffer size");
        words_.resize(size / 4);
        for (auto& w : words_) w = s.u32();
        pos_ = 0;
        used_ = 0;
    }
    bool bit() override {  // past the end: 0, as the encoder's reader does
        if (pos_ >= words_.size()) return false;
        const bool b = (words_[pos_] >> (31 - used_)) & 1u;
        if (++used_ == 32) { ++pos_; used_ = 0; }
        return b;
    }
    std::uint32_t bits(int n) override {
        if (n <= 0) return 0;
        const int remaining = 32 - int(used_);
        if (n <= remaining ? pos_ >= words_.size() : pos_ + 1 >= words_.size()) bad("direct bit buffer exhausted");
        std::uint32_t v = 0;
        for (int i = 0; i < n; ++i) v = (v << 1) | std::uint32_t(bit());
        return v;
    }

private:
    std::vector<std::uint32_t> words_;
    std::size_t pos_ = 0;
    std::uint32_t used_ = 0;
};

class RansBits final : public BitSource {
public:
    void start(Stream& s) override { d_.start(s); }
    bool bit() override { return d_.next(); }
    std::uint32_t bits(int n) override {
        std::uint32_t v = 0;
        for (int i = 0; i < n; ++i) v = (v << 1) + std::uint32_t(d_.next());
        return v;
    }

private:
    BitDecoder d_;
};

class FoldedBits final : public BitSource {
public:
    void start(Stream& s) override {
        for (auto& d : per_bit_) d.start(s);
        bit_.start(s);
    }
    bool bit() override { return bit_.next(); }
    std::uint32_t bits(int n) override {
        if (n > 32) bad("folded bit read longer than 32 bits");
        std::uint32_t v = 0;
        for (int i = 0; i < n; ++i) v = (v << 1) + std::uint32_t(per_bit_[std::size_t(i)].next());
        return v;
    }

private:
    BitDecoder per_bit_[32];
    BitDecoder bit_;
};

int msb(std::uint32_t n) { return 31 - __builtin_clz(n); }

}  // namespace

void decode_kd_tree_points(Stream& s, int level, std::uint32_t dimension, std::uint32_t expected_points, std::vector<std::uint32_t>& out) {
    if (level < 0 || level > 6) bad("kd-tree compression level out of range");
    if (dimension == 0) bad("kd-tree without dimensions");
    out.assign(std::size_t(expected_points) * dimension, 0);
    const std::uint32_t bit_length = s.u32();
    if (bit_length > 32) bad("kd-tree bit length over 32");
    const std::uint32_t num_points = s.u32();
    if (num_points != expected_points) bad("kd-tree point count differs from the header");
    if (num_points == 0) return;

    std::unique_ptr<BitSource> numbers;
    if (level >= 4) numbers = std::make_unique<FoldedBits>();
    else if (level >= 2) numbers = std::make_unique<RansBits>();
    else numbers = std::make_unique<DirectBits>();
    DirectBits remaining_bits, axis_bits, half_bits;
    const bool select_axis = level == 6;
    numbers->start(s);
    remaining_bits.start(s);
    axis_bits.start(s);
    half_bits.start(s);

    struct Status { std::uint32_t remaining, last_axis, depth; };
    struct Undo { std::uint32_t depth, axis, old_level, old_base; };
    std::vector<std::uint32_t> base(dimension, 0), levels(dimension, 0), p(dimension, 0), axes(dimension, 0);
    std::vector<Undo> undo;
    std::vector<Status> stack{{num_points, 0, 0}};
    std::uint32_t depth = 0, decoded = 0;
    auto emit = [&](const std::vector<std::uint32_t>& v) {
        if (decoded >= num_points) bad("kd-tree decodes more points than declared");
        std::copy(v.begin(), v.end(), out.begin() + std::ptrdiff_t(std::size_t(decoded) * dimension));
        ++decoded;
    };
    while (!stack.empty()) {
        const Status st = stack.back();
        stack.pop_back();
        if (st.depth != depth) {
            if (st.depth > depth) bad("kd-tree depth jumps forward");
            while (!undo.empty()) {
                const Undo& u = undo.back();
                if (u.depth < st.depth) break;
                base[u.axis] = u.old_base;
                if (u.depth == st.depth) break;
                levels[u.axis] = u.old_level;
                undo.pop_back();
            }
            depth = st.depth;
        }
        if (st.remaining > num_points) bad("kd-tree node larger than the point count");
        std::uint32_t axis;
        if (!select_axis) {
            axis = st.last_axis + 1 == dimension ? 0 : st.last_axis + 1;
        } else if (st.remaining < 64) {
            axis = 0;
            for (std::uint32_t a = 1; a < dimension; ++a)
                if (levels[axis] > levels[a]) axis = a;
        } else {
            axis = axis_bits.bits(4);
        }
        if (axis >= dimension) bad("kd-tree axis out of range");
        const std::uint32_t lvl = levels[axis];
        if (bit_length - lvl == 0) {
            for (std::uint32_t i = 0; i < st.remaining; ++i) emit(base);
            continue;
        }
        if (st.remaining <= 2) {
            axes[0] = axis;
            for (std::uint32_t i = 1; i < dimension; ++i) axes[i] = axes[i - 1] + 1 == dimension ? 0 : axes[i - 1] + 1;
            for (std::uint32_t i = 0; i < st.remaining; ++i) {
                for (std::uint32_t j = 0; j < dimension; ++j) {
                    const std::uint32_t a = axes[j];
                    const std::uint32_t nbits = bit_length - levels[a];
                    p[a] = base[a] | (nbits ? remaining_bits.bits(int(nbits)) : 0u);
                }
                emit(p);
            }
            continue;
        }
        if (decoded > num_points) bad("kd-tree decodes more points than declared");
        const std::uint32_t nrem_bits = bit_length - lvl;
        const std::uint32_t modifier = 1u << (nrem_bits - 1);
        const std::uint32_t number = numbers->bits(msb(st.remaining));
        std::uint32_t first_half = st.remaining / 2;
        if (first_half < number) bad("kd-tree split count out of range");
        first_half -= number;
        std::uint32_t second_half = st.remaining - first_half;
        if (first_half != second_half && !half_bits.bit()) std::swap(first_half, second_half);
        if (depth == std::numeric_limits<std::uint32_t>::max()) bad("kd-tree too deep");
        undo.push_back({depth, axis, levels[axis], base[axis]});
        levels[axis] += 1;
        base[axis] += modifier;
        if (first_half) stack.push_back({first_half, axis, depth});
        if (second_half) stack.push_back({second_half, axis, depth + 1});
        ++depth;
    }
    if (decoded != num_points) bad("kd-tree decoded fewer points than declared");
}

}  // namespace cadcore::codec::draco
