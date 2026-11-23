// Internal: one accessor per format module; all_decoders() lists them in registration order.
#pragma once

#include <vector>

#include "cadcore/decode.hpp"

namespace cadcore::formats {

const Decoder& stl();
const Decoder& obj();
const Decoder& ply();
const Decoder& off();
const Decoder& gltf();
const Decoder& threemf();
const Decoder& catpart();
const Decoder& catproduct();
const Decoder& threedxml();
const Decoder& cgr();
const Decoder& drc();
const Decoder& dxf();

std::vector<const Decoder*> all_decoders();

// Shared helper: a decoder declared as a table of functions.
class SimpleDecoder final : public Decoder {
public:
    using SniffFn = int (*)(ByteView, std::string_view);
    using DecodeFn = Scene (*)(ByteView, const DecodeOptions&);
    SimpleDecoder(FormatInfo format_info, SniffFn sniff_fn, DecodeFn decode_fn)
        : info_(std::move(format_info)), sniff_(sniff_fn), decode_(decode_fn) {}
    const FormatInfo& info() const override { return info_; }
    int sniff(ByteView d, std::string_view ext) const override { return sniff_(d, ext); }
    Scene decode(ByteView d, const DecodeOptions& o) const override { return decode_(d, o); }

private:
    FormatInfo info_;
    SniffFn sniff_;
    DecodeFn decode_;
};


}  // namespace cadcore::formats
