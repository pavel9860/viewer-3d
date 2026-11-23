// Zip reader written from PKWARE APPNOTE 6.3.x.
//
// Layout used: End Of Central Directory record (signature 0x06054b50) found by
// scanning back from the end (it may be followed by a comment of up to 65535
// bytes); ZIP64 EOCD locator (0x07064b50) / record (0x06064b50) when counts or
// offsets are saturated; central directory headers (0x02014b50) give method,
// CRC-32, sizes and local header offset (ZIP64 extra field 0x0001 supplies the
// 64-bit values); local file headers (0x04034b50) give the data start.
// Methods: 0 stored, 8 deflate. Every extracted entry is checked against its CRC-32.
#include "formats/zip.hpp"

#include <algorithm>
#include <cctype>

#include "cadcore/error.hpp"
#include "codec/inflate.hpp"

namespace cadcore::formats {

namespace {

struct Entry {
    std::string name;
    std::uint16_t flags = 0, method = 0;
    std::uint32_t crc = 0;
    std::uint64_t csize = 0, usize = 0, local_offset = 0;
};

std::string lower(std::string_view s) {
    std::string o(s);
    std::transform(o.begin(), o.end(), o.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    return o;
}

}  // namespace

struct ZipArchive::Impl {
    ByteView data;
    std::vector<Entry> entries;
};

ZipArchive::ZipArchive(ByteView data, unsigned threads) : threads_(threads), impl_(new Impl{data, {}}) {
    try {
        if (data.size < 22) fail(ErrorCode::Malformed, "not a zip archive (too short)");
        // End of central directory: last occurrence within the final 22 + 65535 bytes.
        std::size_t eocd = ByteView::npos;
        const std::size_t lo = data.size > 22 + 65535 ? data.size - 22 - 65535 : 0;
        for (std::size_t p = data.size - 22 + 1; p-- > lo;)
            if (Reader(data, p).u32le() == 0x06054b50) { eocd = p; break; }
        if (eocd == ByteView::npos) fail(ErrorCode::Malformed, "not a zip archive (no end-of-central-directory record)");
        Reader r(data, eocd + 4);
        r.skip(4);  // disk numbers
        std::uint64_t count = r.u16le();
        r.skip(2);
        std::uint64_t cd_size = r.u32le();
        std::uint64_t cd_offset = r.u32le();
        if ((count == 0xFFFF || cd_offset == 0xFFFFFFFF || cd_size == 0xFFFFFFFF) && eocd >= 20 && Reader(data, eocd - 20).u32le() == 0x07064b50) {
            Reader loc(data, eocd - 20 + 8);
            const std::uint64_t z64 = loc.u64le();
            Reader z(data, std::size_t(z64));
            if (z.u32le() != 0x06064b50) fail(ErrorCode::Malformed, "bad ZIP64 end-of-central-directory record", std::size_t(z64));
            z.skip(8 + 2 + 2 + 4 + 4 + 8);
            count = z.u64le();
            cd_size = z.u64le();
            cd_offset = z.u64le();
        }
        if (cd_offset > data.size || cd_size > data.size - cd_offset) fail(ErrorCode::Malformed, "central directory outside the archive");
        if (count > cd_size / 46 + 1) fail(ErrorCode::Malformed, "implausible zip entry count");
        Reader c(data, std::size_t(cd_offset));
        impl_->entries.reserve(std::size_t(count));
        for (std::uint64_t i = 0; i < count; ++i) {
            if (c.u32le() != 0x02014b50) fail(ErrorCode::Malformed, "bad central directory header", c.pos() - 4);
            Entry e;
            c.skip(4);  // version made by, version needed
            e.flags = c.u16le();
            e.method = c.u16le();
            c.skip(4);  // time, date
            e.crc = c.u32le();
            e.csize = c.u32le();
            e.usize = c.u32le();
            const std::uint16_t nlen = c.u16le(), xlen = c.u16le(), clen = c.u16le();
            c.skip(8);  // disk start, internal attrs, external attrs
            e.local_offset = c.u32le();
            const std::uint8_t* name = c.bytes(nlen);
            e.name.assign(reinterpret_cast<const char*>(name), nlen);
            const std::size_t xend = c.pos() + xlen;
            while (c.pos() + 4 <= xend) {  // extra fields: look for ZIP64 (0x0001)
                const std::uint16_t id = c.u16le(), sz = c.u16le();
                const std::size_t next = c.pos() + sz;
                if (id == 0x0001) {
                    if (e.usize == 0xFFFFFFFF && c.pos() + 8 <= next) e.usize = c.u64le();
                    if (e.csize == 0xFFFFFFFF && c.pos() + 8 <= next) e.csize = c.u64le();
                    if (e.local_offset == 0xFFFFFFFF && c.pos() + 8 <= next) e.local_offset = c.u64le();
                }
                c.seek(next);
            }
            c.seek(xend);
            c.skip(clen);
            impl_->entries.push_back(std::move(e));
        }
    } catch (...) {
        delete impl_;
        impl_ = nullptr;
        throw;
    }
}

ZipArchive::~ZipArchive() { delete impl_; }

std::vector<std::string> ZipArchive::names() const {
    std::vector<std::string> out;
    for (const auto& e : impl_->entries) out.push_back(e.name);
    return out;
}

std::optional<Buffer> ZipArchive::read(std::string_view wanted) const {
    while (!wanted.empty() && wanted.front() == '/') wanted.remove_prefix(1);
    const std::string key = lower(wanted);
    for (const Entry& e : impl_->entries) {
        if (lower(e.name) != key) continue;
        if (e.flags & 1) fail(ErrorCode::Unsupported, "encrypted zip entry '" + e.name + "'");
        const ByteView& d = impl_->data;
        Reader l(d, std::size_t(e.local_offset));
        if (l.u32le() != 0x04034b50) fail(ErrorCode::Malformed, "bad local file header for '" + e.name + "'", std::size_t(e.local_offset));
        l.skip(22);
        const std::uint16_t nlen = l.u16le(), xlen = l.u16le();
        l.skip(std::size_t(nlen) + xlen);
        if (e.csize > l.remaining()) fail(ErrorCode::Truncated, "zip entry '" + e.name + "' runs past the archive", l.pos());
        if (e.usize > (std::uint64_t(1) << 32)) fail(ErrorCode::Unsupported, "zip entry larger than 4 GiB");
        const ByteView packed = d.sub(l.pos(), std::size_t(e.csize));
        Buffer out;
        if (e.method == 0) {
            if (e.csize != e.usize) fail(ErrorCode::Malformed, "stored zip entry with different sizes");
            out.assign(packed.data, packed.data + packed.size);
        } else if (e.method == 8) {
            out = codec::inflate_parallel(packed, std::size_t(e.usize), std::size_t(e.usize), threads_);
            if (out.size() != e.usize) fail(ErrorCode::Malformed, "zip entry '" + e.name + "' inflates to the wrong size");
        } else {
            fail(ErrorCode::Unsupported, "zip compression method " + std::to_string(e.method) + " for '" + e.name + "'");
        }
        if (codec::crc32_parallel(out.data(), out.size(), threads_) != e.crc) fail(ErrorCode::Malformed, "CRC-32 mismatch in zip entry '" + e.name + "'");
        return out;
    }
    return std::nullopt;
}

}  // namespace cadcore::formats
