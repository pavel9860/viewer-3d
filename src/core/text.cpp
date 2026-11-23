#include "core/text.hpp"

#include <algorithm>
#include <cstring>

#include "core/parallel.hpp"

namespace cadcore::text {

std::string_view as_utf8(std::string_view raw, std::string& storage) {
    auto b = [&](std::size_t i) { return static_cast<unsigned char>(raw[i]); };
    if (raw.size() >= 3 && b(0) == 0xEF && b(1) == 0xBB && b(2) == 0xBF) return raw.substr(3);
    if (raw.size() < 2) return raw;
    const bool le = b(0) == 0xFF && b(1) == 0xFE, be = b(0) == 0xFE && b(1) == 0xFF;
    if (!le && !be) return raw;
    storage.clear();
    storage.reserve(raw.size() / 2);
    auto unit = [&](std::size_t i) { return le ? char16_t(b(i) | b(i + 1) << 8) : char16_t(b(i) << 8 | b(i + 1)); };
    for (std::size_t i = 2; i + 1 < raw.size(); i += 2) {
        std::uint32_t cp = unit(i);
        if (cp >= 0xD800 && cp < 0xDC00 && i + 3 < raw.size()) {
            const std::uint32_t lo = unit(i + 2);
            if (lo >= 0xDC00 && lo < 0xE000) { cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00); i += 2; }
        }
        if (cp < 0x80) storage += char(cp);
        else if (cp < 0x800) { storage += char(0xC0 | cp >> 6); storage += char(0x80 | (cp & 63)); }
        else if (cp < 0x10000) { storage += char(0xE0 | cp >> 12); storage += char(0x80 | (cp >> 6 & 63)); storage += char(0x80 | (cp & 63)); }
        else { storage += char(0xF0 | cp >> 18); storage += char(0x80 | (cp >> 12 & 63)); storage += char(0x80 | (cp >> 6 & 63)); storage += char(0x80 | (cp & 63)); }
    }
    return storage;
}

std::vector<std::size_t> line_starts(std::string_view t, std::size_t from, unsigned threads) {
    const std::size_t n = t.size() - from;
    const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>(thread_count(threads), n >> 20));
    std::vector<Slot<std::vector<std::size_t>>> part(parts);
    parallel_for(parts, 1, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t k = b0; k < e0; ++k) {
            const char* lo = t.data() + from + n * k / parts;
            const char* hi = t.data() + from + n * (k + 1) / parts;
            auto& v = part[k].value;
            if (k == 0) v.push_back(from);
            for (const char* p = lo; (p = static_cast<const char*>(std::memchr(p, '\n', std::size_t(hi - p)))) != nullptr; ++p)
                if (std::size_t(p - t.data()) + 1 < t.size()) v.push_back(std::size_t(p - t.data()) + 1);
        }
    });
    std::vector<std::size_t> all;
    append_all(all, parts, [&](std::size_t k) -> auto& { return part[k].value; }, threads);
    return all;
}

std::size_t parallel_line_records(const std::vector<std::size_t>& lines, std::size_t pos, std::size_t count, std::size_t text_size, unsigned threads,
                                  std::size_t& parts, const std::function<std::size_t(std::size_t, std::size_t, std::size_t, std::size_t)>& parse) {
    parts = count < 8192 ? 1 : thread_count(threads);
    if (parts < 2) return std::string_view::npos;
    const auto first = std::size_t(std::lower_bound(lines.begin(), lines.end(), pos) - lines.begin());
    if (first >= lines.size() || lines[first] != pos || first + count > lines.size()) return std::string_view::npos;
    auto line_at = [&](std::size_t r) { return first + r < lines.size() ? lines[first + r] : text_size; };
    std::vector<Slot<bool>> ok(parts);
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < parts; ++k)
        tasks.push_back([&, k] {
            const std::size_t a = count * k / parts, b = count * (k + 1) / parts;
            try {
                ok[k].value = parse(k, line_at(a), a, b - a) == line_at(b);
            } catch (const DecodeError&) {
                ok[k].value = false;  // the sequential parse decides
            }
        });
    run_all(tasks, threads);
    for (auto& o : ok)
        if (!o.value) return std::string_view::npos;
    return line_at(count);
}

}  // namespace cadcore::text
