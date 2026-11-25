// Byte buffer whose resize leaves new bytes uninitialised: for outputs that are written in full
// right after sizing (decompression), so no time goes into zero-filling and the pages are first
// touched by the threads that write them.
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

namespace cadcore {

template <class T> struct NoInitAllocator : std::allocator<T> {
    using std::allocator<T>::allocator;
    template <class U> struct rebind { using other = NoInitAllocator<U>; };
    template <class U> void construct(U* p) noexcept { ::new (static_cast<void*>(p)) U; }  // default-init: no zeroing
    template <class U, class... A> void construct(U* p, A&&... a) { ::new (static_cast<void*>(p)) U(std::forward<A>(a)...); }
};

using Buffer = std::vector<std::uint8_t, NoInitAllocator<std::uint8_t>>;

}  // namespace cadcore
