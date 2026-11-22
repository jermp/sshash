#pragma once

#include <cstdint>
#include <type_traits>

#ifdef __linux__
#include <sys/mman.h>
#endif

namespace sshash {

/*
    Advise the kernel to back every large buffer of a structure with huge
    pages, by walking the structure's own `visit` traversal (the one its
    serialization uses, so nothing is missed: the packed text, the control
    codewords, the locate sets, and the internals of the minimizers and
    kmers MPHFs alike). On an index of gigabytes accessed at random, 4KB
    pages make every cache miss drag a TLB miss and a page-table walk with
    it; 2MB pages mostly remove that.

    Linux only (a no-op elsewhere), and advisory: MADV_COLLAPSE (when the
    kernel has it) collapses the ranges synchronously, and MADV_HUGEPAGE
    lets khugepaged collapse them in the background otherwise. Failures are
    ignored.
*/
struct hugepage_advisor  //
{
    template <typename T, typename V, typename = void>
    struct has_visit : std::false_type {};
    template <typename T, typename V>
    struct has_visit<T, V,
                     std::void_t<decltype(std::declval<T const&>().visit(std::declval<V&>()))>>
        : std::true_type {};

    template <typename T, typename = void>
    struct is_contiguous : std::false_type {};
    template <typename T>
    struct is_contiguous<
        T, std::enable_if_t<std::is_pointer_v<decltype(std::declval<T const&>().data())> and
                            std::is_integral_v<decltype(std::declval<T const&>().size())>>>
        : std::true_type {};

    template <typename T>
    void visit(T const& t) {
        if constexpr (std::is_arithmetic_v<T> or std::is_enum_v<T>) {
            /* nothing to advise */
        } else if constexpr (is_contiguous<T>::value) {
            using E = std::remove_const_t<std::remove_pointer_t<decltype(t.data())>>;
            advise(t.data(), t.size() * sizeof(E));
            if constexpr (has_visit<E, hugepage_advisor>::value) {
                for (uint64_t i = 0; i != t.size(); ++i) t.data()[i].visit(*this);
            }
        } else if constexpr (has_visit<T, hugepage_advisor>::value) {
            t.visit(*this);
        } else {
            /* a plain aggregate: nothing to advise */
        }
    }

    uint64_t advised_bytes() const { return m_collapsed_bytes + m_lazy_bytes; }
    uint64_t collapsed_bytes() const { return m_collapsed_bytes; }
    uint64_t lazy_bytes() const { return m_lazy_bytes; }

private:
    uint64_t m_collapsed_bytes = 0;  // synchronously collapsed (MADV_COLLAPSE)
    uint64_t m_lazy_bytes = 0;       // advised only (khugepaged, eventually)

    void advise(void const* ptr, const uint64_t bytes) {
#ifdef __linux__
        constexpr uint64_t huge = uint64_t(1) << 21;
        const uint64_t addr = reinterpret_cast<uint64_t>(ptr);
        const uint64_t begin = (addr + huge - 1) & ~(huge - 1);
        const uint64_t end = (addr + bytes) & ~(huge - 1);
        if (end <= begin) return;
        void* p = reinterpret_cast<void*>(begin);
#ifdef MADV_COLLAPSE
        if (madvise(p, end - begin, MADV_COLLAPSE) == 0) {
            m_collapsed_bytes += end - begin;
            return;
        }
#endif
        if (madvise(p, end - begin, MADV_HUGEPAGE) == 0) m_lazy_bytes += end - begin;
#else
        (void)ptr;
        (void)bytes;
#endif
    }
};

/* Advise huge pages for every large buffer reachable from `t`. Returns the
   number of bytes advised. */
template <typename T>
hugepage_advisor advise_hugepages(T const& t) {
    hugepage_advisor advisor;
    t.visit(advisor);
    return advisor;
}

}  // namespace sshash
