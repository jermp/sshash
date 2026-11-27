#pragma once

#include <cstring>
#if defined(__x86_64__)
#include <x86intrin.h>
#endif

#include "include/util.hpp"

namespace sshash {
namespace lookup_kernels {

/*
    The shared front end of the batched lookup engines: packing a query
    string into a kmer, and computing the minimizer of one isolated kmer.
    Both are compile-time gated on the instruction set and on the kmer/hash
    types, and fall back to the scalar code of util.hpp elsewhere; the
    vector minimizer kernels fall back to the scalar `compute_minimizer`
    also on any hash tie, so its tie-breaking stays the single source of
    truth.
*/

/* Pack the k valid characters of str into a kmer.

   With BMI2, eight characters per step: for both nucleotide encodings the
   code of a character is a function of bits 1-3 of that same byte, so the
   whole word can be shifted and the two code bits of every byte extracted
   in one pext. The tail of k % 8 characters stays scalar, so no byte
   beyond str[k-1] is ever read. */
template <typename Kmer>
[[maybe_unused]] static Kmer pack_string(char const* str, const uint64_t k) {
#if defined(__BMI2__)
    if constexpr (Kmer::bits_per_char == 2 and Kmer::has_reverse_complement) {
        static_assert(Kmer::uint_kmer_bits % 64 == 0);
        constexpr uint64_t num_words = Kmer::uint_kmer_bits / 64;
        uint64_t words[num_words] = {};
        const uint64_t num_chunks = k / 8;
        for (uint64_t j = 0; j != num_chunks; ++j) {
            uint64_t v;
            std::memcpy(&v, str + 8 * j, 8);
#ifdef SSHASH_USE_TRADITIONAL_NUCLEOTIDE_ENCODING
            const uint64_t y = (v >> 1) ^ (v >> 2);
#else
            const uint64_t y = v >> 1;
#endif
            const uint64_t codes = _pext_u64(y, 0x0303030303030303ULL);
            words[(16 * j) >> 6] |= codes << ((16 * j) & 63);
        }
        for (uint64_t i = 8 * num_chunks; i != k; ++i) {
            words[(2 * i) >> 6] |= Kmer::char_to_uint(str[i]) << ((2 * i) & 63);
        }
        Kmer x = 0;
        std::memcpy(&x.bits, words, sizeof(x.bits));
        assert(x == util::string_to_uint_kmer<Kmer>(str, k));
        return x;
    }
#endif
    return util::string_to_uint_kmer<Kmer>(str, k);
}

#if defined(__AVX512F__) && defined(__AVX512DQ__)
/*
    As minimizer_avx2 below, but eight loci per iteration, with the native
    64-bit multiply and unsigned compares into mask registers; loci beyond
    the window in the last iteration are excluded by the lane mask instead
    of a scalar tail.
*/
template <typename Kmer>
[[maybe_unused]] static minimizer_info minimizer_avx512(Kmer const& kmer, Kmer const& kmer_rc,
                                                        const uint64_t k, const uint64_t m,
                                                        hasher_type const& hasher)  //
{
    static_assert(Kmer::uint_kmer_bits % 64 == 0);
    constexpr uint64_t bits_per_char = Kmer::bits_per_char;
    const uint64_t w = k - m + 1;
    constexpr uint64_t num_words = Kmer::uint_kmer_bits / 64;
    alignas(8) uint64_t words[num_words + 2];
    alignas(8) uint64_t words_rc[num_words + 2];
    std::memcpy(words, &kmer.bits, sizeof(kmer.bits));
    std::memcpy(words_rc, &kmer_rc.bits, sizeof(kmer_rc.bits));
    words[num_words] = words[num_words + 1] = 0;
    words_rc[num_words] = words_rc[num_words + 1] = 0;

    const uint64_t mb = m * bits_per_char;
    const uint64_t mask64 = (uint64_t(1) << mb) - 1;
    const __m512i mask = _mm512_set1_epi64(mask64);
    const __m512i mulc = _mm512_set1_epi64(0x517cc1b727220a95ULL);
    const uint64_t magic64 = hasher.hash(0);
    const __m512i magic = _mm512_set1_epi64(magic64);
    const __m512i c63 = _mm512_set1_epi64(63);
    const __m512i c64 = _mm512_set1_epi64(64);
    const __m512i bitr_base = _mm512_set1_epi64(int64_t(2 * (k - m)));

    __m512i bit = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);
    __m512i posv = _mm512_setr_epi64(0, 1, 2, 3, 4, 5, 6, 7);
    __m512i min_hash = _mm512_set1_epi64(-1);  // unsigned max
    __m512i min_value = _mm512_setzero_si512();
    __m512i min_pos = _mm512_setzero_si512();
    __mmask8 tie = 0;

    /* extract eight m-mers at bit positions `at` from `src`; lanes whose
       position falls outside the kmer yield garbage and are masked out of
       the minimum below */
    auto extract = [&](uint64_t const* src, const __m512i at) {
        if constexpr (num_words == 1) {
            return _mm512_and_si512(_mm512_srlv_epi64(_mm512_set1_epi64(int64_t(src[0])), at),
                                    mask);
        } else {
            const __mmask8 in_hi = _mm512_cmpgt_epi64_mask(at, c63);
            const __m512i w0 = _mm512_mask_blend_epi64(in_hi, _mm512_set1_epi64(int64_t(src[0])),
                                                       _mm512_set1_epi64(int64_t(src[1])));
            const __m512i w1 = _mm512_mask_blend_epi64(in_hi, _mm512_set1_epi64(int64_t(src[1])),
                                                       _mm512_set1_epi64(int64_t(src[2])));
            const __m512i sh = _mm512_and_si512(at, c63);
            return _mm512_and_si512(
                _mm512_or_si512(_mm512_srlv_epi64(w0, sh),
                                _mm512_sllv_epi64(w1, _mm512_sub_epi64(c64, sh))),
                mask);
        }
    };

    for (uint64_t i = 0; i < w; i += 8) {
        const __mmask8 valid =
            w - i >= 8 ? __mmask8(0xff) : __mmask8((uint32_t(1) << (w - i)) - 1);
        const __m512i mmer = extract(words, bit);
        const __m512i r = extract(words_rc, _mm512_sub_epi64(bitr_base, bit));
        const __mmask8 rc_lt = _mm512_cmplt_epu64_mask(r, mmer);
        const __m512i value = _mm512_mask_blend_epi64(rc_lt, mmer, r);
        const __m512i hash = _mm512_xor_si512(_mm512_mullo_epi64(value, mulc), magic);
        tie |= _mm512_mask_cmpeq_epi64_mask(valid, hash, min_hash);
        const __mmask8 lt = _mm512_mask_cmplt_epu64_mask(valid, hash, min_hash);
        min_hash = _mm512_mask_blend_epi64(lt, min_hash, hash);
        min_value = _mm512_mask_blend_epi64(lt, min_value, value);
        min_pos = _mm512_mask_blend_epi64(lt, min_pos, posv);
        bit = _mm512_add_epi64(bit, _mm512_set1_epi64(16));
        posv = _mm512_add_epi64(posv, _mm512_set1_epi64(8));
    }

    if (tie != 0) return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);

    alignas(64) uint64_t hs[8], vs[8], ps[8];
    _mm512_store_si512(reinterpret_cast<__m512i*>(hs), min_hash);
    _mm512_store_si512(reinterpret_cast<__m512i*>(vs), min_value);
    _mm512_store_si512(reinterpret_cast<__m512i*>(ps), min_pos);

    uint64_t best_hash = hs[0], best_value = vs[0], best_pos = ps[0];
    uint64_t num_mins = 1;
    for (int l = 1; l != 8; ++l) {
        if (hs[l] < best_hash) {
            best_hash = hs[l];
            best_value = vs[l];
            best_pos = ps[l];
            num_mins = 1;
        } else if (hs[l] == best_hash) {
            num_mins += 1;
        }
    }
    if (num_mins > 1) return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
    minimizer_info info(best_value, best_pos);
    assert([&] {
        auto expected = util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
        return minimizer_info(info.minimizer, info.pos_in_kmer) ==
               minimizer_info(expected.minimizer, expected.pos_in_kmer);
    }());
    return info;
}
#endif

#if defined(__AVX2__)
/*
    The minimizer of one isolated kmer, four loci per iteration: the m-mers
    of loci i..i+3 are extracted with variable shifts from the kmer's packed
    words, their reverse complements from the precomputed kmer_rc at the
    mirror loci, canonicalized by unsigned minimum, hashed with the mixer,
    and folded into a per-lane running minimum; a scalar pass reduces the
    four lanes and the tail loci. Any hash equality along the way raises a
    (conservative) tie flag, on which the computation falls back to the
    generic `compute_minimizer`, whose tie-breaking is the single source of
    truth; true ties are vanishingly rare.
*/
template <typename Kmer>
[[maybe_unused]] static minimizer_info minimizer_avx2(Kmer const& kmer, Kmer const& kmer_rc,
                                                      const uint64_t k, const uint64_t m,
                                                      hasher_type const& hasher)  //
{
    static_assert(Kmer::uint_kmer_bits % 64 == 0);
    constexpr uint64_t bits_per_char = Kmer::bits_per_char;
    const uint64_t w = k - m + 1;
    if (w < 8) return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
    constexpr uint64_t num_words = Kmer::uint_kmer_bits / 64;
    alignas(32) uint64_t words[num_words + 2];
    std::memcpy(words, &kmer.bits, sizeof(kmer.bits));
    words[num_words] = 0;
    words[num_words + 1] = 0;

    const uint64_t mb = m * bits_per_char;
    const uint64_t mask64 = (uint64_t(1) << mb) - 1;
    const __m256i mask = _mm256_set1_epi64x(mask64);
    const __m256i mulc = _mm256_set1_epi64x(0x517cc1b727220a95ULL);
    const uint64_t magic64 = hasher.hash(0);
    const __m256i magic = _mm256_set1_epi64x(magic64);
    const __m256i sign = _mm256_set1_epi64x(int64_t(uint64_t(1) << 63));
    const __m256i c64 = _mm256_set1_epi64x(64);

    /* the m-mer of locus i read on the reverse strand is the m-mer of
       kmer_rc at the mirror locus k - m - i, i.e. at bit 2*(k - m) - 2*i */
    alignas(32) uint64_t words_rc[num_words + 2];
    std::memcpy(words_rc, &kmer_rc.bits, sizeof(kmer_rc.bits));
    words_rc[num_words] = 0;
    words_rc[num_words + 1] = 0;
    const __m256i bitr_base = _mm256_set1_epi64x(int64_t(2 * (k - m)));

    __m256i min_hash = _mm256_set1_epi64x(-1);  // unsigned max
    __m256i min_value = _mm256_setzero_si256();
    __m256i min_pos = _mm256_setzero_si256();
    __m256i tie = _mm256_setzero_si256();

    auto mullo64 = [](__m256i a, __m256i b) {
        const __m256i lo = _mm256_mul_epu32(a, b);
        const __m256i ah = _mm256_srli_epi64(a, 32);
        const __m256i bh = _mm256_srli_epi64(b, 32);
        const __m256i cross = _mm256_add_epi64(_mm256_mul_epu32(ah, b), _mm256_mul_epu32(a, bh));
        return _mm256_add_epi64(lo, _mm256_slli_epi64(cross, 32));
    };

    uint64_t i = 0;
    for (; i + 4 <= w; i += 4) {
        const __m256i bit = _mm256_setr_epi64x(int64_t(2 * i), int64_t(2 * (i + 1)),
                                               int64_t(2 * (i + 2)), int64_t(2 * (i + 3)));

        /* extract four m-mers at bit positions `at` from `src` */
        auto extract = [&](uint64_t const* src, const __m256i at) {
            if constexpr (num_words == 1) {
                return _mm256_and_si256(
                    _mm256_srlv_epi64(_mm256_set1_epi64x(int64_t(src[0])), at), mask);
            } else {
                /* per-lane word pair, stitched with variable shifts; shift
                   counts >= 64 yield zero, which handles the boundaries */
                const __m256i in_hi = _mm256_cmpgt_epi64(at, _mm256_set1_epi64x(63));
                const __m256i w0 = _mm256_blendv_epi8(_mm256_set1_epi64x(int64_t(src[0])),
                                                      _mm256_set1_epi64x(int64_t(src[1])), in_hi);
                const __m256i w1 = _mm256_blendv_epi8(_mm256_set1_epi64x(int64_t(src[1])),
                                                      _mm256_set1_epi64x(int64_t(src[2])), in_hi);
                const __m256i sh = _mm256_and_si256(at, _mm256_set1_epi64x(63));
                return _mm256_and_si256(
                    _mm256_or_si256(_mm256_srlv_epi64(w0, sh),
                                    _mm256_sllv_epi64(w1, _mm256_sub_epi64(c64, sh))),
                    mask);
            }
        };

        const __m256i mmer = extract(words, bit);
        const __m256i r = extract(words_rc, _mm256_sub_epi64(bitr_base, bit));

        /* canonical = unsigned min(mmer, rc) */
        const __m256i rc_lt =
            _mm256_cmpgt_epi64(_mm256_xor_si256(mmer, sign), _mm256_xor_si256(r, sign));
        const __m256i value = _mm256_blendv_epi8(mmer, r, rc_lt);
        const __m256i hash = _mm256_xor_si256(mullo64(value, mulc), magic);

        /* per-lane running minimum, with equality detection */
        tie = _mm256_or_si256(tie, _mm256_cmpeq_epi64(hash, min_hash));
        const __m256i lt = _mm256_cmpgt_epi64(_mm256_xor_si256(min_hash, sign),
                                              _mm256_xor_si256(hash, sign));
        min_hash = _mm256_blendv_epi8(min_hash, hash, lt);
        min_value = _mm256_blendv_epi8(min_value, value, lt);
        const __m256i pos =
            _mm256_setr_epi64x(int64_t(i), int64_t(i + 1), int64_t(i + 2), int64_t(i + 3));
        min_pos = _mm256_blendv_epi8(min_pos, pos, lt);
    }

    if (_mm256_movemask_epi8(tie) != 0) {
        return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
    }

    alignas(32) uint64_t hs[4], vs[4], ps[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(hs), min_hash);
    _mm256_store_si256(reinterpret_cast<__m256i*>(vs), min_value);
    _mm256_store_si256(reinterpret_cast<__m256i*>(ps), min_pos);

    /* the tail loci, scalar, and the final reduction */
    uint64_t best_hash = hs[0], best_value = vs[0], best_pos = ps[0];
    uint64_t num_mins = 1;
    for (int l = 1; l != 4; ++l) {
        if (hs[l] < best_hash) {
            best_hash = hs[l];
            best_value = vs[l];
            best_pos = ps[l];
            num_mins = 1;
        } else if (hs[l] == best_hash) {
            num_mins += 1;
        }
    }
    for (; i != w; ++i) {
        const uint64_t bit = 2 * i;
        uint64_t mmer = words[bit >> 6] >> (bit & 63);
        if ((bit & 63) != 0) mmer |= words[(bit >> 6) + 1] << (64 - (bit & 63));
        mmer &= mask64;
        const uint64_t value = std::min(mmer, Kmer::reverse_complement_mmer(mmer, m));
        const uint64_t hash = (value * 0x517cc1b727220a95ULL) ^ magic64;
        if (hash < best_hash) {
            best_hash = hash;
            best_value = value;
            best_pos = i;
            num_mins = 1;
        } else if (hash == best_hash) {
            num_mins += 1;
        }
    }
    if (num_mins > 1) return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
    minimizer_info info(best_value, best_pos);
    assert([&] {
        auto expected = util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
        return minimizer_info(info.minimizer, info.pos_in_kmer) ==
               minimizer_info(expected.minimizer, expected.pos_in_kmer);
    }());
    return info;
}
#endif

/* The dispatcher: the fastest compiled-in kernel for DNA kmers under the
   mixer hash, the scalar code otherwise. `no_avx512` is a runtime escape
   hatch for machines where 512-bit multiplies downclock the core. */
template <typename Kmer>
[[maybe_unused]] static minimizer_info compute_minimizer(Kmer const& kmer, Kmer const& kmer_rc,
                                                         const uint64_t k, const uint64_t m,
                                                         hasher_type const& hasher,
                                                         [[maybe_unused]] const bool no_avx512)  //
{
    if constexpr (Kmer::has_reverse_complement and Kmer::bits_per_char == 2 and
                  std::is_same_v<hasher_type, mixer_64>) {
#if defined(__AVX512F__) && defined(__AVX512DQ__)
        if (!no_avx512) return minimizer_avx512(kmer, kmer_rc, k, m, hasher);
#endif
#if defined(__AVX2__)
        return minimizer_avx2(kmer, kmer_rc, k, m, hasher);
#endif
    }
    return util::compute_minimizer<Kmer>(kmer, kmer_rc, k, m, hasher);
}

}  // namespace lookup_kernels
}  // namespace sshash
