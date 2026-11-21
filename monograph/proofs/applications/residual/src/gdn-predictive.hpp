#pragma once
#include "gdn-lossless.hpp"
#include <cfenv>
#include <cmath>
#include <limits>
#if defined(__x86_64__) && defined(__GNUC__)
#include <immintrin.h>
#define GDN_PRED_X86 1
#else
#define GDN_PRED_X86 0
#endif

// A low-rank predictor reduces the exact integer residual, not state precision.
// Compile without fast-math. Encoder and decoder use the same FMA order.
namespace gdn_predictive {
constexpr size_t n = 16384, width = 128, max_rank = 32;
constexpr size_t factor_bytes = 516, header_bytes = 8;
inline uint32_t word(float x) { uint32_t u; std::memcpy(&u, &x, 4); return u; }
inline void put16(uint8_t * p, uint16_t x) { p[0] = uint8_t(x); p[1] = uint8_t(x >> 8); }
inline uint16_t get16(const uint8_t * p) { return uint16_t(p[0]) | (uint16_t(p[1]) << 8); }
inline bool simd_available() {
#if GDN_PRED_X86
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
    return false;
#endif
}
inline void update_scalar(float * p, float * r, const float * u, const float * v) {
    for (size_t i = 0; i < width; ++i) {
        for (size_t j = 0; j < width; ++j) {
            const size_t k = i * width + j;
            p[k] = std::fma(u[i], v[j], p[k]);
            if (r) { r[k] = std::fma(-u[i], v[j], r[k]); }
        }
    }
}
inline size_t pivot_scalar(const float * r) {
    size_t index = 0; float best = std::fabs(r[0]);
    for (size_t i = 1; i < n; ++i) {
        if (std::fabs(r[i]) > best) { best = std::fabs(r[i]); index = i; }
    }
    return index;
}
#if GDN_PRED_X86
__attribute__((target("avx2,fma"))) inline void update_avx(float * p, float * r, const float * u, const float * v) {
    for (size_t i = 0; i < width; ++i) {
        const __m256 a = _mm256_set1_ps(u[i]);
        for (size_t j = 0; j < width; j += 8) {
            const size_t k = i * width + j;
            const __m256 b = _mm256_loadu_ps(v + j);
            _mm256_storeu_ps(p + k, _mm256_fmadd_ps(a, b, _mm256_loadu_ps(p + k)));
            if (r) { _mm256_storeu_ps(r + k, _mm256_fnmadd_ps(a, b, _mm256_loadu_ps(r + k))); }
        }
    }
}
__attribute__((target("avx2,fma"))) inline size_t pivot_avx(const float * r) {
    const __m256 mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    __m256 best = _mm256_setzero_ps();
    for (size_t i = 0; i < n; i += 8) { best = _mm256_max_ps(best, _mm256_and_ps(mask, _mm256_loadu_ps(r + i))); }
    float lanes[8]; _mm256_storeu_ps(lanes, best);
    const float maximum = *std::max_element(lanes, lanes + 8);
    best = _mm256_set1_ps(maximum);
    for (size_t i = 0; i < n; i += 8) {
        const unsigned hits = unsigned(_mm256_movemask_ps(_mm256_cmp_ps(best, _mm256_and_ps(mask, _mm256_loadu_ps(r + i)), _CMP_EQ_OQ)));
        if (hits) { return i + unsigned(__builtin_ctz(hits)); }
    }
    return 0;
}
#endif
inline void update(float * p, float * r, const float * u, const float * v, bool simd) {
#if GDN_PRED_X86
    if (simd) { update_avx(p, r, u, v); return; }
#else
    (void) simd;
#endif
    update_scalar(p, r, u, v);
}
inline size_t pivot(const float * r, bool simd) {
#if GDN_PRED_X86
    if (simd) { return pivot_avx(r); }
#else
    (void) simd;
#endif
    return pivot_scalar(r);
}
inline bool quantize_vector(float * x, uint8_t * out) {
    float maximum = 0;
    for (size_t i = 0; i < width; ++i) { maximum = std::max(maximum, std::fabs(x[i])); }
    int exponent = 0;
    if (maximum) {
        std::frexp(maximum, &exponent); exponent -= 15;
        if (std::ldexp(double(maximum), -exponent) > 32767) { ++exponent; }
    }
    if (exponent < -126 || exponent > 100) { return false; }
    put16(out, uint16_t(int16_t(exponent)));
    for (size_t i = 0; i < width; ++i) {
        const int16_t code = int16_t(std::nearbyint(std::ldexp(double(x[i]), -exponent)));
        put16(out + 2 + 2 * i, uint16_t(code));
        x[i] = std::ldexp(float(code), exponent);
    }
    return true;
}
inline void decode_vector(const uint8_t * in, float * x) {
    const uint16_t e = get16(in);
    const int exponent = e <= 32767 ? int(e) : int(e) - 65536;
    // Valid predictors stay in this range; malformed snapshots must not overflow.
    if (exponent < -126 || exponent > 100) { throw std::runtime_error("invalid predictor exponent"); }
    for (size_t i = 0; i < width; ++i) {
        const uint16_t code = get16(in + 2 + 2 * i);
        x[i] = std::ldexp(float(code <= 32767 ? int(code) : int(code) - 65536), exponent);
    }
}
struct bit_writer {
    uint8_t * p; uint64_t bits = 0; unsigned have = 0;
    void put(uint32_t x, unsigned b) {
        if (!b) { return; }
        bits |= uint64_t(x) << have; have += b;
        while (have >= 8) { *p++ = uint8_t(bits); bits >>= 8; have -= 8; }
    }
    void finish() { if (have) { *p++ = uint8_t(bits); bits = 0; have = 0; } }
};
struct bit_reader {
    const uint8_t * p; uint64_t bits = 0; unsigned have = 0;
    uint32_t get(unsigned b) {
        if (!b) { return 0; }
        while (have < b) { bits |= uint64_t(*p++) << have; have += 8; }
        const uint32_t x = uint32_t(bits & ((uint64_t(1) << b) - 1));
        bits >>= b; have -= b; return x;
    }
};
inline size_t residual_bytes(unsigned b, size_t exceptions) {
    return n / 8 + n * b / 8 + (exceptions ? n / 8 : 0) + (exceptions * (32 - b) + 7) / 8;
}
inline bool valid_environment() {
    if (std::fegetround() != FE_TONEAREST) { return false; }
#if GDN_PRED_X86
    if (_mm_getcsr() & ((1u << 15) | (1u << 6))) { return false; }
#endif
    return true;
}
inline size_t encode(const float * x, uint8_t * out, bool use_simd = true) {
    size_t best = gdn_lossless::encode(x, out);
    if (!best || !valid_environment()) { return best; }
    for (size_t i = 0; i < n; ++i) {
        const uint32_t w = word(x[i]); const unsigned e = (w >> 23) & 255;
        if ((w & 0x7fffffff) && (e < 65 || e > 190)) { return best; }
    }
    const bool simd = use_simd && simd_available();
    float residual[n], prediction[n] = {}, u[width], v[width];
    uint32_t z[n]; uint8_t factors[max_rank * factor_bytes];
    std::memcpy(residual, x, sizeof(residual));
    for (unsigned rank = 1; rank <= max_rank; ++rank) {
        const size_t index = pivot(residual, simd), row = index / width, col = index % width;
        const float value = residual[index];
        if (value == 0) { break; }
        for (size_t i = 0; i < width; ++i) { u[i] = residual[i * width + col] / value; v[i] = residual[row * width + i]; }
        uint8_t * factor = factors + (rank - 1) * factor_bytes;
        if (!quantize_vector(u, factor) || !quantize_vector(v, factor + 258)) { break; }
        update(prediction, residual, u, v, simd);
        if (rank & (rank - 1)) { continue; }
        size_t histogram[33] = {};
        for (size_t i = 0; i < n; ++i) {
            const int64_t d = int64_t(word(x[i]) & 0x7fffffff) - int64_t(word(prediction[i]) & 0x7fffffff);
            z[i] = uint32_t(d >= 0 ? 2 * d : -2 * d - 1);
            const unsigned b = z[i] ? 32 - unsigned(__builtin_clz(z[i])) : 0;
            ++histogram[b];
        }
        size_t exceptions = n - histogram[0], length = n * 8, selected_exceptions = 0; unsigned selected_bits = 32;
        for (unsigned b = 0; b <= 32; ++b) {
            const size_t candidate = header_bytes + factor_bytes * rank + residual_bytes(b, exceptions);
            if (candidate < length) { length = candidate; selected_bits = b; selected_exceptions = exceptions; }
            if (b < 32) { exceptions -= histogram[b + 1]; }
        }
        if (length >= best) { continue; }
        best = length; const unsigned b = selected_bits;
        out[0] = 2; out[1] = uint8_t(rank); out[2] = uint8_t(b); out[3] = selected_exceptions != 0;
        put16(out + 4, uint16_t(selected_exceptions)); put16(out + 6, 0);
        std::memcpy(out + header_bytes, factors, factor_bytes * rank);
        uint8_t * signs = out + header_bytes + factor_bytes * rank;
        uint8_t * low = signs + n / 8; uint8_t * bitmap = low + n * b / 8;
        std::memset(signs, 0, n / 8);
        if (selected_exceptions) { std::memset(bitmap, 0, n / 8); }
        bit_writer lows{low}, highs{bitmap + (selected_exceptions ? n / 8 : 0)};
        const uint32_t mask = b == 32 ? 0xffffffffu : (uint32_t(1) << b) - 1;
        for (size_t i = 0; i < n; ++i) {
            signs[i / 8] |= uint8_t((word(x[i]) >> 31) << (i % 8));
            lows.put(z[i] & mask, b);
            if (b < 32 && (z[i] >> b)) {
                bitmap[i / 8] |= uint8_t(1u << (i % 8)); highs.put(z[i] >> b, 32 - b);
            }
        }
        lows.finish(); highs.finish();
        if (size_t(highs.p - out) != best) { throw std::logic_error("predictor length mismatch"); }
    }
    return best;
}
inline void decode(const uint8_t * in, size_t bytes, float * x, bool use_simd = true) {
    if (!bytes || in[0] != 2) { gdn_lossless::decode(in, bytes, x); return; }
    if (bytes < header_bytes || bytes > gdn_lossless::max_encoded) { throw std::runtime_error("invalid predictive size"); }
    const unsigned rank = in[1], b = in[2]; const size_t exceptions = get16(in + 4);
    if (!rank || rank > max_rank || b > 32 || exceptions > n || in[3] != (exceptions != 0) || get16(in + 6) ||
        (b == 32 && exceptions) || bytes != header_bytes + rank * factor_bytes + residual_bytes(b, exceptions)) {
        throw std::runtime_error("invalid predictive header");
    }
    if (!valid_environment()) { throw std::runtime_error("predictive decoder requires round-to-nearest and gradual underflow"); }
    const bool simd = use_simd && simd_available();
    std::memset(x, 0, n * sizeof(float)); float u[width], v[width];
    for (unsigned k = 0; k < rank; ++k) {
        const uint8_t * factor = in + header_bytes + k * factor_bytes;
        decode_vector(factor, u); decode_vector(factor + 258, v);
        update(x, nullptr, u, v, simd);
    }
    const uint8_t * signs = in + header_bytes + rank * factor_bytes;
    const uint8_t * low = signs + n / 8; const uint8_t * bitmap = low + n * b / 8;
    size_t count = 0;
    if (exceptions) {
        for (size_t i = 0; i < n / 8; ++i) { count += unsigned(__builtin_popcount(unsigned(bitmap[i]))); }
        if (count != exceptions) { throw std::runtime_error("invalid exception count"); }
    }
    bit_reader lows{low}, highs{bitmap + (exceptions ? n / 8 : 0)};
    for (size_t i = 0; i < n; ++i) {
        uint32_t z = lows.get(b);
        if (exceptions && ((bitmap[i / 8] >> (i % 8)) & 1)) { z |= highs.get(32 - b) << b; }
        const int64_t d = (z & 1) ? -int64_t(z / 2) - 1 : int64_t(z / 2);
        const int64_t a = int64_t(word(x[i]) & 0x7fffffff) + d;
        if (a < 0 || a > 0x7fffffff) { throw std::runtime_error("invalid magnitude residual"); }
        const uint32_t w = uint32_t(a) | (uint32_t((signs[i / 8] >> (i % 8)) & 1) << 31);
        std::memcpy(x + i, &w, 4);
    }
}
}
