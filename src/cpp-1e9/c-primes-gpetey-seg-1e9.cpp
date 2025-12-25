#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

using u32 = uint32_t;
using u64 = uint64_t;

static inline bool bit_test(const std::vector<u64>& b, u64 i) {
    return (b[i >> 6] >> (i & 63)) & 1ULL;
}
static inline void bit_clear(std::vector<u64>& b, u64 i) {
    b[i >> 6] &= ~(1ULL << (i & 63));
}

// Odd-only simple sieve up to limit (returns primes including 2).
static std::vector<u32> base_sieve(u32 limit) {
    if (limit < 2) return {};
    if (limit == 2) return {2};

    const u32 size = (limit >> 1) + 1;                 // odds index space
    std::vector<u64> bits((size + 63) / 64, ~0ULL);
    bit_clear(bits, 0);                                // 1 is not prime

    const u32 root = (u32)std::sqrt((long double)limit);
    for (u32 i = 1, p = 3; p <= root; ++i, p += 2) {
        if (bit_test(bits, i)) {
            u32 start = (u32)(((u64)p * p) >> 1);
            for (u32 j = start; j < size; j += p) bit_clear(bits, j);
        }
    }

    std::vector<u32> primes;
    primes.reserve(size / 10);
    primes.push_back(2);
    for (u32 i = 1, p = 3; p <= limit; ++i, p += 2)
        if (bit_test(bits, i)) primes.push_back(p);
    return primes;
}

// Count primes up to n using odd-only segmented sieve.
// seg_odds = number of odd candidates per segment.
static u64 segmented_count(u64 n, u32 seg_odds = 1u << 20) {
    if (n < 2) return 0;
    if (n == 2) return 1;

    const u64 root = (u64)std::sqrt((long double)n);
    auto primes = base_sieve((u32)root);

    // Odd numbers from 3..n inclusive:
    const u64 total_odds = (n >= 3) ? (((n - 3) >> 1) + 1) : 0;
    const u64 num_segs  = (total_odds + seg_odds - 1) / seg_odds;

    u64 total = 1; // prime 2

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:total)
    #endif
    for (u64 s = 0; s < num_segs; ++s) {
        const u64 low  = 3 + 2 * s * (u64)seg_odds;     // odd
        u64 high = low + 2 * (u64)seg_odds;             // exclusive-ish bound
        const u64 n1 = n + 1;
        if (high > n1) high = n1;
        if (high <= low) continue;

        const u64 len_odds = (high - low) >> 1;         // count of odds in [low, high)
        if (!len_odds) continue;

        std::vector<u64> bits((len_odds + 63) / 64, ~0ULL);

        // Mark composites using odd base primes (skip 2).
        for (size_t k = 1; k < primes.size(); ++k) {
            const u64 p = primes[k];
            const u64 step = p << 1;

            // start = max(p*p, ceil(low/p)*p), then make it odd
            u64 start = (u64)((__int128)p * p);
            if (start < low) start = ((low + p - 1) / p) * p;
            if ((start & 1) == 0) start += p;

            if (start >= high) continue;

            for (u64 m = start; m < high; m += step) {
                const u64 idx = (m - low) >> 1;
                bit_clear(bits, idx);
            }
        }

        // Popcount all set bits (mask tail).
        u64 c = 0;
        const u64 tail = len_odds & 63;
        for (size_t i = 0; i < bits.size(); ++i) {
            u64 w = bits[i];
            if (tail && i + 1 == bits.size())
                w &= (tail == 64 ? ~0ULL : ((1ULL << tail) - 1));
            c += (u64)__builtin_popcountll(w);
        }

        total += c;
    }

    return total;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    u64 n = 1000000000ULL;               // default: 1e9
    u32 seg_odds = 1u << 20;             // ~1,048,576 odds => ~128 KiB bitset

    if (argc > 1) n = std::stoull(argv[1]);
    if (argc > 2) seg_odds = (u32)std::stoul(argv[2]);

    const u64 count = segmented_count(n, seg_odds);

    std::cout << "π(" << n << ") = " << count << "\n";
    #ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << "\n";
    #endif
}
