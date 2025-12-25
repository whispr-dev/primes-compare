#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = uint64_t;
using u32 = uint32_t;

std::vector<u32> base_sieve(u32 n) {
    u32 h = n / 2 + 1;
    std::vector<u64> b((h + 63) >> 6, ~0ULL);
    b[0] ^= 1;
    for (u32 i = 1, L = (u32)std::sqrt(n) / 2; i <= L; ++i)
        if (b[i >> 6] >> (i & 63) & 1)
            for (u32 j = 2*i*(i+1), s = 2*i+1; j < h; j += s)
                b[j >> 6] &= ~(1ULL << (j & 63));
    std::vector<u32> P{2};
    for (u32 i = 0; i < b.size(); ++i)
        for (auto w = b[i]; w; w &= w - 1) {
            u32 v = ((i << 6) + __builtin_ctzll(w)) * 2 + 1;
            if (v > 1 && v <= n) P.push_back(v);
        }
    return P;
}

int main() {
    u64 n = 1'000'000'000ULL;
    constexpr u32 S = 1 << 17;
    auto B = base_sieve((u32)std::sqrt((double)n) + 1);
    std::vector<u64> seg((S + 63) >> 6);
    std::array<u64, 5> ring{};
    u64 cnt = 0, pos = 0;

    auto add = [&](u64 p) { ++cnt; ring[pos++ % 5] = p; };
    if (n >= 2) add(2);

    for (u64 lo = 3; lo <= n; lo += S << 1) {
        u64 hi = std::min(lo + (S << 1) - 2, n);
        std::fill(seg.begin(), seg.end(), ~0ULL);
        for (size_t i = 1; i < B.size(); ++i) {
            u64 p = B[i], s = std::max(p * p, ((lo + p - 1) / p) * p);
            if (!(s & 1)) s += p;
            for (u64 j = s; j <= hi; j += p << 1)
                seg[(j - lo) >> 7] &= ~(1ULL << (((j - lo) >> 1) & 63));
        }
        for (size_t i = 0; i < seg.size(); ++i)
            for (auto w = seg[i]; w; w &= w - 1) {
                u64 v = lo + (((i << 6) + __builtin_ctzll(w)) << 1);
                if (v <= n) add(v);
            }
    }

    std::cout << "Found " << cnt << " primes up to " << n << ".\nLast 5: ";
    auto k = std::min(cnt, 5ULL);
    for (size_t i = 0; i < k; ++i) std::cout << ring[(pos + 5 - k + i) % 5] << ' ';
    std::cout << '\n';
}