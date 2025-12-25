#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

std::vector<uint64_t> sieve(uint64_t n) {
    if (n < 2) return {};
    uint64_t half = n / 2;
    std::vector<uint64_t> b((half >> 6) + 1, ~0ULL), r{2};
    b[0] ^= 1;
    for (uint64_t i = 1, lim = std::sqrt(n) / 2; i <= lim; ++i)
        if (b[i >> 6] >> (i & 63) & 1)
            for (uint64_t j = 2*i*(i+1), s = 2*i+1; j <= half; j += s)
                b[j >> 6] &= ~(1ULL << (j & 63));
    for (size_t i = 0; i < b.size(); ++i)
        for (auto w = b[i]; w; w &= w-1)
            if (auto p = ((i<<6)+__builtin_ctzll(w))*2+1; p <= n) r.push_back(p);
    return r;
}

int main() {
    auto p = sieve(1'000'000'000);  // Generate primes up to 1 billion
    std::cout << p.size() << " primes, last 5: ";
    for (size_t i = p.size()-5; i < p.size(); ++i) std::cout << p[i] << ' ';
}