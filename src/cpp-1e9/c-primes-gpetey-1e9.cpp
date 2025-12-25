#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>

using u32 = uint32_t;
using u64 = uint64_t;

#define TEST_BIT(a, i)  ((a)[(i) >> 6] &  (1ULL << ((i) & 63)))
#define CLEAR_BIT(a, i) ((a)[(i) >> 6] &= ~(1ULL << ((i) & 63)))

static std::vector<u32> sieve_odd_bitset(u32 n) {
    if (n < 2) return {};
    if (n == 2) return {2};

    const u32 size = (n >> 1) + 1;                 // index i => value = 2*i + 1
    std::vector<u64> bits((size + 63) / 64, ~0ULL);
    CLEAR_BIT(bits, 0);                            // 1 is not prime

    const u32 limit = (u32)std::sqrt((long double)n);
    for (u32 i = 1, p = 3; p <= limit; ++i, p += 2) {
        if (TEST_BIT(bits, i)) {
            u32 start = (u32)(((u64)p * p) >> 1);  // index of p^2
            for (u32 j = start; j < size; j += p)  CLEAR_BIT(bits, j);
        }
    }

    std::vector<u32> primes;
    primes.reserve(size / 10);
    primes.push_back(2);
    for (u32 i = 1, p = 3; p <= n; ++i, p += 2)
        if (TEST_BIT(bits, i)) primes.push_back(p);

    return primes;
}

int main() {
    constexpr u32 n = 1000000000ULL; // 1 billion limit
    auto primes = sieve_odd_bitset(n);

    std::cout << "Found " << primes.size() << " primes up to " << n << ".\n";
    if (!primes.empty()) {
        std::cout << "Last few primes: ";
        for (size_t i = primes.size() > 5 ? primes.size() - 5 : 0; i < primes.size(); ++i)
            std::cout << primes[i] << ' ';
        std::cout << "\n";
    }
}
