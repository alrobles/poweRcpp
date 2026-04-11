/**
 * @file RandomGen.cpp
 * @brief Thread-safe random number generator implementation.
 *
 * Uses thread_local storage so that each thread in an RcppParallel
 * context gets its own independently-seeded Mersenne Twister, avoiding
 * data races.  A single shared thread_local generator is used by both
 * GetInt() and GetUniform01() so that all variates on a given thread
 * come from the same RNG stream.
 */
#include "RandomGen.h"
#include <thread>
#include <functional>

/* One Mersenne Twister per thread, seeded from random_device XOR'd with the
 * thread ID hash to reduce the chance of identical seeds across threads. */
static std::mt19937& thread_rng()
{
    thread_local static std::mt19937 gen(
        static_cast<unsigned>(
            std::random_device{}() ^
            (std::hash<std::thread::id>{}(std::this_thread::get_id()) << 1u)));
    return gen;
}

/* Seeding only affects subsequent calls on the calling thread. */
void RandomGen::Seed()
{
    /* No-op for thread-local variant; each thread auto-seeds on first use. */
}

int RandomGen::GetInt(int max)
{
    std::uniform_int_distribution<> dist(0, max);
    return dist(thread_rng());
}

double RandomGen::GetUniform01()
{
    std::uniform_real_distribution<> dist(0.0, 1.0);
    return dist(thread_rng());
}
