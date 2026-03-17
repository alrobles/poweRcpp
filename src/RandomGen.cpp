/**
 * @file RandomGen.cpp
 * @brief Thread-safe random number generator implementation.
 *
 * Uses thread_local storage so that each thread in an RcppParallel
 * context gets its own independently-seeded Mersenne Twister, avoiding
 * data races.
 */
#include "RandomGen.h"
#include <thread>
#include <functional>
using namespace std;

/* Seeding only affects subsequent calls on the calling thread. */
void RandomGen::Seed()
{
    /* No-op for thread-local variant; each thread auto-seeds on first use. */
}

int RandomGen::GetInt(int max)
{
    thread_local static mt19937 gen(
        static_cast<unsigned>(
            random_device{}() ^
            (hash<thread::id>{}(this_thread::get_id()) << 1u)));
    uniform_int_distribution<> dist(0, max);
    return dist(gen);
}

double RandomGen::GetUniform01()
{
    thread_local static mt19937 gen(
        static_cast<unsigned>(
            random_device{}() ^
            (hash<thread::id>{}(this_thread::get_id()) << 1u)));
    uniform_real_distribution<> dist(0.0, 1.0);
    return dist(gen);
}
