/**
 * @file RandomGen.h
 * @brief Thread-safe random number generator interface.
 *
 * Uses thread_local Mersenne Twisters so that parallel RcppParallel workers
 * each have an independent, independently-seeded RNG.
 */
#pragma once
#include <random>

class RandomGen
{
public:
    /**
     * Compatibility stub – with thread_local RNGs each thread seeds itself
     * automatically on first use.
     */
    static void Seed();
    static int GetInt(int max);
    static double GetUniform01();
};
