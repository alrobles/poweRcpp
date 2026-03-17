/**
 * @file VectorUtilities.h
 * @brief Template utility functions for std::vector operations.
 *
 * Adapted from the alrobles/powerlaw repository.
 */
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>

namespace VectorUtilities
{
    template<typename T>
    void RemoveLower(std::vector<T>& v, T n)
    {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&](const T& val){ return val < n; }), v.end());
    }

    template<typename T>
    void RemoveLowerOrEqual(std::vector<T>& v, T n)
    {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&](const T& val){ return val <= n; }), v.end());
    }

    template<typename T>
    void RemoveGreater(std::vector<T>& v, T n)
    {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&](const T& val){ return val > n; }), v.end());
    }

    template<typename T>
    void RemoveGreaterOrEqual(std::vector<T>& v, T n)
    {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&](const T& val){ return val >= n; }), v.end());
    }

    template<typename T>
    void Sort(std::vector<T>& v)
    {
        std::sort(v.begin(), v.end());
    }

    template<typename T>
    void Insert(std::vector<T>& v1, const std::vector<T>& v2)
    {
        v1.insert(v1.end(), v2.begin(), v2.end());
    }

    template<typename T>
    int IndexOf(const std::vector<T>& v, T n)
    {
        return static_cast<int>(
            std::upper_bound(v.begin(), v.end(), n) - v.begin());
    }

    template<typename T>
    int NumberOfGreater(const std::vector<T>& v, T n)
    {
        return static_cast<int>(
            std::count_if(v.begin(), v.end(),
                          [&](const T& val){ return val > n; }));
    }

    template<typename T>
    int NumberOfLower(const std::vector<T>& v, T n)
    {
        return static_cast<int>(
            std::count_if(v.begin(), v.end(),
                          [&](const T& val){ return val < n; }));
    }

    template<typename T>
    int NumberOfGreaterOrEqual(const std::vector<T>& v, T n)
    {
        return static_cast<int>(
            std::count_if(v.begin(), v.end(),
                          [&](const T& val){ return val >= n; }));
    }

    template<typename T>
    int NumberOfLowerOrEqual(const std::vector<T>& v, T n)
    {
        return static_cast<int>(
            std::count_if(v.begin(), v.end(),
                          [&](const T& val){ return val <= n; }));
    }

    template<typename T>
    T Max(const std::vector<T>& v)
    {
        return *std::max_element(v.begin(), v.end());
    }

    template<typename T>
    int IndexOfMax(const std::vector<T>& v)
    {
        return static_cast<int>(
            std::max_element(v.begin(), v.end()) - v.begin());
    }

    template<typename T>
    T Min(const std::vector<T>& v)
    {
        return *std::min_element(v.begin(), v.end());
    }

    template<typename T>
    int IndexOfMin(const std::vector<T>& v)
    {
        return static_cast<int>(
            std::min_element(v.begin(), v.end()) - v.begin());
    }
}
