#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <utility>

class Timer
{
public:
    using clock = std::chrono::steady_clock;
    using time_point = clock::time_point;

    explicit Timer(std::string name = "")
        : name_(std::move(name)), start_(clock::now()) {
    }

    void reset()
    {
        start_ = clock::now();
    }

    [[nodiscard]] double elapsed_ns() const
    {
        return std::chrono::duration<double, std::nano>(
            clock::now() - start_).count();
    }

    [[nodiscard]] double elapsed_us() const
    {
        return std::chrono::duration<double, std::micro>(
            clock::now() - start_).count();
    }

    [[nodiscard]] double elapsed_ms() const
    {
        return std::chrono::duration<double, std::milli>(
            clock::now() - start_).count();
    }

    [[nodiscard]] double elapsed_s() const
    {
        return std::chrono::duration<double>(
            clock::now() - start_).count();
    }

    void print_ms(const std::string& prefix = "") const
    {
        if (!prefix.empty())
            std::cout << prefix;
        else if (!name_.empty())
            std::cout << "[" << name_ << "] ";

        std::cout << elapsed_ms() << " ms" << std::endl;
    }

    void print_s(const std::string& prefix = "") const
    {
        if (!prefix.empty())
            std::cout << prefix;
        else if (!name_.empty())
            std::cout << "[" << name_ << "] ";

        std::cout << elapsed_s() << " s" << std::endl;
    }

private:
    std::string name_;
    time_point start_;
};

class ScopedTimer
{
public:
    using clock = std::chrono::steady_clock;

    explicit ScopedTimer(std::string name)
        : name_(std::move(name)), start_(clock::now()) {
    }

    ~ScopedTimer()
    {
        const double ms = std::chrono::duration<double, std::milli>(
            clock::now() - start_).count();
        std::cout << "[" << name_ << "] " << ms << " ms" << std::endl;
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
    std::string name_;
    clock::time_point start_;
};