#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace shipsim
{
    inline void setupGlobalThreads(int reservedCores = 2)
    {
        const int hw = static_cast<int>(std::thread::hardware_concurrency());
        const int nThread = std::max(1, hw - std::max(0, reservedCores));

#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(nThread);
        omp_set_nested(0);
#endif
        Eigen::setNbThreads(nThread);

        std::cout << "[Thread] hardware=" << hw
                  << ", OpenMP=" << nThread
                  << ", Eigen=" << nThread << "\n";
    }

    class EigenSingleThreadGuard
    {
    public:
        EigenSingleThreadGuard() : old_(Eigen::nbThreads())
        {
            Eigen::setNbThreads(1);
        }
        ~EigenSingleThreadGuard()
        {
            Eigen::setNbThreads(old_);
        }
        EigenSingleThreadGuard(const EigenSingleThreadGuard&) = delete;
        EigenSingleThreadGuard& operator=(const EigenSingleThreadGuard&) = delete;
    private:
        int old_;
    };
}
