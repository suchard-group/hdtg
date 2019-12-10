//
// Created by Marc Suchard on 2019-12-04.
//

#ifndef ZIG_ZAG_TIMING_H
#define ZIG_ZAG_TIMING_H

// Set-up C++11 clock profiling support since Travis-CI does not yet have std::chrono::steady_clock

#include <ctime>
#include <chrono>

namespace zz {
    namespace chrono {
        typedef std::chrono::high_resolution_clock steady_clock; // Travis does not support steady_clock
//        typedef std::chrono::nanoseconds TimingUnits; // Travis-CI does not support using aliases
        typedef std::chrono::microseconds TimingUnits;

        using std::chrono::system_clock;
        using std::chrono::duration_cast;
        using std::chrono::duration;

    } // namespace chrono
} /* namespace bsccs */

#endif //ZIG_ZAG_TIMING_H
