#ifndef SPHSIMULATION_TEST_H
#define SPHSIMULATION_TEST_H

#include "CPUSim.hpp"

namespace Vltava {
    class Test {
    public:
        static void runTests();
        static void basicPlaceTest();
        static void tupleBoundsCheck();
        static void basicNeighbourhoodTest();
        static void cpuGpuPlaceCompare();
    };
}


#endif //SPHSIMULATION_TEST_H
