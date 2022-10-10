#include "Test.h"
#include "CPUSim.hpp"

namespace Vltava {

    void Test::runTests() {
        testPlacefunc();
    }

    void Test::testPlacefunc() {
        CPUSim sim1;

        int cellSize = 10;
        int list_size = 12;

        Particle p1 = {
                glm::vec3(0,0,0),
                0.1f,
                glm::vec3(0,0,0),
                0.536f,
                0.0f,
                0.0f,
                0.005236f
        };

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                1.0f,
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };
    }
}