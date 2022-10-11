#include "Test.h"
#include "CPUSim.hpp"

namespace Vltava {

    void Test::runTests() {
        basicPlaceTest();
        tupleBoundsCheck();
        basicNeighbourhoodTest();
    }

    // Grid tests
    //------------------------------------------------------------------------------------------------------------------
    // Testing placement
    void Test::basicPlaceTest() {
        std::cout << "doing basicPlaceTest" << std::endl;

        CPUSim sim1;

        Particle p1 = {
                glm::vec3(0.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
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

        std::vector<Particle> particles;
        particles.push_back(p1);

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // number of cells on a given axis = 10; list_size = 12
        /* --> cellwidth = 0.2
         * --> cellTuple = abs(floor((a.x - x)/0.2), floor((a.y - y)/0.2), floor((a.z - z)/0.2))
         * --> cellTuple = (-1 - 0.45)/0.2 = 7 etc.
         * --> cellTuple = (7, 8, 9)
         */
        glm::vec3 cellTuple = sim1.determineGridTuple(0);
        assert((cellTuple.x == 7) && "CellTuple x isn't correct");
        assert((cellTuple.y == 8) && "CellTuple y isn't correct");
        assert((cellTuple.z == 5) && "CellTuple z isn't correct");

        // Index should be: x * cellz * celly + y * cellz + z
        // Therefore idx = (7 * 10 * 10 + 8 * 10 + 5) * 12 = 785 * 20 = 15 700
        int idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == 15700) && "getStartIdxOfCell returned a wrong number.");

        // Since we only put in 1 particle, size = 1
        int size = sim1.grid_data[idx];
        assert((size == 1) && "size isn't correct");
    }

    void Test::tupleBoundsCheck() {
        std::cout << "doing tupleBoundsCheck" << std::endl;
        CPUSim sim1;

        Particle p1 = {
                glm::vec3(-2.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
        };

        Particle p2 = {
                glm::vec3(2.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
        };

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                2.0f,
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };

        std::vector<Particle> particles;
        particles.push_back(p1);
        particles.push_back(p2);

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // number of cells on a given axis = 10; list_size = 12
        /* --> cellwidth = 0.2
         * --> cellTuple = abs(floor((a.x - x)/0.2), floor((a.y - y)/0.2), floor((a.z - z)/0.2))
         * --> cellTuple = abs( (-1 - (-2.45))/0.2 ) = 7 etc. --> problem
         * --> cellTuple = should be (7, 8, 9), but -2.45 is out of the bounds so
         * --> cellTuple = (-1, -1, -1)
         */
        glm::vec3 cellTuple = sim1.determineGridTuple(0);
        assert((cellTuple.x == -1) && "CellTuple x isn't correct");
        assert((cellTuple.y == -1) && "CellTuple y isn't correct");
        assert((cellTuple.z == -1) && "CellTuple z isn't correct");

        // There shouldn't be a calculated index.
        int idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == -1) && "getStartIdxOfCell returned a wrong number.");

        cellTuple = sim1.determineGridTuple(1);
        assert((cellTuple.x == -1) && "2.CellTuple x isn't correct");
        assert((cellTuple.y == -1) && "2.CellTuple y isn't correct");
        assert((cellTuple.z == -1) && "2.CellTuple z isn't correct");

        idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == -1) && "2.getStartIdxOfCell returned a wrong number.");
    }

    void Test::basicNeighbourhoodTest() {
        std::cout << "doing basicNeighbourhoodTest" << std::endl;

        CPUSim sim1;
        std::vector<Particle> particles;

        float s = 0.2;
        int r = -1;
        int z = -1;
        for (int i = 0; i < 1000; i++) {
            Particle data;

            if (i % 10 == 0)
                r++;

            if (i % (100) == 0)
                z++;

            //data.x = glm::vec3((i%10) * s,(r%10) * s,z * s);
            data.x = glm::vec3((i%10) * s,(r%10) * s,z * s);
            data.x -= glm::vec3(1);
            data.h = 0.1; // Only sets size for visualization atm
            data.v = glm::vec3(0,0,0);
            data.m = 0.005236f;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 0;
            data.padding = 0;

            particles.push_back(data);
        }

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                (float) particles.size(),
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // s = 0.2; gridA = (-1, -1, -1)
        // x = (555%10 = 5) * 0.2 - 1, y = (555%100 = 5) * 0.2 - 1, z = (floor(555/100)) * 0.2 - 1
        // x = y = z = 0 -> middle
        glm::vec3 tuple = sim1.determineGridTuple(555); // middle -> should be 5,5,5
        assert((tuple.x == 5) && "[Neighbourhood] tuple x isnt correct");
        assert((tuple.y == 5) && "[Neighbourhood] tuple y isnt correct");
        assert((tuple.z == 5) && "[Neighbourhood] tuple z isnt correct");

        int neighbournr = 0;
        Neighbourhood nh = sim1.getNeighbouringCells(tuple);
        for (auto& current: nh.neighbour) {
            int idx = sim1.getStartIdxOfCell(current);
            int size = sim1.grid_data[idx];
            neighbournr += size;
        }
        assert((neighbournr == 27) && "[Neighbourhood] the number of neighbours isn't correct");
    }
}