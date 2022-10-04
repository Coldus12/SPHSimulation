#include "StdWindow.hpp"
#include "Managed/Managed.hpp"
#include "CPUSim.hpp"

#define cpuTest true

int main(int, char**) {
#if cpuTest
    Vltava::CPUSim cpuSim;
    cpuSim.setSimProps();
    /*for (int i = 0; i < 2; i++) {
        cpuSim.run(50);
        //cpuSim.printData();

        std::cout << i <<std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
    }*/
    cpuSim.run(3);
    cpuSim.run(1);
    cpuSim.printData();
    /*std::cout << "--------------------------------------" << std::endl;
    cpuSim.run(1);
    cpuSim.printData();*/
    /*cpuSim.run(86);
    std::cout << "----------------------------------------------------" << std::endl;
    cpuSim.printData();
    std::cout << "\n\nUtolso" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    cpuSim.run(1);
    cpuSim.printData();*/
#else
    try {
        Vltava::StdWindow window(1280, 720);
    } catch ( vk::SystemError & err )
    {
        std::cout << "vk::SystemError: " << err.what() << std::endl;
        exit( -1 );
    }
    catch ( std::exception & err )
    {
        std::cout << "std::exception: " << err.what() << std::endl;
        exit( -1 );
    }
    catch ( ... )
    {
        std::cout << "unknown error\n";
        exit( -1 );
    }
#endif
    return 0;
}