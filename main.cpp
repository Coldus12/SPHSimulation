#include "StdWindow.hpp"
#include "Managed/Managed.hpp"
#include "CPUSim.h"

int main(int, char**) {
    Vltava::CPUSim cpuSim;
    cpuSim.setSimProps();
    cpuSim.run(100);
    cpuSim.printData();
    //std::cout << "----------------------------------------------------" << std::endl;
    //std::cout << "Real data:" << std::endl;
    //Vltava::StdWindow window(1280, 720);


    //try {
        //Vltava::StdWindow window(1280, 720);
    /*} catch ( vk::SystemError & err )
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
    }*/
    return 0;
}