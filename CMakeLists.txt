cmake_minimum_required(VERSION 3.19)
project(SPHSimulation)

set(CMAKE_CXX_STANDARD 17)

add_executable(SPHSimulation main.cpp StdWindow.cpp StdWindow.hpp Vertex.cpp Vertex.hpp Model.cpp Model.hpp VltavaFunctions.hpp VltavaFunctions.cpp Buffer.cpp Buffer.hpp ComputeShader.cpp ComputeShader.h)
#add_subdirectory( imgui )
#add_subdirectory( vma )

if (WIN32)
    add_subdirectory( glfw )
    add_subdirectory( glm )
endif(WIN32)

find_package(Vulkan REQUIRED)
target_include_directories(
        ${PROJECT_NAME}
        PUBLIC ${Vulkan_INCLUDE_DIRS}
)

if (WIN32)
    target_link_libraries(${PROJECT_NAME} Vulkan::Vulkan glfw glm)
endif(WIN32)

if (UNIX)
    target_link_libraries(${PROJECT_NAME} Vulkan::Vulkan glfw3 glm)

    find_package(glm CONFIG REQUIRED)
    target_link_libraries(${PROJECT_NAME} glm)

    find_package(glfw3 REQUIRED)
    target_link_libraries(${PROJECT_NAME} glfw3)
endif(UNIX)