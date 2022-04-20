cmake_minimum_required(VERSION 3.19)
project(SPHSimulation)

set(CMAKE_CXX_STANDARD 17)

add_executable(SPHSimulation main.cpp StdWindow.cpp StdWindow.hpp Vertex.cpp Vertex.hpp Model.cpp Model.hpp VltavaFunctions.hpp VltavaFunctions.cpp Buffer.cpp Buffer.hpp ComputeShader.cpp ComputeShader.hpp Material.cpp Material.hpp ParticleModel.cpp ParticleModel.hpp Managed/MInstance.cpp Managed/MInstance.hpp Managed/MSurface.cpp Managed/MSurface.hpp Managed/MPhysDev.cpp Managed/MPhysDev.hpp Managed/MLogDev.cpp Managed/MLogDev.hpp Managed/MSwapChain.cpp Managed/MSwapChain.hpp Managed/MImageView.hpp Managed/MRenderPass.hpp Managed/MFramebuffer.hpp Managed/MCommandPool.hpp Managed/MCommandBuffers.hpp Managed/MSemaphore.hpp Managed/MFence.hpp Managed/Managed.hpp Managed/MBuffer.hpp Managed/MDeviceMemory.hpp Managed/MPipelineLayout.hpp Managed/MPipeline.hpp Managed/MDescriptorSetLayout.hpp Managed/MDescriptorPool.hpp Managed/MDescriptorSets.hpp Managed/MShaderModule.hpp Managed/MDescriptorSet.hpp)
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
    target_link_libraries(${PROJECT_NAME} Vulkan::Vulkan glfw glm -static-libgcc -static-libstdc++ -static -lpthread)
endif(WIN32)

if (UNIX)
    target_link_libraries(${PROJECT_NAME} Vulkan::Vulkan glfw3 glm)

    find_package(glm CONFIG REQUIRED)
    target_link_libraries(${PROJECT_NAME} glm)

    find_package(glfw3 REQUIRED)
    target_link_libraries(${PROJECT_NAME} glfw3)
endif(UNIX)