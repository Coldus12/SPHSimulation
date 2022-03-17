target_sources(
        ${PROJECT_NAME}
        PUBLIC
            VmaUsage.hpp
            VmaUsage.cpp
            vk_mem_alloc.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})