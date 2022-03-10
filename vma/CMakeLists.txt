target_sources(
        ${PROJECT_NAME}
        PUBLIC
            vk_mem_alloc.hpp
            vk_mem_alloc.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})