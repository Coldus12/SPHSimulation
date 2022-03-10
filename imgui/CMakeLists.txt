target_sources(
        ${PROJECT_NAME}
        PRIVATE
            imgui.cpp
            imgui_demo.cpp
            imgui_draw.cpp
            imgui_tables.cpp
            imgui_widgets.cpp
            imstb_rectpack.h
            imstb_textedit.h
            imstb_truetype.h
            imgui_impl_glfw.h
            imgui_impl_glfw.cpp
            imgui_impl_vulkan.h
            imgui_impl_vulkan.cpp
            imconfig.h
        PUBLIC
            imgui.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#add_library(imgui.h imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp imstb_rectpack.h imstb_textedit.h imstb_truetype.h imgui_impl_glfw.h imgui_impl_glfw.cpp imgui_impl_vulkan.h imgui_impl_vulkan.cpp imconfig.h)