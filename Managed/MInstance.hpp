#ifndef SPHSIMULATION_MINSTANCE_HPP
#define SPHSIMULATION_MINSTANCE_HPP


#include <memory>
#include <GLFW/glfw3.h>
#include <vulkan/vulkan.hpp>
//#include "../VltavaFunctions.hpp"

namespace Vltava {
    class MInstance {
    public:
        MInstance(std::string name = "test", bool enableValidationLayers = true);
        ~MInstance();

        vk::Instance getHandle();
        std::vector<const char *> getValidationLayers();
    private:
        vk::Instance instance;
        vk::DebugUtilsMessengerEXT debugUtilsMessenger;
        bool enableValidationLayers;

        const std::vector<const char *> validationLayers = {
                "VK_LAYER_KHRONOS_validation"
        };

        bool checkValidationLayerSupport();
        void setupDebugMessenger();
        std::vector<const char *> getRequiredExtensions();
    };
}

#endif //SPHSIMULATION_MINSTANCE_HPP