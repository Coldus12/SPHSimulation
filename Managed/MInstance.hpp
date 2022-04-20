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
        std::unique_ptr<vk::Instance> instance;
        std::unique_ptr<vk::DebugUtilsMessengerEXT> debugUtilsMessenger;
        bool enableValidationLayers;

        const std::vector<const char *> validationLayers = {
                "VK_LAYER_KHRONOS_validation"
        };

        PFN_vkCreateDebugUtilsMessengerEXT pfnVkCreateDebugUtilsMessengerEXT = VK_NULL_HANDLE;
        PFN_vkDestroyDebugUtilsMessengerEXT pfnVkDestroyDebugUtilsMessengerEXT = VK_NULL_HANDLE;

        static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(
                VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                VkDebugUtilsMessageTypeFlagsEXT type,
                const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData,
                void *pUserData);

        VKAPI_ATTR VkResult VKAPI_CALL vkCreateDebugUtilsMessengerEXT(
                VkInstance instance,
                const VkDebugUtilsMessengerCreateInfoEXT *pCreateInfo,
                const VkAllocationCallbacks *pAllocator,
                VkDebugUtilsMessengerEXT *pDebugMessenger);

        VKAPI_ATTR void VKAPI_CALL vkDestroyDebugUtilsMessengerEXT(
                VkInstance instance,
                VkDebugUtilsMessengerEXT debugMessenger,
                const VkAllocationCallbacks *pAllocator);

        bool checkValidationLayerSupport();
        void setupDebugMessenger();
        std::vector<const char *> getRequiredExtensions();
    };
}

#endif //SPHSIMULATION_MINSTANCE_HPP