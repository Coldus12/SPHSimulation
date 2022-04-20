#include "VltavaFunctions.hpp"

namespace Vltava {
    uint32_t findMemoryType(uint32_t typeFilter, const vk::PhysicalDeviceMemoryProperties &properties, vk::MemoryPropertyFlags flags) {
        for (uint32_t i = 0; i < properties.memoryTypeCount; i++) {
            if ((typeFilter & (1 << i)) && (properties.memoryTypes[i].propertyFlags & flags) == flags) {
                return i;
            }
        }

        throw std::runtime_error("failed to find suitable memory type!");
    }

    std::vector<char> readFile(const std::string &filename) {
        std::ifstream file(filename, std::ios::ate | std::ios::binary);
        if (!file.is_open())
            throw std::runtime_error("Failed to open file!");

        size_t fileSize = (size_t) file.tellg();
        std::vector<char> buffer(fileSize);

        file.seekg(0);
        file.read(buffer.data(), fileSize);

        file.close();

        return buffer;
    }

    // Function for easier debugging
    void sanityCheck(const std::string &str) {
        std::cout << "[Sanity check] " << str << std::endl;
    }

    /*PFN_vkCreateDebugUtilsMessengerEXT pfnVkCreateDebugUtilsMessengerEXT = VK_NULL_HANDLE;
    PFN_vkDestroyDebugUtilsMessengerEXT pfnVkDestroyDebugUtilsMessengerEXT = VK_NULL_HANDLE;

    VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                                            VkDebugUtilsMessageTypeFlagsEXT type,
                                                            const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData,
                                                            void *pUserData) {

        std::cerr << "[Validation layer] " << pCallbackData->pMessage << std::endl;
        return VK_FALSE;
    }

    VKAPI_ATTR VkResult VKAPI_CALL vkCreateDebugUtilsMessengerEXT(
            VkInstance instance,
            const VkDebugUtilsMessengerCreateInfoEXT *pCreateInfo,
            const VkAllocationCallbacks *pAllocator,
            VkDebugUtilsMessengerEXT *pDebugMessenger) {

        return pfnVkCreateDebugUtilsMessengerEXT(instance, pCreateInfo, pAllocator, pDebugMessenger);
    }

    VKAPI_ATTR void VKAPI_CALL vkDestroyDebugUtilsMessengerEXT(
            VkInstance instance,
            VkDebugUtilsMessengerEXT debugMessenger,
            const VkAllocationCallbacks *pAllocator) {

        return pfnVkDestroyDebugUtilsMessengerEXT(instance, debugMessenger, pAllocator);
    }*/
}
