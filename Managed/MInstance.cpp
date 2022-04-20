#include "MInstance.hpp"
#include "../VltavaFunctions.hpp"

namespace Vltava {
    VKAPI_ATTR VkBool32 VKAPI_CALL MInstance::debugCallback(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                                            VkDebugUtilsMessageTypeFlagsEXT type,
                                                            const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData,
                                                            void *pUserData) {

        std::cerr << "[Validation layer] " << pCallbackData->pMessage << std::endl;
        return VK_FALSE;
    }

    VKAPI_ATTR VkResult VKAPI_CALL MInstance::vkCreateDebugUtilsMessengerEXT(VkInstance instance,
                                                                             const VkDebugUtilsMessengerCreateInfoEXT *pCreateInfo,
                                                                             const VkAllocationCallbacks *pAllocator,
                                                                             VkDebugUtilsMessengerEXT *pDebugMessenger) {

        return pfnVkCreateDebugUtilsMessengerEXT(instance, pCreateInfo, pAllocator, pDebugMessenger);
    }

    VKAPI_ATTR void VKAPI_CALL MInstance::vkDestroyDebugUtilsMessengerEXT(VkInstance instance,
                                                                          VkDebugUtilsMessengerEXT debugMessenger,
                                                                          const VkAllocationCallbacks *pAllocator) {

        return pfnVkDestroyDebugUtilsMessengerEXT(instance, debugMessenger, pAllocator);
    }


    MInstance::MInstance(std::string name, bool enableValidationLayers): enableValidationLayers(enableValidationLayers) {
        vk::ApplicationInfo appInfo(
                name.c_str(),
                VK_MAKE_VERSION(1, 0, 0),
                "No Engine",
                VK_MAKE_VERSION(1, 0, 0),
                VK_API_VERSION_1_0
        );

        vk::InstanceCreateInfo createInfo({}, &appInfo);

        auto glfwExtensions = getRequiredExtensions();

        if (enableValidationLayers) {
            if (!checkValidationLayerSupport())
                throw std::runtime_error("Validation layers requested, but not available!");

            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
            createInfo.ppEnabledLayerNames = validationLayers.data();

            vk::DebugUtilsMessengerCreateInfoEXT debugCreateInfo(
                    {},
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eError,
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance,
                    debugCallback
            );

            createInfo.pNext = (vk::DebugUtilsMessengerCreateInfoEXT *) &debugCreateInfo;
        } else {
            createInfo.enabledLayerCount = 0;
            //createInfo.pNext = nullptr;
        }

        createInfo.enabledExtensionCount = static_cast<uint32_t>(glfwExtensions.size());
        createInfo.ppEnabledExtensionNames = glfwExtensions.data();

        createInfo.pNext = nullptr;

        instance = std::make_unique<vk::Instance>(vk::createInstance(createInfo));
        setupDebugMessenger();
    }

    MInstance::~MInstance() {
        if (enableValidationLayers)
            instance->destroyDebugUtilsMessengerEXT(*debugUtilsMessenger, nullptr, vk::DispatchLoaderDynamic());

        instance->destroy();
    }

    bool MInstance::checkValidationLayerSupport() {
        uint32_t layerCount;
        vk::enumerateInstanceLayerProperties(&layerCount, nullptr);

        std::vector<vk::LayerProperties> availableLayers(layerCount);
        vk::enumerateInstanceLayerProperties(&layerCount, availableLayers.data());

        for (auto layerName: validationLayers) {
            bool layerFound = false;

            for (const auto &layerProperties: validationLayers) {
                if (strcmp(layerProperties, layerName) == 0) {
                    layerFound = true;
                    break;
                }
            }

            if (!layerFound)
                return false;
        }

        return true;
    }

    std::vector<const char *> MInstance::getRequiredExtensions() {
        uint32_t extensionCount = 0;
        const char **p_ext = glfwGetRequiredInstanceExtensions(&extensionCount);
        std::vector<const char *> extensions(p_ext, p_ext + extensionCount);

        if (enableValidationLayers)
            extensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);

        return extensions;
    }

    vk::Instance MInstance::getHandle() {
        return *instance;
    }

    void MInstance::setupDebugMessenger() {
        if (!enableValidationLayers) return;

        pfnVkCreateDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>( instance->getProcAddr( "vkCreateDebugUtilsMessengerEXT" ) );
        if ( !pfnVkCreateDebugUtilsMessengerEXT )
        {
            std::cout << "GetInstanceProcAddr: Unable to find pfnVkCreateDebugUtilsMessengerEXT function." << std::endl;
            exit( 1 );
        }

        pfnVkDestroyDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>( instance->getProcAddr( "vkDestroyDebugUtilsMessengerEXT" ) );
        if ( !pfnVkDestroyDebugUtilsMessengerEXT )
        {
            std::cout << "GetInstanceProcAddr: Unable to find pfnVkDestroyDebugUtilsMessengerEXT function." << std::endl;
            exit( 1 );
        }

        vk::DebugUtilsMessengerCreateInfoEXT debugCreateInfo(
                {},
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eError,
                vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
                vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
                vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance,
                debugCallback
        );

        //vk::DebugUtilsMessengerEXT debugUtilsMessenger = instance->createDebugUtilsMessengerEXT(debugCreateInfo);
        //instance->createDebugUtilsMessengerEXT(debugCreateInfo, nullptr, vk::DispatchLoaderDynamic(*instance, ))
        debugUtilsMessenger = std::make_unique<vk::DebugUtilsMessengerEXT>(instance->createDebugUtilsMessengerEXT(debugCreateInfo, nullptr, vk::DispatchLoaderDynamic()));

        //instance->createDebugUtilsMessengerEXTUnique(debugCreateInfo);
    }

    std::vector<const char *> MInstance::getValidationLayers() {
        return validationLayers;
    }
}