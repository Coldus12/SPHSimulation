#include "MInstance.hpp"
#include "../VltavaFunctions.hpp"

/** Taken from: https://github.com/KhronosGroup/Vulkan-Hpp/blob/master/samples/CreateDebugUtilsMessenger/CreateDebugUtilsMessenger.cpp
 * */
PFN_vkCreateDebugUtilsMessengerEXT  pfnVkCreateDebugUtilsMessengerEXT;
PFN_vkDestroyDebugUtilsMessengerEXT pfnVkDestroyDebugUtilsMessengerEXT;

VKAPI_ATTR VkResult VKAPI_CALL vkCreateDebugUtilsMessengerEXT( VkInstance                                 instance,
                                                               const VkDebugUtilsMessengerCreateInfoEXT * pCreateInfo,
                                                               const VkAllocationCallbacks *              pAllocator,
                                                               VkDebugUtilsMessengerEXT *                 pMessenger )
{
    return pfnVkCreateDebugUtilsMessengerEXT( instance, pCreateInfo, pAllocator, pMessenger );
}

VKAPI_ATTR void VKAPI_CALL vkDestroyDebugUtilsMessengerEXT( VkInstance instance, VkDebugUtilsMessengerEXT messenger, VkAllocationCallbacks const * pAllocator )
{
    return pfnVkDestroyDebugUtilsMessengerEXT( instance, messenger, pAllocator );
}

VKAPI_ATTR VkBool32 VKAPI_CALL debugMessageFunc( VkDebugUtilsMessageSeverityFlagBitsEXT       messageSeverity,
                                                 VkDebugUtilsMessageTypeFlagsEXT              messageTypes,
                                                 VkDebugUtilsMessengerCallbackDataEXT const * pCallbackData,
                                                 void * /*pUserData*/ )
{
    std::ostringstream message;

    message << vk::to_string( static_cast<vk::DebugUtilsMessageSeverityFlagBitsEXT>( messageSeverity ) ) << ": "
            << vk::to_string( static_cast<vk::DebugUtilsMessageTypeFlagsEXT>( messageTypes ) ) << ":\n";
    message << "\t"
            << "messageIDName   = <" << pCallbackData->pMessageIdName << ">\n";
    message << "\t"
            << "messageIdNumber = " << pCallbackData->messageIdNumber << "\n";
    message << "\t"
            << "message         = <" << pCallbackData->pMessage << ">\n";
    if ( 0 < pCallbackData->queueLabelCount )
    {
        message << "\t"
                << "Queue Labels:\n";
        for ( uint32_t i = 0; i < pCallbackData->queueLabelCount; i++ )
        {
            message << "\t\t"
                    << "labelName = <" << pCallbackData->pQueueLabels[i].pLabelName << ">\n";
        }
    }
    if ( 0 < pCallbackData->cmdBufLabelCount )
    {
        message << "\t"
                << "CommandBuffer Labels:\n";
        for ( uint32_t i = 0; i < pCallbackData->cmdBufLabelCount; i++ )
        {
            message << "\t\t"
                    << "labelName = <" << pCallbackData->pCmdBufLabels[i].pLabelName << ">\n";
        }
    }
    if ( 0 < pCallbackData->objectCount )
    {
        message << "\t"
                << "Objects:\n";
        for ( uint32_t i = 0; i < pCallbackData->objectCount; i++ )
        {
            message << "\t\t"
                    << "Object " << i << "\n";
            message << "\t\t\t"
                    << "objectType   = " << vk::to_string( static_cast<vk::ObjectType>( pCallbackData->pObjects[i].objectType ) ) << "\n";
            message << "\t\t\t"
                    << "objectHandle = " << pCallbackData->pObjects[i].objectHandle << "\n";
            if ( pCallbackData->pObjects[i].pObjectName )
            {
                message << "\t\t\t"
                        << "objectName   = <" << pCallbackData->pObjects[i].pObjectName << ">\n";
            }
        }
    }

#ifdef _WIN32
    MessageBox( NULL, message.str().c_str(), "Alert", MB_OK );
#else
    std::cout << message.str() << std::endl;
#endif
    return false;
}

namespace Vltava {
    MInstance::MInstance(std::string name, bool enableValidationLayers): enableValidationLayers(enableValidationLayers) {
        vk::ApplicationInfo appInfo(
                name.c_str(),
                VK_MAKE_VERSION(1, 0, 0),
                "No Engine",
                VK_MAKE_VERSION(1, 0, 0),
                VK_API_VERSION_1_1
        );

        vk::InstanceCreateInfo createInfo({}, &appInfo);

        auto glfwExtensions = getRequiredExtensions();

        if (enableValidationLayers) {
            if (!checkValidationLayerSupport())
                throw std::runtime_error("Validation layers requested, but not available!");

            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
            createInfo.ppEnabledLayerNames = validationLayers.data();

            /*vk::DebugUtilsMessengerCreateInfoEXT debugCreateInfo(
                    {},
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eError,
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance,
                    &debugMessageFunc
            );

            createInfo.pNext = (vk::DebugUtilsMessengerCreateInfoEXT *) &debugCreateInfo;*/
        } else {
            createInfo.enabledLayerCount = 0;
            //createInfo.pNext = nullptr;
        }

        createInfo.enabledExtensionCount = static_cast<uint32_t>(glfwExtensions.size());
        createInfo.ppEnabledExtensionNames = glfwExtensions.data();

        createInfo.pNext = nullptr;

        instance = vk::createInstance(createInfo);
        setupDebugMessenger();
    }

    MInstance::~MInstance() {
        if (enableValidationLayers)
            instance.destroyDebugUtilsMessengerEXT(debugUtilsMessenger);

        instance.destroy();
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
        return instance;
    }

    void MInstance::setupDebugMessenger() {
        if (!enableValidationLayers) return;

        pfnVkCreateDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>( instance.getProcAddr( "vkCreateDebugUtilsMessengerEXT" ) );
        if ( !pfnVkCreateDebugUtilsMessengerEXT )
        {
            std::cout << "GetInstanceProcAddr: Unable to find pfnVkCreateDebugUtilsMessengerEXT function." << std::endl;
            exit( 1 );
        }

        pfnVkDestroyDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>( instance.getProcAddr( "vkDestroyDebugUtilsMessengerEXT" ) );
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
                &debugMessageFunc
        );

        debugUtilsMessenger = instance.createDebugUtilsMessengerEXT(debugCreateInfo);
    }

    std::vector<const char *> MInstance::getValidationLayers() {
        return validationLayers;
    }
}