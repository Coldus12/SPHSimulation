#ifndef SPHSIMULATION_MSHADERMODULE_HPP
#define SPHSIMULATION_MSHADERMODULE_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MShaderModule {
    public:
        MShaderModule(vk::Device dev, vk::ShaderModuleCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::ShaderModule>(dev.createShaderModule(info));
        }

        ~MShaderModule() {
            dev.destroyShaderModule(*handle);
        }

        vk::ShaderModule getHandle() {
            return *handle;
        }

        vk::ShaderModule* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::ShaderModule> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MSHADERMODULE_HPP
