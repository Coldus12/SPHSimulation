#ifndef SPHSIMULATION_MSEMAPHORE_HPP
#define SPHSIMULATION_MSEMAPHORE_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MSemaphore {
    public:
        MSemaphore(vk::Device dev, vk::SemaphoreCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::Semaphore>(dev.createSemaphore(info));
        }

        MSemaphore(MSemaphore&& sem)  noexcept {
            handle = std::move(sem.handle);
            this->dev = sem.dev;
        }

        ~MSemaphore() {
            dev.destroySemaphore(*handle);
        }

        vk::Semaphore* getAddress() {
            return handle.get();
        }

        vk::Semaphore getHandle() {
            return *handle;
        }
    private:
        std::unique_ptr<vk::Semaphore> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MSEMAPHORE_HPP
