#ifndef SPHSIMULATION_MIMAGEVIEW_HPP
#define SPHSIMULATION_MIMAGEVIEW_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MImageView {
    public:
        MImageView(vk::Device dev, vk::ImageViewCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::ImageView>(dev.createImageView(info));
        }

        MImageView(MImageView&& view)  noexcept {
            handle = std::move(view.handle);
            this->dev = view.dev;
        }

        ~MImageView() {
            dev.destroyImageView(*handle);
        }

        vk::ImageView getHandle() {
            return *handle;
        }

        vk::ImageView* getAddress() {
            return handle.get();
        }
    private:
        std::unique_ptr<vk::ImageView> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MIMAGEVIEW_HPP
