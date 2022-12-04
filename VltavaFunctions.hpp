#ifndef SPHSIMULATION_VLTAVAFUNCTIONS_HPP
#define SPHSIMULATION_VLTAVAFUNCTIONS_HPP

#include <iostream>
#include <fstream>
#include <chrono>
#include "vulkan/vulkan.hpp"
#include "Managed/Managed.hpp"

namespace Vltava {
    class VulkanResources {
    public:
        static VulkanResources& getInstance() {
            static VulkanResources instance;
            return instance;
        }

        VulkanResources(VulkanResources const&) = delete;
        void operator=(VulkanResources const&) = delete;

        std::unique_ptr<MInstance> instance;
        std::unique_ptr<MSurface> surface;
        std::unique_ptr<MPhysDev> physDev;
        std::unique_ptr<MLogDev> logDev;
        std::unique_ptr<MSwapChain> swapChain;
        std::unique_ptr<vk::RenderPass> renderPass;
        std::unique_ptr<vk::CommandPool> graphicalCmdPool;
        std::unique_ptr<vk::CommandPool> computeCmdPool;

        std::unique_ptr<vk::Queue> graphicsQueue;
        std::unique_ptr<vk::Queue> presentQueue;
        std::unique_ptr<vk::Queue> computeQueue;

        vk::Extent2D extent;
        int FRAMES_IN_FLIGHT = 0;
    private:
        VulkanResources() = default;
    };

    class CPUTimeQuery {
    public:
        void start() {
            start_time = std::chrono::system_clock::now();
        }

        void stop() {
            stop_time = std::chrono::system_clock::now();
        }

        int getElapsedTimeInMillis() {
            std::chrono::duration<float> elapsed_millis = stop_time-start_time;
            return std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_millis).count();
        }

    private:
        std::chrono::time_point<std::chrono::system_clock> start_time, stop_time;
    };

    class VulkanTimeQuery {
    public:
        explicit VulkanTimeQuery(int queryCount) {
            vk::QueryPoolCreateInfo queryInfo;
            queryInfo.queryType = vk::QueryType::eTimestamp;
            queryInfo.queryCount = queryCount;

            pool = VulkanResources::getInstance().logDev->getHandle().createQueryPool(queryInfo);
        }

        vk::QueryPool pool;

        float getTime(int stampNr1, int stampNr2) {
            std::array<uint64_t, 2> timestamps = {0,0};
            vk::Result res = VulkanResources::getInstance().logDev->getHandle().getQueryPoolResults(pool, stampNr1, stampNr2+1, sizeof(timestamps), &timestamps[stampNr1], sizeof(uint64_t), vk::QueryResultFlagBits::e64);
            return ((float) (timestamps[1] - timestamps[0]))/1000000.F;
        }

        ~VulkanTimeQuery() {
            cleanup();
        }
    private:
        VulkanTimeQuery();
        void cleanup() {
            VulkanResources::getInstance().logDev->getHandle().destroyQueryPool(pool);
        }
    };

    uint32_t findMemoryType(uint32_t typeFilter, vk::MemoryPropertyFlags flags);
    std::vector<char> readFile(const std::string &filename);
    void sanityCheck(const std::string &str);
}

#endif //SPHSIMULATION_VLTAVAFUNCTIONS_HPP
