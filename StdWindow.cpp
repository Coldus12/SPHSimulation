#include "StdWindow.hpp"
#include "VltavaFunctions.hpp"
#include "ComputeShader.hpp"
#include "ParticleModel.hpp"
#include "CPUSim.h"

#define IMGUI_ENABLED true

// ImGui log - taken from imgui/showDemoWindow.cpp
struct ImLog
{
    ImGuiTextBuffer buff;
    ImVector<int> lineOffsets; // Index to lines offset. We maintain this with AddLog() calls.
    bool autoScroll;  // Keep scrolling if already at the bottom.

    ImLog() {
        autoScroll = true;
        clear();
    }

    void clear() {
        buff.clear();
        lineOffsets.clear();
        lineOffsets.push_back(0);
    }

    void addLog(const char* fmt, ...) IM_FMTARGS(2)
    {
        int old_size = buff.size();
        va_list args;
        va_start(args, fmt);
        buff.appendfv(fmt, args);
        va_end(args);
        for (int new_size = buff.size(); old_size < new_size; old_size++)
            if (buff[old_size] == '\n')
                lineOffsets.push_back(old_size + 1);
    }

    void draw(const char* title, bool* p_open = NULL) {
        if (!ImGui::Begin(title, p_open)) {
            ImGui::End();
            return;
        }

        // Main window
        if (ImGui::Button("Options"))
            ImGui::OpenPopup("Options");
        ImGui::SameLine();
        bool clear = ImGui::Button("Clear");
        ImGui::SameLine();
        bool copy = ImGui::Button("Copy");
        ImGui::SameLine();

        ImGui::Separator();
        ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

        if (clear)
            this->clear();
        if (copy)
            ImGui::LogToClipboard();

        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 0));
        const char* buf = buff.begin();
        const char* buf_end = buff.end();
        ImGuiListClipper clipper;
        clipper.Begin(lineOffsets.Size);
        while (clipper.Step()) {
            for (int line_no = clipper.DisplayStart; line_no < clipper.DisplayEnd; line_no++) {
                const char* line_start = buf + lineOffsets[line_no];
                const char* line_end = (line_no + 1 < lineOffsets.Size) ? (buf + lineOffsets[line_no + 1] - 1) : buf_end;
                ImGui::TextUnformatted(line_start, line_end);
            }
        }
        clipper.End();
        ImGui::PopStyleVar();

        if (autoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
            ImGui::SetScrollHereY(1.0f);

        ImGui::EndChild();
        ImGui::End();
    }
};

static ImLog my_log;

namespace Vltava {
    bool StdWindow::run = false;
    bool StdWindow::rot = true;

    // Main loop
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::mainloop() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
#if IMGUI_ENABLED
            ImGui_ImplVulkan_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            //ImGui::ShowDemoWindow();
            ImGui::Begin("Test");

            ImGui::Text("Reset data??");
            if (ImGui::Button("Yes"))
                resetData();

            ImGui::InputInt("Number of iterations", &nrOfIter);
            ImGui::InputInt("Number of particles: ", &particleNr);

            ImGui::Checkbox("Show log?", &show_log);
            ImGui::Checkbox("Write log?", &write_log);
            ImGui::Checkbox("Write log into cli?", &console_log);

            if (show_log)
                my_log.draw("Sup");

            ImGui::End();
#endif
            runComp();
            drawFrame();
        }

        VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    void StdWindow::runComp() {
        if (run) {
            dispatchCompute(nrOfP, 1, 1);

            if (write_log) {
                auto spheres = sBuffers[1].getData<Particle>();
                std::string str = "";
                str += "Compute data - size: " + std::to_string(spheres.size()) + "\n";
                str += "----------------------------------------------\n";

                for (int i = 0; i < nrOfP; i++) {
                    str += "Density: " + std::to_string(spheres[i].rho) +
                           "; Pressure: " + std::to_string(spheres[i].p) +
                           " ; Position: " + std::to_string(spheres[i].x.x) + " " + std::to_string(spheres[i].x.y) +
                           " " + std::to_string(spheres[i].x.z) +
                           " ; Mass: " + std::to_string(spheres[i].m) +
                           " ; Velocity: " + std::to_string(spheres[i].v.x) + " " + std::to_string(spheres[i].v.y) +
                           " " + std::to_string(spheres[i].v.z) + ";\n";
                }
                str += "\n----------------------------------------------\n";

                my_log.addLog(str.c_str());

                if (console_log)
                    std::cout << str << std::endl;
            }

            run = false;
        }
    }

    void StdWindow::resetData() {
        comp1.reset();
        comp2.reset();
        sBuffers.clear();
        uBuffers.clear();
        initCompute();
        model->changeModel(all_particle_nr, &uBuffers, &sBuffers);
    }

    // Try changing data to std::vec, and only allocate gpu memory after you are done filling said vector.
    void StdWindow::setComputeData() {
        // Water rest density = 1 kg/m^3
        //
        // Particle h = 0.05 meter
        // Particle volume = 5.2359 * 10^-4 m^3 = 0.005236 m^3 (4/3 * 0.05^3 * PI)
        // Particle mass = volume * density = 0.005236 kg
        //
        // Smoothing length := 3 * 2 * h = 0.1 * 3 meter (for now)
        int pnrAlongAxis = round(pow(particleNr, 1.0/3.0));
        std::vector<Particle> particles;

        nrOfP = pnrAlongAxis * pnrAlongAxis * pnrAlongAxis;
        particleNr = nrOfP;

        // Solid box
        int bsize = 16;

        float s=0.1;
        //vk::DeviceSize size = sizeof(Particle) * nrOfP * bsize * bsize;
        //auto* data = new Particle[nrOfP + bsize * bsize];

        int r = -1;
        int z = -1;

        float mass = 0.005236f;

        // Actual simulated particles
        //--------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < nrOfP; i++) {
            Particle data;

            if (i % pnrAlongAxis == 0)
                r++;

            if (i % (pnrAlongAxis * pnrAlongAxis) == 0)
                z++;

            data.x = glm::vec3((i%pnrAlongAxis) * s,(r%pnrAlongAxis) * s, z * s);
            data.h = 1;
            data.v = glm::vec3(0,0,0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 0;
            data.padding = 0;

            particles.push_back(data);
        }

        float dist = 0.05;
        r = -1;
        z = -1;

        // Container
        //--------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < bsize*bsize; i++) {
            if (i%bsize == 0)
                r++;

            Particle data;

            // Bottom
            data.x = glm::vec3((i%bsize) * dist,(r%bsize) * dist, 0);
            data.h = 1;
            data.v = glm::vec3(0,0,0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 1;
            data.padding = 0;

            particles.push_back(data);

            // Wall1
            data.x = glm::vec3((i%bsize) * dist,-dist*0, (r%bsize) * dist);
            particles.push_back(data);

            // Wall2
            data.x = glm::vec3(-dist*0,(i%bsize) * dist, (r%bsize) * dist);
            particles.push_back(data);

            // Wall3
            data.x = glm::vec3((i%bsize) * dist, (bsize+1*0) * dist, (bsize - r%bsize) * dist);
            particles.push_back(data);

            // Wall4
            data.x = glm::vec3((bsize+1*0) * dist, ( i%bsize) * dist, (r%bsize) * dist);
            particles.push_back(data);
        }

        all_particle_nr = particles.size();

        vk::DeviceSize size = sizeof(Particle) * all_particle_nr;

        SimProps props{
                1.0f,
                0.1f,
                (float) all_particle_nr,
                0.1f
        };
        uBuffers[0].writeToBuffer(&props, sizeof(props));

        Buffer inBuffer(
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        inBuffer.bind(0);

        Buffer outBuffer(
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        outBuffer.setSize(all_particle_nr * sizeof(Particle));
        outBuffer.bind(0);

        sBuffers.push_back(std::move(inBuffer));
        sBuffers.push_back(std::move(outBuffer));

        for (auto& b: sBuffers)
            b.writeToBuffer(particles.data(), size);

        //delete[] data;
    }

    void StdWindow::initCompute() {
        Buffer UBO(
                sizeof(SimProps),
                vk::BufferUsageFlagBits::eUniformBuffer,
                vk::MemoryPropertyFlagBits::eHostCoherent | vk::MemoryPropertyFlagBits::eHostVisible
        );
        uBuffers.push_back(std::move(UBO));

        setComputeData();

        comp1 = std::make_unique<ComputeShader>("shaders/comp.spv");
        comp1->setBuffers(&uBuffers, &sBuffers);
        comp1->createPipeline();

        comp2 = std::make_unique<ComputeShader>("shaders/comp_it.spv");
        comp2->setBuffers(&uBuffers, &sBuffers);
        comp2->createPipeline();
    }

    //------------------------------------------------------------------------------------------------------------------

    void StdWindow::dispatchCompute(int groupCountX, int groupCountY, int groupCountZ) {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        computeCmdBuffer.begin(beginInfo);

        std::vector<vk::BufferMemoryBarrier> membarriers;
        for (auto& buffer: sBuffers) {
            membarriers.emplace_back(
                    vk::AccessFlagBits::eShaderWrite,
                    vk::AccessFlagBits::eShaderRead,
                    computeQueueFamily,
                    computeQueueFamily,
                    buffer.getBufferHandle(),
                    0,
                    VK_WHOLE_SIZE
            );
        }

        for (int i = 0; i < nrOfIter; i++) {
            comp1->bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(groupCountX, groupCountY, groupCountZ);
            computeCmdBuffer.pipelineBarrier(
                    vk::PipelineStageFlagBits::eComputeShader,
                    vk::PipelineStageFlagBits::eComputeShader,
                    {},
                    {},
                    membarriers,
                    {}
            );

            comp2->bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(groupCountX, groupCountY, groupCountZ);
            computeCmdBuffer.pipelineBarrier(
                    vk::PipelineStageFlagBits::eComputeShader,
                    vk::PipelineStageFlagBits::eComputeShader,
                    {},
                    {},
                    membarriers,
                    {}
            );
        }

        computeCmdBuffer.end();

        vk::SubmitInfo submitInfo(
                {},
                {},
                computeCmdBuffer,
                {}
        );

        VulkanResources::getInstance().computeQueue->submit(submitInfo);
        //VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    //------------------------------------------------------------------------------------------------------------------

    void StdWindow::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_W && action == GLFW_PRESS) {
            run = true;
        }

        if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
            rot = !rot;
        }
    }

    // Constructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::StdWindow(int width, int height) : width(width), height(height) {
        VulkanResources::getInstance().FRAMES_IN_FLIGHT = 2;

        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        //glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

        window = glfwCreateWindow(width, height, "Window", nullptr, nullptr);
        glfwSetWindowUserPointer(window, this);
        glfwSetFramebufferSizeCallback(window, frameBufferResizeCallback);

        glfwSetKeyCallback(window, key_callback);

        initVulkan();
#if IMGUI_ENABLED
        initImgui();
#endif
        mainloop();
    }

    void StdWindow::initVulkan() {
        createInstance("Test");
        createSurface();
        selectPhysicalDevice();
        selectQueues();
        createLogicalDevice();
        createSwapChain();
        createImageViews();
        createRenderPass();
        createFramebuffers();
        createCommandPool();
        createCommandBuffers();
        createSyncObjects();

        initCompute();

        loadModel();
    }

    // Cleaning up
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::cleanup() {
#if IMGUI_ENABLED
        vkDestroyDescriptorPool(
                VulkanResources::getInstance().logDev->getHandle(),
                imguiPool,
                nullptr
        );
        ImGui_ImplVulkan_Shutdown();
#endif

        cleanupSwapChain();

        for (int i = 0; i < VulkanResources::getInstance().FRAMES_IN_FLIGHT; i++) {
            VulkanResources::getInstance().logDev->getHandle().destroySemaphore(renderFinishedSemaphores[i]);
            VulkanResources::getInstance().logDev->getHandle().destroySemaphore(imageAvailableSemaphores[i]);
            VulkanResources::getInstance().logDev->getHandle().destroyFence(inFlightFences[i]);
        }

        VulkanResources::getInstance().logDev->getHandle().destroyCommandPool(*VulkanResources::getInstance().graphicalCmdPool);
        VulkanResources::getInstance().logDev->getHandle().destroyCommandPool(*VulkanResources::getInstance().computeCmdPool);

        //glfwDestroyWindow(window);
        //glfwTerminate();
    }

    // Destructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::~StdWindow() {
        cleanup();
    }

    // DearImgui initialization
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::initImgui() {
        //1: create descriptor pool for IMGUI
        // the size of the pool is very oversize, but it's copied from imgui demo itself.
        VkDescriptorPoolSize pool_sizes[] =
                {
                        { VK_DESCRIPTOR_TYPE_SAMPLER, 1000 },
                        { VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1000 },
                        { VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_TEXEL_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_TEXEL_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 1000 },
                        { VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 1000 }
                };

        VkDescriptorPoolCreateInfo pool_info = {};
        pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
        pool_info.maxSets = 1000;
        pool_info.poolSizeCount = std::size(pool_sizes);
        pool_info.pPoolSizes = pool_sizes;

        vkCreateDescriptorPool(
                VulkanResources::getInstance().logDev->getHandle(),
                &pool_info,
                nullptr,
                &imguiPool
        );

        // 2: initialize imgui library

        //this initializes the core structures of imgui
        ImGui::CreateContext();

        //this initializes imgui for GLFW
        ImGui_ImplGlfw_InitForVulkan(window, true);

        //this initializes imgui for Vulkan
        ImGui_ImplVulkan_InitInfo init_info = {};
        init_info.Instance = VulkanResources::getInstance().instance->getHandle();
        init_info.PhysicalDevice = VulkanResources::getInstance().physDev->getHandle();
        init_info.Device = VulkanResources::getInstance().logDev->getHandle();
        init_info.Queue = *VulkanResources::getInstance().graphicsQueue;
        init_info.DescriptorPool = imguiPool;
        init_info.MinImageCount = 3;
        init_info.ImageCount = 3;
        init_info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;

        ImGui_ImplVulkan_Init(&init_info, *VulkanResources::getInstance().renderPass);

        //execute a gpu command to upload imgui font textures
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        commandBuffers[currentFrame].begin(beginInfo);
        ImGui_ImplVulkan_CreateFontsTexture(commandBuffers[currentFrame]);
        commandBuffers[currentFrame].end();

        vk::SubmitInfo submitInfo({},{}, commandBuffers[currentFrame]);
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        VulkanResources::getInstance().graphicsQueue->submit(infos, inFlightFences[currentFrame]);
        VulkanResources::getInstance().logDev->getHandle().waitIdle();

        //clear font textures from cpu data
        ImGui_ImplVulkan_DestroyFontUploadObjects();
    }

    // Instance creation
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createInstance(std::string app_name) {
        VulkanResources::getInstance().instance = std::make_unique<MInstance>(
                app_name,
                enableValidationLayers
        );
    }

    // Creating a surface
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSurface() {
        VulkanResources::getInstance().surface = std::make_unique<MSurface>(
                VulkanResources::getInstance().instance->getHandle(),
                window
        );
    }

    // Selecting the best gpu for the task
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectPhysicalDevice() {
        VulkanResources::getInstance().physDev = std::make_unique<MPhysDev>(
                VulkanResources::getInstance().instance->getHandle(),
                deviceExtensions
        );
    }

    // Selecting queue(s)
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectQueues() {
        std::vector<vk::QueueFamilyProperties> queueFamilies = VulkanResources::getInstance().physDev->getHandle().getQueueFamilyProperties();

        uint32_t i = 0;
        for (const auto &queueFamily: queueFamilies) {
            if (VulkanResources::getInstance().physDev->getHandle().getSurfaceSupportKHR(i, VulkanResources::getInstance().surface->getHandle())) {
                presentQueueFamily = i;
            }

            if (queueFamily.queueFlags & vk::QueueFlagBits::eGraphics) {
                graphicsQueueFamily = i;
            }

            if (queueFamily.queueFlags & vk::QueueFlagBits::eCompute) {
                computeQueueFamily = i;
            }
            i++;
        }
    }

    // Creating a logical device
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createLogicalDevice() {
        std::set<uint32_t> uniqueQueueFamilies({graphicsQueueFamily, presentQueueFamily, computeQueueFamily});

        VulkanResources::getInstance().logDev = std::make_unique<MLogDev>(
                uniqueQueueFamilies,
                VulkanResources::getInstance().physDev->getHandle(),
                deviceExtensions,
                &validationLayers
        );

        // Queue creation
        //--------------------------------------------------------------------------------------------------------------
        //Clearing previous queues
        VulkanResources::getInstance().graphicsQueue.reset();
        VulkanResources::getInstance().presentQueue.reset();
        VulkanResources::getInstance().computeQueue.reset();

        VulkanResources::getInstance().graphicsQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(graphicsQueueFamily, 0));
        VulkanResources::getInstance().presentQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(presentQueueFamily, 0));
        VulkanResources::getInstance().computeQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(computeQueueFamily, 0));
    }

    // SwapChain
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSwapChain() {
        VulkanResources::getInstance().swapChain = std::make_unique<MSwapChain>(
                window,
                VulkanResources::getInstance().physDev->getHandle(),
                VulkanResources::getInstance().logDev->getHandle(),
                VulkanResources::getInstance().surface->getHandle(),
                vk::Format::eB8G8R8A8Srgb,
                vk::ColorSpaceKHR::eSrgbNonlinear,
                vk::PresentModeKHR::eMailbox,
                graphicsQueueFamily,
                presentQueueFamily
        );

        swapChainImages = VulkanResources::getInstance().logDev->getHandle().getSwapchainImagesKHR(VulkanResources::getInstance().swapChain->getHandle());
        imageViews.reserve(swapChainImages.size());

        swapChainImageFormat = VulkanResources::getInstance().swapChain->getFormat();
        VulkanResources::getInstance().extent = VulkanResources::getInstance().swapChain->getExtent();
    }

    void StdWindow::recreateSwapChain() {
        int w = 0, h = 0;
        glfwGetFramebufferSize(window, &w, &h);
        while (w == 0 || h == 0) {
            glfwGetFramebufferSize(window, &w, &h);
            glfwWaitEvents();
        }

        VulkanResources::getInstance().logDev->getHandle().waitIdle();
        cleanupSwapChain();
        createSwapChain();
        createImageViews();
        createRenderPass();
        model->recreatePipeline();
        createFramebuffers();
    }

    void StdWindow::cleanupSwapChain() {
        for (auto framebuffer : swapChainFramebuffers) {
            VulkanResources::getInstance().logDev->getHandle().destroyFramebuffer(framebuffer);
        }
        swapChainFramebuffers.clear();

        VulkanResources::getInstance().logDev->getHandle().destroyRenderPass(*VulkanResources::getInstance().renderPass);

        for (auto imgView: imageViews) {
            VulkanResources::getInstance().logDev->getHandle().destroyImageView(imgView);
        }
        imageViews.clear();

        VulkanResources::getInstance().swapChain.reset();
    }

    // Framebuffer resize callback function
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::frameBufferResizeCallback(GLFWwindow *window, int width, int height) {
        auto app = reinterpret_cast<StdWindow *>(glfwGetWindowUserPointer(window));
        app->framebufferResized = true;
    }

    // ImageView
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createImageViews() {
        imageViews.reserve(swapChainImages.size());

        vk::ImageViewCreateInfo createInfo(
                {},
                {},
                vk::ImageViewType::e2D,
                swapChainImageFormat,
                {vk::ComponentSwizzle::eIdentity, vk::ComponentSwizzle::eIdentity, vk::ComponentSwizzle::eIdentity,
                 vk::ComponentSwizzle::eIdentity},
                {vk::ImageAspectFlagBits::eColor, 0, 1, 0, 1}
        );

        for (auto image: swapChainImages) {
            createInfo.image = image;
            imageViews.push_back(VulkanResources::getInstance().logDev->getHandle().createImageView(createInfo));
        }
    }

    // RenderPass
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createRenderPass() {
        vk::SubpassDependency dependency(
                VK_SUBPASS_EXTERNAL,
                0,
                vk::PipelineStageFlagBits::eColorAttachmentOutput,
                vk::PipelineStageFlagBits::eColorAttachmentOutput,
                {},
                vk::AccessFlagBits::eColorAttachmentWrite
        );

        vk::AttachmentDescription colorAttachment(
                {},
                swapChainImageFormat,
                vk::SampleCountFlagBits::e1,
                vk::AttachmentLoadOp::eClear,
                vk::AttachmentStoreOp::eStore,
                vk::AttachmentLoadOp::eDontCare,
                vk::AttachmentStoreOp::eDontCare,
                vk::ImageLayout::eUndefined,
                vk::ImageLayout::ePresentSrcKHR
        );

        vk::AttachmentReference colorReference(0, vk::ImageLayout::eColorAttachmentOptimal);

        vk::SubpassDescription subpass({}, vk::PipelineBindPoint::eGraphics, {}, colorReference);
        vk::RenderPassCreateInfo createInfo(
                {},
                1,
                &colorAttachment,
                1,
                &subpass,
                1,
                &dependency
        );

        VulkanResources::getInstance().renderPass = std::make_unique<vk::RenderPass>(
                VulkanResources::getInstance().logDev->getHandle().createRenderPass(createInfo)
        );
    }

    // Framebuffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createFramebuffers() {
        swapChainFramebuffers.reserve(imageViews.size());

        for (int i = 0; i < imageViews.size(); i++) {
            vk::FramebufferCreateInfo framebufferInfo(
                    {},
                    *VulkanResources::getInstance().renderPass,
                    1,
                    &imageViews[i],
                    VulkanResources::getInstance().extent.width,
                    VulkanResources::getInstance().extent.height,
                    1
            );

            swapChainFramebuffers.push_back(VulkanResources::getInstance().logDev->getHandle().createFramebuffer(framebufferInfo));
        }
    }

    // Commandpool
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandPool() {
        vk::CommandPoolCreateInfo poolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueFamily);
        VulkanResources::getInstance().graphicalCmdPool = std::make_unique<vk::CommandPool>(
                VulkanResources::getInstance().logDev->getHandle().createCommandPool(poolInfo)
        );

        vk::CommandPoolCreateInfo cmdPoolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, computeQueueFamily);
        VulkanResources::getInstance().computeCmdPool = std::make_unique<vk::CommandPool>(
                VulkanResources::getInstance().logDev->getHandle().createCommandPool(cmdPoolInfo)
        );
    }

    // Command buffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandBuffers() {
        vk::CommandBufferAllocateInfo allocInfo(*VulkanResources::getInstance().graphicalCmdPool, vk::CommandBufferLevel::ePrimary, VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        commandBuffers = VulkanResources::getInstance().logDev->getHandle().allocateCommandBuffers(allocInfo);

        vk::CommandBufferAllocateInfo cmdBufferInfo(*VulkanResources::getInstance().computeCmdPool, vk::CommandBufferLevel::ePrimary, 1);
        computeCmdBuffer = VulkanResources::getInstance().logDev->getHandle().allocateCommandBuffers(cmdBufferInfo).front();
    }

    // Sync objects
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSyncObjects() {
        imageAvailableSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        renderFinishedSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        inFlightFences.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);

        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        vk::SemaphoreCreateInfo semaphoreCreateInfo;

        for (int i = 0; i < VulkanResources::getInstance().FRAMES_IN_FLIGHT; i++) {
            imageAvailableSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            renderFinishedSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            inFlightFences.push_back(VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo));
        }
    }

    // Load model
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::loadModel() {
        model = std::make_unique<ParticleModel>(&uBuffers, &sBuffers);
    }

    // Draw frame
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::drawFrame() {
#if IMGUI_ENABLED
        ImGui::Render();
#endif
        //VulkanResources::getInstance().logDev->getHandle().waitForFences(inFlightFences[currentFrame], true, UINT64_MAX);
        VulkanResources::getInstance().logDev->getHandle().waitForFences(inFlightFences[currentFrame], true, 10);

        vk::Result result;
        uint32_t imageIndex;

        std::tie(result, imageIndex) = VulkanResources::getInstance().logDev->getHandle().acquireNextImageKHR(
                VulkanResources::getInstance().swapChain->getHandle(),
                UINT64_MAX,
                imageAvailableSemaphores[currentFrame]
        );

        if (result == vk::Result::eErrorOutOfDateKHR) {
            recreateSwapChain();
            return;
        } else if (result != vk::Result::eSuccess && result != vk::Result::eSuboptimalKHR) {
            throw std::runtime_error("Failed to acquire swap chain image!");
        }

        VulkanResources::getInstance().logDev->getHandle().resetFences(inFlightFences[currentFrame]);

        commandBuffers[currentFrame].reset();
        recordCommandBuffer(imageIndex);

        vk::PipelineStageFlags flags(vk::PipelineStageFlagBits::eColorAttachmentOutput);
        vk::SubmitInfo submitInfo(
                (uint32_t) 1,
                &imageAvailableSemaphores[currentFrame],
                &flags,
                (uint32_t) 1,
                &commandBuffers[currentFrame],
                (uint32_t) 1,
                &renderFinishedSemaphores[currentFrame]
        );
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        VulkanResources::getInstance().graphicsQueue->submit(infos, inFlightFences[currentFrame]);

        std::array<vk::SwapchainKHR, 1> swapChains = {VulkanResources::getInstance().swapChain->getHandle()};
        vk::PresentInfoKHR presentInfo(renderFinishedSemaphores[currentFrame], swapChains, imageIndex);

        /** For some reason instead of returning the vk::eErrorOutOfDateKHR result,
         * this just throws an exception with the this error. But this try/catch block
         * seems to fix the problem.
         * */
        try {
            result = VulkanResources::getInstance().presentQueue->presentKHR(presentInfo);
        } catch (vk::SystemError &err) {
            framebufferResized = false;
            recreateSwapChain();
        }

        if (result == vk::Result::eErrorOutOfDateKHR || result == vk::Result::eSuboptimalKHR || framebufferResized) {
            framebufferResized = false;
            recreateSwapChain();
        } else if (result != vk::Result::eSuccess) {
            throw std::runtime_error("Failed to present swap chain image!");
        }

        currentFrame = (currentFrame + 1) % VulkanResources::getInstance().FRAMES_IN_FLIGHT;
    }

    void StdWindow::recordCommandBuffer(uint32_t imageIndex) {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        commandBuffers[currentFrame].begin(beginInfo);

        std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
        vk::ClearValue clrVal((vk::ClearColorValue(color)));

        vk::RenderPassBeginInfo renderPassInfo(
                *VulkanResources::getInstance().renderPass,
                swapChainFramebuffers[imageIndex],
                {{0, 0}, VulkanResources::getInstance().extent},
                1, &clrVal
        );

        commandBuffers[currentFrame].beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);

        // Drawing the models
        model->draw(commandBuffers[currentFrame],currentFrame);
        //model->draw(commandBuffers->at(currentFrame));

#if IMGUI_ENABLED
        // Recording imgui draw cmds into the commandBuffer
        ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), commandBuffers[currentFrame]);
#endif

        commandBuffers[currentFrame].endRenderPass();
        commandBuffers[currentFrame].end();
    }
}