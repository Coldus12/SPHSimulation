#ifndef SPHSIMULATION_IMLOG_HPP
#define SPHSIMULATION_IMLOG_HPP

#include "imgui_impl_glfw.h"
#include "imgui_impl_vulkan.h"

namespace Vltava {

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

} // Vltava

#endif //SPHSIMULATION_IMLOG_HPP
