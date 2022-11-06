#include "SPH.h"

#define PI 3.1415926538

namespace Vltava {
    // CPU kernel / gradKernel
    //------------------------------------------------------------------------------------------------------------------
    float SPH::kernel(int i, int j) {
        auto& particles = first ? particles1 : particles2;

        float q = glm::length(particles[i].x - particles[j].x)/props.kernelh;
        float m_k = 8.0 / (PI * pow(props.kernelh, 3));
        float m_l = 1.0/(pow(props.kernelh, 3)) * 3.0/(2.0*PI);

        float ret = 0;
        if (q <= 1.0) {
            if (q <= 0.5) {
                ret = m_k * (6.0 * (pow(q, 3.0) - pow(q, 2.0)) + 1.0);
            } else {
                ret = m_k * (2.0 * pow(1 - q, 3));
            }
        }

        return ret;
    }

    // One of the errors: vec3 i == vec3 j ---> normalize(i-j) -> nan
    glm::vec3 SPH::gradKernel(int i, int j) {
        float m_k = 48.0 / (pow(props.kernelh, 3) * PI);
        float m_l = 1.0/(pow(props.kernelh, 3)) * 3.0/(2.0*PI);

        auto& particles = first ? particles1 : particles2;

        glm::vec3 r = particles[i].x - particles[j].x;
        float rlength = length(r);
        float q = rlength / props.kernelh;
        glm::vec3 ret = glm::vec3(0);

        if (q > 0.0001 && q <= 1.0) {
            glm::vec3 gradq = r / (rlength * props.kernelh);

            if (q <= 0.5) {
                gradq *= m_k * q * (3.0 * q - 2);
                ret = gradq;
            } else {
                float factor = 1.0 - q;
                ret = m_k * (-factor * factor) * gradq;
            }
        }

        return ret;
    }

    // Log function
    //------------------------------------------------------------------------------------------------------------------
    std::string SPH::log() {
        std::string str = "";

        if (!logCpuData) {
            auto spheres = sBuffers->at(1).getData<Particle>();
            str += "Compute data - size: " + std::to_string(spheres.size()) + "\n";
            str += "----------------------------------------------\n";

            for (auto & sphere : spheres) {
                str += "Density: " + std::to_string(sphere.rho) +
                       "; Pressure: " + std::to_string(sphere.p) +
                       " ; Position: " + std::to_string(sphere.x.x) + " " + std::to_string(sphere.x.y) + " " + std::to_string(sphere.x.z) +
                       " ; Mass: " + std::to_string(sphere.m) +/* " ; Padding = " +
                       std::to_string(spheres[i].padding) + " ; diff = " +
                       std::to_string(spheres[i].padding - spheres[i].rho) +*/
                       " ; Velocity: " + std::to_string(sphere.v.x) + " " + std::to_string(sphere.v.y) + " " + std::to_string(sphere.v.z) + ";\n";
            }
            str += "\n----------------------------------------------\n";
        } else {
            for (int i = 0; i < particles1.size(); i++) {
                str += std::to_string(i) + "; Density1: " + std::to_string(particles1[i].rho) + " pressure1: " + std::to_string(particles1[i].p) + "\n";
                str += std::to_string(i) + "; Density2: " + std::to_string(particles2[i].rho) + " pressure2: " + std::to_string(particles2[i].p) + "\n\n";
            }
        }

        return str;
    }

    // CPU container
    //------------------------------------------------------------------------------------------------------------------
    void SPH::container(int idx) {
        float left = (props.gridA.x < props.gridB.x ? props.gridA.x : props.gridB.x) + 0.1;
        float right = (props.gridA.x > props.gridB.x ? props.gridA.x : props.gridB.x) - 0.1;

        float front = (props.gridA.y > props.gridB.y ? props.gridA.y : props.gridB.y) - 0.1;
        float back = (props.gridA.y < props.gridB.y ? props.gridA.y : props.gridB.y) + 0.1;

        float bottom = (props.gridA.z < props.gridB.z ? props.gridA.z : props.gridB.z) + 0.1;
        //float bottom = -0.1;
        auto& particles = first ? particles2 : particles1;

        // Bottom
        if (particles[idx].x.z <= bottom) {
            particles[idx].v.z = -particles[idx].v.z;
            particles[idx].x.z = bottom;
        }

        // Left
        if (particles[idx].x.x <= left) {
            particles[idx].v.x = -particles[idx].v.x;
            particles[idx].x.x = left;
        }

        // Right
        if (particles[idx].x.x >= right) {
            particles[idx].v.x = -particles[idx].v.x;
            particles[idx].x.x = right;
        }

        // Front
        if (particles[idx].x.y >= front) {
            particles[idx].v.y = -particles[idx].v.y;
            particles[idx].x.y = front;
        }

        // Back
        if (particles[idx].x.y <= back) {
            particles[idx].v.y = -particles[idx].v.y;
            particles[idx].x.y = back;
        }
    }

    // CPU neighbourhood stuff
    //------------------------------------------------------------------------------------------------------------------
    glm::vec3 SPH::determineGridTuple(int particleIdx) {
        auto& particles = first ? particles1 : particles2;
        auto& p = particles.at(particleIdx);

        // Checking bounds
        glm::vec3 lower(
                props.gridA.x < props.gridB.x ? props.gridA.x : props.gridB.x,
                props.gridA.y < props.gridB.y ? props.gridA.y : props.gridB.y,
                props.gridA.z < props.gridB.z ? props.gridA.z : props.gridB.z
        );

        glm::vec3 upper(
                props.gridA.x > props.gridB.x ? props.gridA.x : props.gridB.x,
                props.gridA.y > props.gridB.y ? props.gridA.y : props.gridB.y,
                props.gridA.z > props.gridB.z ? props.gridA.z : props.gridB.z
        );

        if ((p.x.x < lower.x) || (p.x.y < lower.y) || (p.x.z < lower.z) ||
            (p.x.x > upper.x) || (p.x.y > upper.y) || (p.x.z > upper.z)) {
            return {-1,-1,-1};
        }

        glm::vec3 diff = glm::vec3(props.gridA.x, props.gridA.y, props.gridA.z) - particles.at(particleIdx).x;
        diff /= props.kernelh;

        return floor(abs(diff));
    }

    int SPH::getStartIdxOfCell(glm::vec3 tuple) {
        if (!checkBounds(tuple)) return -1;

        return int((tuple.x * cellz * celly + tuple.y * cellz + tuple.z)*list_size);
    }

    Neighbourhood SPH::getNeighbouringCells(glm::vec3 cellTuple) {
        Neighbourhood ret;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ret.neighbour[i * 9 + j * 3 + k] = glm::vec3(cellTuple.x + (i-1), cellTuple.y + (j-1), cellTuple.z + (k-1));
                }
            }
        }

        return ret;
    }

    bool SPH::checkBounds(glm::vec3 tuple) {
        if (tuple.x < 0 || tuple.x >= cellx)
            return false;

        if (tuple.y < 0 || tuple.y >= celly)
            return false;

        if (tuple.z < 0 || tuple.z >= cellz)
            return false;

        return true;
    }

    void SPH::place() {
        auto& particles = first ? particles1 : particles2;

        grid_data.clear();
        grid_data.resize(1, 0);
        grid_data.resize(list_size*cellx*celly*cellz, 0);

        for (int i = 0; i < particles.size(); i++) {
            placeParticleIntoCell(i);
        }
    }

    void SPH::placeParticleIntoCell(int particleIdx) {
        glm::vec3 cellTuple = determineGridTuple(particleIdx);
#if ppic
        std::cout << "[ppic] particleIdx" << particleIdx << " placed into (" << cellTuple.x << "," << cellTuple.y << "," << cellTuple.z << ")\n";
#endif
        //if (checkBounds(cellTuple)) {
        int startIdx = getStartIdxOfCell(cellTuple);
        if (startIdx >= 0) {
            int realIdx = /*atomicAdd(*/grid_data.at(startIdx);
            grid_data.at(startIdx) += 1;
            realIdx += startIdx + 1;

            // Placing particleIdx into cell
            grid_data.at(realIdx) = particleIdx;
        }
        //}
    }
}