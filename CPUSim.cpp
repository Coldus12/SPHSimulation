#include "CPUSim.hpp"

#define PI 3.1415926538

namespace Vltava {

// Neighbourhood stuff
//----------------------------------------------------------------------------------------------------------------------
    struct Neighbourhood {
        glm::vec3 neighbour[27]; // 3^3
    };

    /*vec3 determineGridTuple(int particleIdx) {
        vec3 diff = ubo.gridA.xyz - in_data.p[particleIdx].x;
        diff /= ubo.kernelh;

        return floor(abs(diff));
    }

    int getStartIdxOfCell(vec3 tuple) {
        return int((tuple.x * cellz * celly + tuple.y * cellz + tuple.z)*list_size);
    }

    Neighbourhood getNeighbouringCells(vec3 cellTuple) {
        Neighbourhood ret;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ret.neighbour[i * 9 + j * 3 + k] = vec3(cellTuple.x + (i-1), cellTuple.y + (j-1), cellTuple.z + (k-1));
                }
            }
        }

        return ret;
    }

    bool checkBounds(vec3 tuple) {
        if (tuple.x < 0 || tuple.x >= cellx)
            return false;

        if (tuple.y < 0 || tuple.y >= celly)
            return false;

        if (tuple.z < 0 || tuple.z >= cellz)
            return false;

        return true;
    }

    void placeParticleIntoCell(int particleIdx) {
        vec3 cellTuple = determineGridTuple(particleIdx);

        int startIdx = getStartIdxOfCell(cellTuple);
        uint realIdx = atomicAdd(grid_data.numbers[startIdx], 1);
        realIdx += startIdx + 1;

        // Placing particleIdx into cell
        grid_data.numbers[realIdx] = particleIdx;
        //grid_data.numbers[particleIdx] = int(realIdx);
    }*/
//----------------------------------------------------------------------------------------------------------------------

    CPUSim::CPUSim() {
        // Water rest density = 1 kg/m^3
        //
        // Particle h = 0.05 meter
        // Particle volume = 5.2359 * 10^-4 m^3 = 0.005236 m^3 (4/3 * 0.05^3 * PI)
        // Particle mass = volume * density = 0.005236 kg
        //
        // Smoothing length := 2 * h = 0.1 meter (for now)

        int nrOfP = 64;
        particles1.reserve(nrOfP);
        particles2.reserve(nrOfP);
        float s=0.2; // meter

        int pnrAlongAxis = round(pow(nrOfP, 1.0/3.0));

        nrOfP = pnrAlongAxis * pnrAlongAxis * pnrAlongAxis;

        int r = -1;
        int z = -1;

        float mass = 0.005236f; // kg

        for (int i = 0; i < nrOfP; i++) {
            Particle data;
            if (i % pnrAlongAxis == 0)
                r++;

            if (i % (pnrAlongAxis * pnrAlongAxis) == 0)
                z++;

            data.x = glm::vec3((i % pnrAlongAxis) * s, (r % pnrAlongAxis) * s, z * s);

            data.h = 1;
            data.v = glm::vec3(0, 0, 0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 0;
            data.p = 0;

            particles1.push_back(data);
            particles2.push_back(data);
        }
    }

    void CPUSim::run(int iterNr) {
        for (int i = 0; i < iterNr; i++) {
            calculateRhoAndP();
            iter();
        }
    }

    void CPUSim::setSimProps() {
        simProps.desired_density = 1.0f;
        simProps.k = 0.001f;
        simProps.nr_of_particles = 64;
        simProps.kernelh = 0.1f;

        simProps.gridA = glm::vec4(-1,-1,-1,0);
        simProps.gridB = glm::vec4(1,1,1,0);

        cellx = (int) ceil(abs((simProps.gridB.x - simProps.gridA.x)/simProps.kernelh)); // Number of cells in x direction
        celly = (int) ceil(abs((simProps.gridB.y - simProps.gridA.y)/simProps.kernelh)); // Number of cells in y direction
        cellz = (int) ceil(abs((simProps.gridB.z - simProps.gridA.z)/simProps.kernelh)); // Number of cells in z direction

        std::cout << cellx << " " << celly << " " << cellz << std::endl;
    }

    void CPUSim::printData() {
        auto& particles = first ? particles1 : particles2;

        for (auto& p: particles) {
            std::cout << "Density: " << p.rho <<
                      "; Pressure: " << p.p <<
                      " ; Position: " << p.x.x << " " << p.x.y << " " << p.x.z <<
                      " ; Mass: " << p.m <<
                      " ; Velocity: " << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
        }
    }

    void CPUSim::calculateRhoAndP() {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                density += particles[j].m * kernel(i, j);
            }

            particles[i].rho = density;

            // Pressure calculation
            float ratio = density / simProps.desired_density;
            float p = simProps.k * (pow(ratio, 7) - 1.0);
            particles[i].p = p < 0 ? 0 : p;
        }
    }

    void CPUSim::iter() {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;

        for (int i = 0; i < p1.size(); i++) {
            // Calculate pressure
            glm::vec3 pressure = glm::vec3(0);
            glm::vec3 viscosity = glm::vec3(0);

            for (int j = 0; j < p1.size(); j++) {
                if (j == i)
                    continue;

                // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                float val = 0;
                if (p1[j].rho != 0 && p1[i].rho != 0)
                    val = p1[j].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[j].p / pow(p1[j].rho, 2)));

                pressure += val * gradKernel(i, j);

                glm::vec3 xij = p1[i].x - p1[j].x;

                float pval = 0;
                if (p1[j].rho != 0)
                    pval = (p1[j].m / p1[j].rho) * (dot(xij, gradKernel(i, j)) /  (dot(xij, xij) + 0.01 * simProps.kernelh));

                glm::vec3 vij = p1[i].v - p1[j].v;
                viscosity += pval * vij;
            }

            pressure *= -p1[i].m;

            float nu = 0.01;
            viscosity *= 2 * nu * p1[i].m;

            glm::vec3 gravity(0,0,-9.81 * p1[i].m);
            glm::vec3 acc = (pressure + viscosity + gravity) / p1[i].m;

            float dt = 0.01;
            glm::vec3 viNext = p1[i].v;
            glm::vec3 xiNext = p1[i].x;
            viNext += acc * (float) dt/1.0f;
            xiNext += viNext * dt;

            p2[i].x = xiNext;
            p2[i].v = viNext;
            p2[i].rho = p1[i].rho;
            p2[i].p = p1[i].p;
        }

        first = !first;
    }

    float CPUSim::kernel(int i, int j) {
        auto& particles = first ? particles1 : particles2;

        float q = glm::length(particles[i].x - particles[j].x)/simProps.kernelh;
        /*float oneoverhd = 1.0/(pow(simProps.kernelh, 3));
        float val = 0;

        if (0 <= q && q < 1) {
            val = 2.0/3.0 - pow(q, 2) + pow(q,3)/2;
        } else if (1 <= q && q < 2) {
            val = pow(2-q, 3)/6.0;
        } else {
            val = 0;
        }

        val *= 3.0/(2.0*PI);
        val *= oneoverhd;
        return val;*/
        float m_k = 8.0 / (PI * pow(simProps.kernelh, 3));
        float m_l = 1.0/(pow(simProps.kernelh, 3)) * 3.0/(2.0*PI);

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
    glm::vec3 CPUSim::gradKernel(int i, int j) {
        float m_k = 48.0 / (pow(simProps.kernelh, 3) * PI);
        float m_l = 1.0/(pow(simProps.kernelh, 3)) * 3.0/(2.0*PI);

        auto& particles = first ? particles1 : particles2;

        glm::vec3 r = particles[i].x - particles[j].x;
        float rlength = length(r);
        float q = rlength / simProps.kernelh;
        glm::vec3 ret = glm::vec3(0);

        if (q > 0.0001 && q <= 1.0) {
            glm::vec3 gradq = r / (rlength * simProps.kernelh);

            if (q <= 0.5) {
                gradq *= m_k * q * (3.0 * q - 2);
                ret = gradq;
            } else {
                float factor = 1.0 - q;
                ret = m_k * (-factor * factor) * gradq;
            }
        }

        return ret;

        /*float xijlength = length(particles[i].x - particles[j].x);
        if (xijlength == 0)
            return glm::vec3(0);

        //glm::vec3 xij = normalize(particles[i].x - particles[j].x);
        glm::vec3 dir = (particles[i].x - particles[j].x) / xijlength;

        float val = 0;
        float q = xijlength/simProps.kernelh;
        if (0 <= q && q < 1) {
            val = 3.0/2.0 * pow(q,2) - 2 * q;
        } else if (1 <= q && q < 2) {
            val = -pow(2 - q, 2)/2.0;
        } else {
            val = 0;
        }

        val *= 3.0/(2*PI) * 1.0/(pow(simProps.kernelh, 3));
        dir *= val;

        return dir;*/
    }
}