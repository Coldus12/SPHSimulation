#version 450
#define PI 3.1415926538

struct Particle {
    vec3 x;                     // position
    float h;                    // radius
    vec3 v;                     // velocity
    float m;                    // mass

    float rho;                  // density
    float p;                    // pressure

    float padding1;
    float padding2;
};

layout(set = 0, binding = 0) uniform SimulationProperties {
    float desired_density;
    float k;                    // normalization constant / stiffness constant
    float nr_of_particles;
    float kernelh;
    float aspect;
} SimProps;

layout(set = 0, binding = 1, std430) readonly buffer inBuffer {
    Particle p[];
} in_data;

layout(set = 0, binding = 2, std430) buffer outBuffer {
    Particle p[];
} out_data;

layout(location = 0) in vec3 fragColor;

void main() {
    
}