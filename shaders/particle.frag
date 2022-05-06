#version 450
#define PI 3.1415926538
#define SMOOTHL 0.02

// Particle
//--------------------------
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

// Uniform Buffer Objects
//--------------------------
layout(set = 0, binding = 0) uniform ModelViewProj {
    mat4 model;
    mat4 view;
    mat4 projection;

    vec2 localPosition1;
    vec2 localPosition2;
    vec2 localPosition3;
    vec2 localPosition4;
} mvp;

layout(set = 0, binding = 1) uniform SimulationProperties {
    float desired_density;
    float k;                    // normalization constant / stiffness constant
    float nr_of_particles;
    float kernelh;
} SimProps;

// Storage buffers
//--------------------------
layout(set = 0, binding = 2, std430) readonly buffer inBuffer {
    Particle p[];
} storage_in;

layout(set = 0, binding = 3, std430) buffer outBuffer {
    Particle p[];
} storage_out;

// Input data
//--------------------------
layout(location = 0) in vec2 localPos;
layout(location = 1) in float r;
layout(location = 2) in float c;

// Output data
//--------------------------
layout(location = 0) out vec4 outColor;

void main() {
    vec3 i = (1 - smoothstep(r - SMOOTHL, r, length(localPos))) * vec3(c,3*c/2,2*c/3);
    outColor = vec4(i, 1);
    if ((1 - smoothstep(r - SMOOTHL, r, length(localPos))) < 0.1)
        discard;
}