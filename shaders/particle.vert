#version 450
#define PI 3.1415926538

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
layout(location = 0) in vec2 inPosition;

// Output data
//--------------------------
layout(location = 0) out vec2 localPos;
layout(location = 1) out float r;
layout(location = 2) out float c;

mat4 translateTo(vec3 p) {
    return mat4(
        vec4(1, 0, 0, 0),
        vec4(0, 1, 0, 0),
        vec4(0, 0, 1, 0),
        vec4(p.x, p.y, p.z, 1)
    );
}

void main() {
    switch(int(mod(gl_VertexIndex, 4))) {
        case 0: localPos = mvp.localPosition1; break;
        case 1: localPos = mvp.localPosition2; break;
        case 2: localPos = mvp.localPosition3; break;
        case 3: localPos = mvp.localPosition4; break;
        default: break;
    }

    c = gl_VertexIndex / (SimProps.nr_of_particles * 4);

    int currentParticleNr = int(floor(gl_VertexIndex / 4));
    Particle p = storage_in.p[currentParticleNr];

    vec3 particlePos = p.x;
    r = 0.2;

    gl_Position = mvp.projection * mvp.view * mvp.model * translateTo(particlePos) * vec4(inPosition, 0.0, 1.0);
}