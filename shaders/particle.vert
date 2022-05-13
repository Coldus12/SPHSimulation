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

    float staticP;
    float padding;
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
layout(location = 3) out float staticP;

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

    vec3 eye = vec3(0, 10, 10);
    vec3 forward = normalize(particlePos - eye);
    vec3 right = normalize(cross(vec3(0,0,1),forward));
    vec3 up = normalize(cross(forward, right));
    vec3 pos;

    switch(int(mod(gl_VertexIndex, 4))) {
        case 0: pos =  - up - right; break;
        case 1: pos =  + up - right; break;
        case 2: pos =  + up + right; break;
        case 3: pos =  - up + right; break;
        default: break;
    }

    pos *= storage_in.p[currentParticleNr].h / 2.0;
    //pos *= 0.1;

    staticP = storage_in.p[currentParticleNr].staticP;
    //gl_Position = mvp.projection * mvp.view * mvp.model * translateTo(particlePos) * vec4(inPosition, 0.0, 1.0);
    gl_Position = mvp.projection * mvp.view * mvp.model * vec4(particlePos.xyz, 1.0) + mvp.projection * mvp.view * vec4(pos.xyz, 1.0);
}