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

/*layout(set = 0, binding = 0) uniform ModelViewProj {
    mat4 model;
    mat4 view;
    mat4 projection;
} mvp;*/

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

/*

// Unused input data
//----------------------------------------------------------------------------------------------------------------------
layout(location = 0) in vec3 unusedX;
layout(location = 1) in float unusedH;
layout(location = 2) in vec3 unusedV;
layout(location = 3) in float unusedM;

layout(location = 4) in float unusedRHO;
layout(location = 5) in float unusedP;

layout(location = 6) in float unusedPad1;
layout(location = 7) in float unusedPad2;
//----------------------------------------------------------------------------------------------------------------------

layout(location = 0) out vec3 fragColor;
layout(location = 1) out float drawRadius;*/


layout(location = 0) in vec2 inPosition;
layout(location = 1) in vec3 inColor;

layout(location = 0) out vec3 fragColor;

void main() {
    /*vec3 particlePos = out_data.p[gl_VertexIndex].x;

    vec3 pointOnSurface = particlePos;
    pointOnSurface.x += 1;
    vec4 transformedSurfP = mvp.projection * mvp.view * mvp.model * vec4(pointOnSurface, 1.0);
    vec4 transformedPos = mvp.projection * mvp.view * mvp.model * vec4(particlePos, 1.0);

    gl_Position = transformedPos;
    fragColor = vec3(0.1, 0.1, 0.9);
    drawRadius = length(transformedPos - transformedSurfP);*/

    gl_Position = vec4(inPosition, 0.0, 1.0);
    fragColor = inColor;
}