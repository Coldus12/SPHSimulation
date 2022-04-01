#version 450

layout(set = 0, binding = 0) uniform UniformBufferObject {
    mat4 model;
    mat4 view;
    mat4 projection;
    float aspect;
} mvp;

layout(location = 0) in vec2 inPosition;
layout(location = 1) in vec3 inColor;

layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragCoord;

void main() {
    gl_Position = mvp.projection * mvp.view * mvp.model * vec4(inPosition, 0.0, 1.0);
    //gl_Position = vec4(inPosition, 0.0, 1.0);
    fragColor = inColor;
    fragCoord = (mvp.projection * mvp.view * mvp.model * vec4(inPosition, 0.0, 1.0)).xyz;
}