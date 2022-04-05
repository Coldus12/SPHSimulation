#version 450

layout(set = 0, binding = 0) uniform UniformBufferObject {
    mat4 model;
    mat4 view;
    mat4 projection;
    vec2 localPosition1;
    vec2 localPosition2;
    vec2 localPosition3;
    vec2 localPosition4;
    float aspect;
} mvp;

layout(location = 0) in vec2 inPosition;
layout(location = 1) in vec3 inColor;

layout(location = 0) out vec3 fragColor;
layout(location = 1) out vec3 fragCoord;
layout(location = 2) out vec2 localPos;

void main() {
    gl_Position = mvp.projection * mvp.view * mvp.model * vec4(inPosition, 0.0, 1.0);
    //gl_Position = vec4(inPosition, 0.0, 1.0);
    fragColor = inColor;
    fragCoord = (mvp.projection * mvp.view * mvp.model * vec4(inPosition, 0.0, 1.0)).xyz;
    //localPos = mvp.localPositions[int(mod(gl_VertexIndex, 4))];

    switch(int(mod(gl_VertexIndex, 4))) {
        case 0: localPos = mvp.localPosition1; break;
        case 1: localPos = mvp.localPosition2; break;
        case 2: localPos = mvp.localPosition3; break;
        case 3: localPos = mvp.localPosition4; break;
        default: break;
    }
}