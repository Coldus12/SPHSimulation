#version 450

#define MAX_STEPS 128
#define SURF_DIST 0.1
#define MAX_DIST 128.0

layout(set = 0, binding = 0) uniform UniformBufferObject {
    mat4 model;
    mat4 view;
    mat4 projection;
    float aspect;
} mvp;

layout(location = 0) out vec4 outColor;
layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragCoord;

float sphMap(vec3 p) {
    float d = distance(p, vec3(-1, 1, 5)) - 1.0;
    d = min(d, distance(p, vec3(2, 1, 3)) - 0.2);
    d = min(d, distance(p, vec3(0, 1.5, 6)) - 0.32);
    d = min(d, distance(p, vec3(-3, 0.2, 2)) - 0.2);
    d = min(d, p.y + 1.0);

    return d;
}

float lineMap(vec3 c) {
    vec3 p1 = vec3(-1, 1, 5);
    vec3 p2 = vec3(2, 1, 3);
    vec3 p = p1 - p2;

    float t = clamp((dot(p,c - p2)/dot(p,p)), 0.0, 1.0);

    float d = distance(c, t * p + p2) - 0.1;
    return d;
}

float map(vec3 p) {
    float d = lineMap(p);
    d = min(d, sphMap(p));
    return d;
}

vec3 getNormal(vec3 p) {
    vec2 e = vec2(0.01, 0);
    float d = map(p);
    vec3 n = d - vec3(
        map(p - e.xyy),
        map(p - e.yxy),
        map(p - e.yyx)
    );
    return normalize(n);
}

float rayMarch(vec3 ro, vec3 rd) {
    float dO=0.0;

    for(int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd*dO;
        float dS = map(p);
        dO += dS;
        if(dO>MAX_DIST || dS<SURF_DIST) break;
    }

    return dO;
}

void main()
{
    outColor = vec4(0, 0, 0, 1);
    //vec2 uv = (fragCoord.xy - 0.5 * mvp.res)/mvp.res.y;
    //vec2 uv = (fragCoord.xy - vec2(0.5, 0.5));
    vec2 uv = fragCoord.xy;
    uv.x *= mvp.aspect;

    vec3 ro = vec3(0,1,-3);
    vec3 rd = normalize(vec3(uv.x, uv.y, 1));

    float d = rayMarch(ro, rd);

    vec3 p = ro + d * rd;
    vec3 normal = getNormal(p);

    // Output to screen
    outColor = vec4(normal,1.0);
}