#version 450

#define MAX_STEPS 128
#define SURF_DIST 0.1
#define MAX_DIST 128.0

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

layout(location = 0) out vec4 outColor;
layout(location = 0) in vec3 fragColor;
layout(location = 1) in vec3 fragCoord;
layout(location = 2) in vec2 localPos;

void main() {
    /*vec2 uv = localPos;
    uv.y *= -1;
    vec3 i = (1-smoothstep(0.18, 0.2, length(uv))) * fragColor;
    outColor = vec4(i, 1);*/


    //outColor = vec4(fragColor, 1);

//    vec2 x = localPos + vec2(0.4, 0.4);
    //vec2 avg = vec2(0,0);
    //for (int i = 0; i < 4; i++)
    //    avg += mvp.localPositions[i];
    //vec2 x = avg / 4.0;

//    x.y *= -1;
    //vec3 i = (1 - smoothstep(0.0, 1.0, length(localPos))) * fragColor;
    //vec3 i = 0.5 * fragColor;
    //outColor = vec4(i, 1);

    //vec2 avg = (mvp.localPosition1 + mvp.localPosition2 + mvp.localPosition3 + mvp.localPosition4) / 4;
    vec3 i = (1 - smoothstep(0.18, 0.2, length(localPos))) * fragColor;
    outColor = vec4(i, 1);
}

/*struct Sphere{
    vec3 pos;
    float r;
    vec3 col;
};

struct Line{
    vec3 p1;
    vec3 p2;
    vec3 col;
};

struct Hit{
    float dist;
    vec3 col;
    vec3 pos;
    vec3 normal;
};*/

/*float sphMap(vec3 p, Sphere s) {
    float d = distance(p, vec3(-1, 1, 5)) - 1.0;
    d = min(d, distance(p, vec3(2, 1, 3)) - 0.2);
    d = min(d, distance(p, vec3(0, 1.5, 6)) - 0.32);
    d = min(d, distance(p, vec3(-3, 0.2, 2)) - 0.2);
    d = min(d, p.y + 1.0);
    return d;
}*/

/*Hit sphMap(vec3 p, Sphere s) {
    Hit ret = Hit(
        distance(p, s.pos) - s.r,
        s.col,
        p,
        vec3(0)//getNormal(p)
    );

    return ret;
}*/

/*float lineMap(vec3 c) {
    vec3 p1 = vec3(-1, 1, 5);
    vec3 p2 = vec3(2, 1, 3);
    vec3 p = p1 - p2;

    float t = clamp((dot(p,c - p2)/dot(p,p)), 0.0, 1.0);

    float d = distance(c, t * p + p2) - 0.1;
    return d;
}*/

/*Hit lineMap(vec3 p, Line l) {
    vec3 lineDir = l.p1 - l.p2;
    float t = clamp((dot(lineDir, p - l.p2)/dot(lineDir, lineDir)), 0.0, 1.0);

    Hit ret = Hit(
        distance(p, t * lineDir + l.p2),
        l.col,
        p,
        vec3(0)//getNormal(p)
    );

    return ret;
}*/

/*float map(vec3 p) {
    float d = lineMap(p);
    d = min(d, sphMap(p));
    return d;
}*/

/*Hit map(vec3 p, Line l[], int ls, Sphere s[], int ss) {
    Hit l = lineMap(p);
    return l;
}*/

/*vec3 getNormal(vec3 p) {
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

float getLight(vec3 p) {
    vec3 light = vec3(0, 5, 2);

    vec3 normal = getNormal(p);
    vec3 lightDir = normalize(light - p);
    float dif = dot(normal, lightDir);

    float d = rayMarch(p + normal * SURF_DIST * 2.0, lightDir);
    if (d < length(light - p)) dif *= 0.1;
    return dif;
}*/

/*
void main()
{
    Sphere s1 = Sphere(
        vec3(-1, 1, 5),
        1.0,
        vec3(0, 0, 1.0)
    );

    Sphere s2 = Sphere(
        vec3(2, 1, 3),
        0.2,
        vec3(0, 1.0, 1.0)
    );

    Sphere s3 = Sphere(
        vec3(0, 1.5, 6),
        0.32,
        vec3(0.25, 0.25, 0.80)
    );

    Sphere s4 = Sphere(
        vec3(-3, 0.2, 2),
        0.2,
        vec3(0.9, 0.9, 0.1)
    );

    outColor = vec4(0, 0, 0, 1);
    //vec2 uv = (fragCoord.xy - 0.5 * mvp.res)/mvp.res.y;
    //vec2 uv = (fragCoord.xy - vec2(0.5, 0.5));
    vec2 uv = fragCoord.xy;
    uv.x *= mvp.aspect;
    uv.y *= -1;

    /*vec3 ro = vec3(0,1,-3);
    vec3 rd = normalize(vec3(uv.x, uv.y, 1));

    float d = rayMarch(ro, rd);

    vec3 p = ro + d * rd;
    float dif = getLight(p);
    vec3 col = vec3(dif);
    //col = pow(col, vec3(.4545));
    //vec3 normal = getNormal(p);

    // Output to screen
    outColor = vec4(col, 1.0);
}
*/