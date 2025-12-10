#version 330 core

layout(location = 0) in vec3 inPos;
layout(location = 1) in vec3 inNor;
layout(location = 2) in vec2 uv;
layout(location = 3) in vec3 inTangent;
layout(location = 4) in vec3 inBitangent;


uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_proj;

out vec3 v_posWorld;
out vec3 v_normWorld;
out vec2 v_uv;
out vec3 v_tangentW;
out vec3 v_bitangentW;

void main() {
    // World-space position
    vec4 posWorld = u_model * vec4(inPos, 1.0);
    v_posWorld = posWorld.xyz;

    // World-space normal: N' = (M^-1)^T * N
    mat3 normalMatrix = mat3(transpose(inverse(u_model)));
    v_normWorld = normalize(normalMatrix * inNor);

    vec3 T_world = normalize(normalMatrix * inTangent);
    vec3 B_world = normalize(normalMatrix * inBitangent);

    v_tangentW = T_world;
    v_bitangentW = B_world;

    v_uv = uv;

    gl_Position = u_proj * u_view * posWorld;
}
