#version 330 core

layout(location = 0) in vec3 inPos;
layout(location = 1) in vec3 inNor;

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_proj;

out vec3 v_posWorld;
out vec3 v_normWorld;

void main() {
    // World-space position
    vec4 posWorld = u_model * vec4(inPos, 1.0);
    v_posWorld = posWorld.xyz;

    // World-space normal: N' = (M^-1)^T * N
    mat3 normalMatrix = mat3(transpose(inverse(u_model)));
    v_normWorld = normalize(normalMatrix * inNor);

    gl_Position = u_proj * u_view * posWorld;
}
