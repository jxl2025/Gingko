#version 330 core

in vec3 v_posWorld;
in vec3 v_normWorld;
in vec2 v_uv;//new bump
in vec3 v_tangentW;
in vec3 v_bitangentW;

out vec4 fragColor;

// Global shading coeffs (SceneGlobalData)
uniform float u_globalKa;
uniform float u_globalKd;
uniform float u_globalKs;

uniform sampler2D u_bumpMap; // new bump

// Camera position in world space
uniform vec3 u_camPos;

// Material (SceneMaterial)
uniform vec3 u_kA;       // ambient color (rgb)
uniform vec3 u_kD;       // diffuse color (rgb)
uniform vec3 u_kS;       // specular color (rgb)
uniform float u_shininess;

// NOTE: LightType enum in scenedata.h
// enum class LightType { LIGHT_POINT, LIGHT_DIRECTIONAL, LIGHT_SPOT };
// so we must match these integer values:
const int LIGHT_POINT       = 0;
const int LIGHT_DIRECTIONAL = 1;
const int LIGHT_SPOT        = 2;

const int MAX_LIGHTS = 8;

struct Light {
    int   type;
    vec3  color;
    vec3  function;   // attenuation coefficients: a0, a1, a2
    vec3  pos;        // world-space position (for point/spot)
    vec3  dir;        // world-space direction (for dir/spot); direction of light rays
    float angle;      // spot cone angle (radians)
    float penumbra;   // spot penumbra (radians)
};

uniform int   u_numLights;
uniform Light u_lights[MAX_LIGHTS];

float computeAttenuation(vec3 func, float dist) {
    // 1 / (a0 + a1 * d + a2 * d^2)    (no clamping)
    float denom = func.x + func.y * dist + func.z * dist * dist;
    if (denom <= 0.0) return 1.0; // fallback if badly-formed scene
    return 1.0 / denom;
}

// Spot angular falloff, similar to ray projects.
// Ldir: direction *from point towards light*.
float computeSpotFactor(Light L, vec3 Ldir) {
    // L.dir is direction of emitted rays (from light to scene).
    // Angle between -L.dir and Ldir gives spotlight angle.
    vec3 spotDir = normalize(-L.dir);
    float cosTheta = dot(spotDir, normalize(Ldir));

    float outer = L.angle;                // full cutoff
    float inner = max(0.0, L.angle - L.penumbra); // fully lit inside

    float cosOuter = cos(outer);
    float cosInner = cos(inner);

    if (cosTheta <= cosOuter) return 0.0;
    if (L.penumbra <= 0.0 || cosTheta >= cosInner) return 1.0;

    float t = (cosTheta - cosOuter) / (cosInner - cosOuter);
    return clamp(t, 0.0, 1.0);
}

void main() {
    vec3 N = normalize(v_normWorld);
    vec3 V = normalize(u_camPos - v_posWorld);

    // BUMP MAPPING – height map derivative method
    float height = texture(u_bumpMap, v_uv).r;
    float hU = texture(u_bumpMap, v_uv + vec2(0.001, 0.0)).r;
    float hV = texture(u_bumpMap, v_uv + vec2(0.0, 0.001)).r;

    float scale = 10.0;

    float dU = (hU - height) * scale;
    float dV = (hV - height) * scale;
    vec3 T = normalize(v_tangentW);
    vec3 B = normalize(v_bitangentW);
    vec3 bumpedNormal = normalize(N + dU * T + dV * B); // final shading normal

    N = bumpedNormal;

    // end bump

    // Global ambient, independent of lights
    vec3 color = u_globalKa * u_kA;

    int nLights = min(u_numLights, MAX_LIGHTS);

    for (int i = 0; i < nLights; ++i) {
        Light L = u_lights[i];

        vec3 Ldir;    // direction from point to light
        float dist = 1.0;
        float attenuation = 1.0;
        float spotFactor = 1.0;

        if (L.type == LIGHT_DIRECTIONAL) {
            // dir is direction of rays; from point-to-light is opposite
            Ldir = normalize(-L.dir);
        } else {
            vec3 toLight = L.pos - v_posWorld;
            dist = length(toLight);
            if (dist > 0.0) Ldir = toLight / dist;
            else           Ldir = vec3(0.0, 0.0, 0.0);

            attenuation = computeAttenuation(L.function, dist);

            if (L.type == LIGHT_SPOT) {
                spotFactor = computeSpotFactor(L, Ldir);
            }
        }

        float NdotL = max(dot(N, Ldir), 0.0);
        if (NdotL <= 0.0) continue;

        // Diffuse term
        vec3 diffuse = u_globalKd * u_kD * NdotL;

        // Specular term (Phong)
        vec3 R = reflect(-Ldir, N);
        float RdotV = max(dot(R, V), 0.0);
        float specAmt = (u_shininess > 0.0) ? pow(RdotV, u_shininess) : 0.0;
        vec3 specular = u_globalKs * u_kS * specAmt;

        vec3 lightColor = L.color;
        float weight = attenuation * spotFactor;

        color += (diffuse + specular) * lightColor * weight;
    }

    fragColor = vec4(color, 1.0);
    // color = 0.5 * (N + vec3(1.0)); // remap [-1,1] -> [0,1]
    // fragColor = vec4(color, 1.0);

    // float H = texture(u_bumpMap, v_uv).r;
    // fragColor = vec4(H, H, H, 1.0);
    // vec3 outN = 0.5*(N + vec3(1.0));
    // fragColor = vec4(outN,1.0);


}
