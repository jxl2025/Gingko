#include "Cylinder.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/constants.hpp>

static inline glm::vec3 cylPos(float r, float y, float theta) {
    return { r * std::cos(theta), y, r * std::sin(theta) };
}

void Cylinder::insertVec3(std::vector<float> &dst, glm::vec3 &pos, const glm::vec3 &n,
                          const glm::vec2 &uv,
                          const glm::vec3 &t,
                          const glm::vec3 &b) {
    // data.push_back(v.x);
    // data.push_back(v.y);
    // data.push_back(v.z);
    dst.push_back(pos.x); dst.push_back(pos.y); dst.push_back(pos.z);
    dst.push_back(n.x);   dst.push_back(n.y);   dst.push_back(n.z);
    dst.push_back(uv.x);  dst.push_back(uv.y);
    dst.push_back(t.x);   dst.push_back(t.y);   dst.push_back(t.z);
    dst.push_back(b.x);   dst.push_back(b.y);   dst.push_back(b.z);

}
void Cylinder::pushVertex(std::vector<float> &dst, const glm::vec3 &pos, const glm::vec3 &n,
                          const glm::vec2 &uv,
                          const glm::vec3 &t,
                          const glm::vec3 &b) {
    dst.push_back(pos.x); dst.push_back(pos.y); dst.push_back(pos.z);
    dst.push_back(n.x);   dst.push_back(n.y);   dst.push_back(n.z);
    dst.push_back(uv.x);  dst.push_back(uv.y);
    dst.push_back(t.x);   dst.push_back(t.y);   dst.push_back(t.z);
    dst.push_back(b.x);   dst.push_back(b.y);   dst.push_back(b.z);

}

// CAP UV + TANGENT/BITANGENT
static inline glm::vec2 capUV(const glm::vec3 &p) {
    return glm::vec2(p.x + 0.5f, p.z + 0.5f);
}

static inline glm::vec3 capTan() { return {1,0,0}; }
static inline glm::vec3 capBit() { return {0,0,1}; }


// ---------- Caps ----------
//
// y = +0.5 → top cap, outward normal +Y
// y = -0.5 → bottom cap, outward normal -Y
//
// We follow the same structure as Cone's base: rings in radius and
// wedges in theta, with CCW winding as viewed from the *outside*
// along the normal direction.
void Cylinder::makeCapSlice(float y, float th0, float th1) {
    // const int   rings = std::max(1, m_param1);
    // const float rStep = 0.5f / float(rings);

    // const glm::vec3 n = (y > 0.0f)
    //                         ? glm::vec3(0.0f,  1.0f, 0.0f)   // top: outward +Y
    //                         : glm::vec3(0.0f, -1.0f, 0.0f);   // bottom: outward -Y

    // auto push = [&](const glm::vec3 &p) {
    //     insertVec3(m_vertexData, p);
    //     insertVec3(m_vertexData, n);
    // };

    // for (int i = 0; i < rings; ++i) {
    //     float r0 = i * rStep;
    //     float r1 = (i + 1) * rStep;

    //     glm::vec3 TL = cylPos(r0, y, th0);
    //     glm::vec3 TR = cylPos(r0, y, th1);
    //     glm::vec3 BL = cylPos(r1, y, th0);
    //     glm::vec3 BR = cylPos(r1, y, th1);

    //     if (y > 0.0f) {
    //         // Top cap: outside is +Y; seen from above (+Y), want CCW.
    //         // Cone's base logic is CCW as seen from *outside* too,
    //         // but there "outside" is -Y. We flip winding relative to that.
    //         //
    //         // Tri 1: TL, TR, BL
    //         // Tri 2: TR, BR, BL
    //         push(TL); push(TR); push(BL);
    //         push(TR); push(BR); push(BL);
    //     } else {
    //         // Bottom cap: outside is -Y; seen from below (-Y), want CCW.
    //         // This is identical to your Cone base.
    //         //
    //         // Tri 1: TL, BL, TR
    //         // Tri 2: TR, BL, BR
    //         push(TL); push(BL); push(TR);
    //         push(TR); push(BL); push(BR);
    //     }
    // }
    int rings = std::max(1, m_param1);
    float rStep = 0.5f / float(rings);

    glm::vec3 normal = (y > 0.f)
                           ? glm::vec3(0,1,0)
                           : glm::vec3(0,-1,0);

    glm::vec3 tan = capTan();
    glm::vec3 bit = capBit();

    auto push = [&](const glm::vec3 &p) {
        glm::vec2 uv = capUV(p);
        pushVertex(m_vertexData, p, normal, uv, tan, bit);
    };

    for (int i = 0; i < rings; ++i) {
        float r0 = i * rStep;
        float r1 = (i+1) * rStep;

        glm::vec3 TL = cylPos(r0, y, th0);
        glm::vec3 TR = cylPos(r0, y, th1);
        glm::vec3 BL = cylPos(r1, y, th0);
        glm::vec3 BR = cylPos(r1, y, th1);

        if (y > 0.f) {
            // top cap CCW from above (+Y)
            push(TL); push(TR); push(BL);
            push(TR); push(BR); push(BL);
        } else {
            // bottom cap CCW from below (-Y)
            push(TL); push(BL); push(TR);
            push(TR); push(BL); push(BR);
        }
    }

}

// SIDE NORMAL / UV / TANGENT / BITANGENT HELPERS
static inline glm::vec3 sideNormal(const glm::vec3 &p) {
    glm::vec3 n = {p.x, 0, p.z};
    float L = glm::length(n);
    return (L < 1e-6f) ? glm::vec3(0,0,0) : n / L;
}

static inline glm::vec2 sideUV(const glm::vec3 &p, float theta) {
    float u = theta / (2.f * M_PI);
    if (u < 0.f) u += 1.f;
    float v = p.y + 0.5f;
    return {u, v};
}

static inline glm::vec3 sideTangent(float theta) {
    return glm::normalize(glm::vec3(-std::sin(theta), 0, std::cos(theta)));
}

static inline glm::vec3 sideBitangent() {
    return {0,1,0};
}


// ---------- Side ----------
//
// We mimic Sphere's approach: build a quad per band/wedge, then
// correct the winding using geometric normal vs an outward direction.
void Cylinder::makeSideSlice(float th0, float th1) {
    // const int   stacks = std::max(1, m_param1);
    // const float dy     = 1.0f / float(stacks);   // y in [-0.5, 0.5]
    // const float r      = 0.5f;

    // auto sideNormal = [](const glm::vec3 &p) {
    //     // radial normal, ignoring y
    //     glm::vec3 n = { p.x, 0.0f, p.z };
    //     float len = glm::length(n);
    //     if (len < 1e-6f) return glm::vec3(0.0f, 0.0f, 0.0f);
    //     return n / len;
    // };

    // auto pushVN = [&](const glm::vec3 &p) {
    //     insertVec3(m_vertexData, p);
    //     insertVec3(m_vertexData, sideNormal(p));
    // };

    // auto emitTri = [&](glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
    //     // Geometric normal
    //     glm::vec3 g = glm::cross(v1 - v0, v2 - v0);
    //     glm::vec3 c = (v0 + v1 + v2) / 3.0f;    // centroid
    //     glm::vec3 out = sideNormal(c);

    //     // If geometric normal points inward (toward -out), flip.
    //     if (glm::dot(g, out) < 0.0f) std::swap(v1, v2);

    //     pushVN(v0); pushVN(v1); pushVN(v2);
    // };

    // for (int i = 0; i < stacks; ++i) {
    //     float yTop = -0.5f + i * dy;
    //     float yBot = -0.5f + (i + 1) * dy;

    //     glm::vec3 TL = cylPos(r, yTop, th0);
    //     glm::vec3 TR = cylPos(r, yTop, th1);
    //     glm::vec3 BL = cylPos(r, yBot, th0);
    //     glm::vec3 BR = cylPos(r, yBot, th1);

    //     // Two tris from the quad; emitTri fixes winding.
    //     emitTri(TL, BL, TR);
    //     emitTri(TR, BL, BR);
    // }
    int stacks = std::max(1, m_param1);
    float dy = 1.f / float(stacks);
    float r = 0.5f;

    auto pushV = [&](const glm::vec3 &p, float th) {
        glm::vec3 n  = sideNormal(p);
        glm::vec2 uv = sideUV(p, th);
        glm::vec3 t  = sideTangent(th);
        glm::vec3 b  = sideBitangent();
        pushVertex(m_vertexData, p, n, uv, t, b);
    };

    auto emitTri = [&](glm::vec3 v0, glm::vec3 v1, glm::vec3 v2,
                       float t0, float t1, float t2)
    {
        glm::vec3 g = glm::cross(v1 - v0, v2 - v0);
        glm::vec3 c = (v0 + v1 + v2) * (1/3.f);
        glm::vec3 out = sideNormal(c);

        if (glm::dot(g, out) < 0.f) std::swap(v1, v2);

        pushV(v0, t0);
        pushV(v1, t1);
        pushV(v2, t2);
    };

    for (int i = 0; i < stacks; ++i) {
        float yTop = -0.5f + i*dy;
        float yBot = -0.5f + (i+1)*dy;

        glm::vec3 TL = cylPos(r, yTop, th0);
        glm::vec3 TR = cylPos(r, yTop, th1);
        glm::vec3 BL = cylPos(r, yBot, th0);
        glm::vec3 BR = cylPos(r, yBot, th1);

        emitTri(TL, BL, TR, th0, th0, th1);
        emitTri(TR, BL, BR, th1, th0, th1);
    }
}

// ---------- Assemble entire cylinder ----------

void Cylinder::setVertexData() {
    m_vertexData.clear();

    // Enforce *effective* minimums here so even if the slider gives 0 or 1,
    // we still get a non-degenerate shape (triangular prism at worst).
    int slices = std::max(3, m_param2);              // ≥ 3 wedges
    float dth  = glm::two_pi<float>() / float(slices);

    for (int k = 0; k < slices; ++k) {
        float th0 = k * dth;
        float th1 = (k + 1) * dth;

        makeCapSlice(+0.5f, th0, th1);
        makeSideSlice(th0, th1);
        makeCapSlice(-0.5f, th0, th1);
    }
}

void Cylinder::updateParams(int param1, int param2) {
    // Clamp to reasonable minimums so the shape is always volumetric.
    m_param1 = std::max(1, param1);   // vertical bands
    m_param2 = std::max(3, param2);   // angular wedges → at least a triangular prism
    setVertexData();
}
