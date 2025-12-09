#include "Cylinder.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/constants.hpp>

static inline glm::vec3 cylPos(float r, float y, float theta) {
    return { r * std::cos(theta), y, r * std::sin(theta) };
}

void Cylinder::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}

// ---------- Caps ----------
//
// y = +0.5 → top cap, outward normal +Y
// y = -0.5 → bottom cap, outward normal -Y
//
// We follow the same structure as Cone's base: rings in radius and
// wedges in theta, with CCW winding as viewed from the *outside*
// along the normal direction.
void Cylinder::makeCapSlice(float y, float th0, float th1) {
    const int   rings = std::max(1, m_param1);
    const float rStep = 0.5f / float(rings);

    const glm::vec3 n = (y > 0.0f)
                            ? glm::vec3(0.0f,  1.0f, 0.0f)   // top: outward +Y
                            : glm::vec3(0.0f, -1.0f, 0.0f);   // bottom: outward -Y

    auto push = [&](const glm::vec3 &p) {
        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, n);
    };

    for (int i = 0; i < rings; ++i) {
        float r0 = i * rStep;
        float r1 = (i + 1) * rStep;

        glm::vec3 TL = cylPos(r0, y, th0);
        glm::vec3 TR = cylPos(r0, y, th1);
        glm::vec3 BL = cylPos(r1, y, th0);
        glm::vec3 BR = cylPos(r1, y, th1);

        if (y > 0.0f) {
            // Top cap: outside is +Y; seen from above (+Y), want CCW.
            // Cone's base logic is CCW as seen from *outside* too,
            // but there "outside" is -Y. We flip winding relative to that.
            //
            // Tri 1: TL, TR, BL
            // Tri 2: TR, BR, BL
            push(TL); push(TR); push(BL);
            push(TR); push(BR); push(BL);
        } else {
            // Bottom cap: outside is -Y; seen from below (-Y), want CCW.
            // This is identical to your Cone base.
            //
            // Tri 1: TL, BL, TR
            // Tri 2: TR, BL, BR
            push(TL); push(BL); push(TR);
            push(TR); push(BL); push(BR);
        }
    }
}

// ---------- Side ----------
//
// We mimic Sphere's approach: build a quad per band/wedge, then
// correct the winding using geometric normal vs an outward direction.
void Cylinder::makeSideSlice(float th0, float th1) {
    const int   stacks = std::max(1, m_param1);
    const float dy     = 1.0f / float(stacks);   // y in [-0.5, 0.5]
    const float r      = 0.5f;

    auto sideNormal = [](const glm::vec3 &p) {
        // radial normal, ignoring y
        glm::vec3 n = { p.x, 0.0f, p.z };
        float len = glm::length(n);
        if (len < 1e-6f) return glm::vec3(0.0f, 0.0f, 0.0f);
        return n / len;
    };

    auto pushVN = [&](const glm::vec3 &p) {
        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, sideNormal(p));
    };

    auto emitTri = [&](glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
        // Geometric normal
        glm::vec3 g = glm::cross(v1 - v0, v2 - v0);
        glm::vec3 c = (v0 + v1 + v2) / 3.0f;    // centroid
        glm::vec3 out = sideNormal(c);

        // If geometric normal points inward (toward -out), flip.
        if (glm::dot(g, out) < 0.0f) std::swap(v1, v2);

        pushVN(v0); pushVN(v1); pushVN(v2);
    };

    for (int i = 0; i < stacks; ++i) {
        float yTop = -0.5f + i * dy;
        float yBot = -0.5f + (i + 1) * dy;

        glm::vec3 TL = cylPos(r, yTop, th0);
        glm::vec3 TR = cylPos(r, yTop, th1);
        glm::vec3 BL = cylPos(r, yBot, th0);
        glm::vec3 BR = cylPos(r, yBot, th1);

        // Two tris from the quad; emitTri fixes winding.
        emitTri(TL, BL, TR);
        emitTri(TR, BL, BR);
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
