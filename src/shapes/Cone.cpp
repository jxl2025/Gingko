#include "Cone.h"
#include "glm/gtc/constants.hpp"

void Cone::updateParams(int param1, int param2) {
    m_vertexData = std::vector<float>();
    m_param1 = param1;
    m_param2 = param2;
    setVertexData();
}

// Task 8 (cap)
static inline glm::vec3 capPos(float r, float theta) {
    // base is at y = -0.5, radius in [0..0.5]
    return { r * std::cos(theta), -0.5f, r * std::sin(theta) };
}

void Cone::makeCapTile(float r0, float r1, float th0, float th1) {
    // Constant outward normal for the base (faces downward): (0, -1, 0)
    const glm::vec3 n(0.0f, -1.0f, 0.0f);

    glm::vec3 TL = capPos(r0, th0); // "top" means smaller v index; visually ist the inner ring
    glm::vec3 TR = capPos(r0, th1);
    glm::vec3 BL = capPos(r1, th0);
    glm::vec3 BR = capPos(r1, th1);

    // Quad triangulation; one of these becomes degenerate when r0 == 0 (ok per spec)
    auto push = [&](const glm::vec3 &p){
        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, n);
    };

    // Keep CCW as seen from below (outside of the base). With y = -0.5, outside
    // is toward -y; the winding below achieves outward facing for default GL_CCW.
    push(TL); push(BL); push(TR);
    push(TR); push(BL); push(BR);
}

void Cone::makeCapSlice(float th0, float th1) {
    const int rings = std::max(1, m_param1);
    const float rStep = 0.5f / float(rings);
    for (int i = 0; i < rings; ++i) {
        float r0 = i * rStep;
        float r1 = (i + 1) * rStep;
        makeCapTile(r0, r1, th0, th1);
    }
}

// Task 9 (slope
// Implicit-cone slope normal from the handout
static inline glm::vec3 coneImplicitNormal(const glm::vec3 &p) {
    float xNorm =  2.0f * p.x;
    float yNorm = -(1.0f / 4.0f) * (2.0f * p.y - 1.0f);
    float zNorm =  2.0f * p.z;
    return glm::normalize(glm::vec3{xNorm, yNorm, zNorm});
}

// Tip normal: perpendicular to the implicit cone, but horizontal component
// aligned with the wedge face direction (mid-angle), per handout guidance.
static inline glm::vec3 coneTipNormal(float thetaMid) {
    // Build a direction whose xz part points along the face direction at thetaMid,
    // and y points "up" consistent with the implicit normal's sign pattern.
    // We can approximate by normalizing the implicit normal evaluated slightly off the tip.
    const float eps = 1e-4f;
    glm::vec3 nearTip = {
        (0.5f * eps) * std::cos(thetaMid),
        0.5f - eps,  // just below tip
        (0.5f * eps) * std::sin(thetaMid)
    };
    return coneImplicitNormal(nearTip);
}

static inline glm::vec3 slopePos(float y, float theta) {
    // Linear radius falloff: r(y) from 0.5 at y=-0.5 to 0 at y=+0.5
    float t = (y + 0.5f);            // t in [0,1] bottom→top
    float r = 0.5f * (1.0f - t);     // r = 0.5 at bottom, 0 at top
    return { r * std::cos(theta), y, r * std::sin(theta) };
}

// normal will look wierd when the resolution too low
void Cone::makeSlopeTile(const glm::vec3 &TL, const glm::vec3 &TR,
                         const glm::vec3 &BL, const glm::vec3 &BR,
                         float thetaMid) {
    auto vn = [&](const glm::vec3 &p) {
        // Tip special-case: y == +0.5 exactly
        if (std::abs(p.y - 0.5f) < 1e-7f) {
            return coneTipNormal(thetaMid);
        }
        return coneImplicitNormal(p);
    };

    auto push = [&](const glm::vec3 &p){
        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, vn(p));
    };

    // Two tris; keep outward CCW (outside of cone is away from axis)
    // Tri 1: TL, BL, TR   |  Tri 2: TR, BL, BR
    // If your GL front-face is CCW (default), this is outward for our parameterization.
    push(TL); push(BL); push(TR);
    push(TR); push(BL); push(BR);
}

void Cone::makeSlopeSlice(float th0, float th1) {
    const int bands = std::max(3, m_param1);      // vertical tiles
    const float dy  = 1.0f / float(bands);        // y from -0.5 .. +0.5
    const float thetaMid = 0.5f * (th0 + th1);

    for (int i = 0; i < bands; ++i) {
        float yTop = -0.5f + i * dy;
        float yBot = -0.5f + (i + 1) * dy;

        glm::vec3 TL = slopePos(yTop, th0);
        glm::vec3 TR = slopePos(yTop, th1);
        glm::vec3 BL = slopePos(yBot, th0);
        glm::vec3 BR = slopePos(yBot, th1);

        makeSlopeTile(TL, TR, BL, BR, thetaMid);
    }
}

//Task 10
void Cone::makeWedge(float th0, float th1) {
    // Build both parts for a single wedge
    makeCapSlice(th0, th1);
    makeSlopeSlice(th0, th1);
}

void Cone::makeCone() {
    const int wedges = std::max(3, m_param2); // at least 3, per spec
    const float dth  = glm::two_pi<float>() / float(wedges);
    for (int k = 0; k < wedges; ++k) {
        float th0 = k * dth;
        float th1 = (k + 1) * dth;
        makeWedge(th0, th1);
    }
}

void Cone::setVertexData() {
    // TODO for Project 5: Lights, Camera
    makeCone();

}


// Inserts a glm::vec3 into a vector of floats.
// This will come in handy if you want to take advantage of vectors to build your shape!
void Cone::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}
