#include "Sphere.h"
#include "glm/ext/scalar_constants.hpp"
#include "glm/gtc/constants.hpp"

void Sphere::updateParams(int param1, int param2) {
    m_vertexData = std::vector<float>();
    m_param1 = param1;
    m_param2 = param2;
    setVertexData();
}

void Sphere::makeTile(glm::vec3 topLeft,
                      glm::vec3 topRight,
                      glm::vec3 bottomLeft,
                      glm::vec3 bottomRight) {
    // Task 5: Implement the makeTile() function for a Sphere
    // Note: this function is very similar to the makeTile() function for Cube,
    //       but the normals are calculated in a different way!
    auto pushVN = [&](const glm::vec3 &p) {
        insertVec3(m_vertexData, p);
        insertVec3(m_vertexData, glm::normalize(p));
    };

    // Emit one triangle; if its geometric normal points inward (toward origin),
    // swap v1 and v2 so the winding becomes CCW as seen from *outside*.
    auto emitTri = [&](glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
        glm::vec3 g = glm::cross(v1 - v0, v2 - v0);   // geometric normal
        glm::vec3 c = (v0 + v1 + v2) / 3.0f;          // triangle centroid
        if (glm::dot(g, c) < 0.0f) std::swap(v1, v2); // flip to face outward
        pushVN(v0); pushVN(v1); pushVN(v2);
    };

    // Two triangles from the quad. Order doesn’t matter now — emitTri fixes it.
    emitTri(topLeft,  bottomLeft, topRight);
    emitTri(topRight, bottomLeft, bottomRight);
}

void Sphere::makeWedge(float currentTheta, float nextTheta) {
    // Task 6: create a single wedge of the sphere using the
    //         makeTile() function you implemented in Task 5
    // Note: think about how param 1 comes into play here!
    auto wrap2pi = [](float a) {
        float t = fmod(a, glm::two_pi<float>());
        if (t < 0.0f) t += glm::two_pi<float>();
        return t;
    };
    float c = wrap2pi(currentTheta);
    float n = wrap2pi(nextTheta);
    float d = n - c;
    if (d <= 0.0f) d += glm::two_pi<float>();
    const int wedges = std::max(3, m_param2);
    if (d < 1e-6f || std::abs(d - glm::two_pi<float>()) < 1e-6f) {
        d = glm::two_pi<float>() / float(wedges);
        n = c + d;
    }
    // Sphere grid in latitude (param1 bands)
    const float r = 0.5f;
    const int bands = std::max(2, m_param1);
    const float dphi = glm::pi<float>() / float(bands);

    auto sph = [&](float phi, float theta) -> glm::vec3 {
        // spherical-to-Cartesian, minus for z)
        float x = r * glm::sin(phi) * glm::cos(theta);
        float y = r * glm::cos(phi);
        float z = -r * glm::sin(phi) * glm::sin(theta);
        return {x, y, z};
    };

    for (int i = 0; i < bands; ++i) {
        float phiTop =  i      * dphi;
        float phiBot = (i + 1) * dphi;
        glm::vec3 TL = sph(phiTop, c);
        glm::vec3 TR = sph(phiTop, n);
        glm::vec3 BL = sph(phiBot, c);
        glm::vec3 BR = sph(phiBot, n);
        makeTile(TL, TR, BL, BR);
    }
}

void Sphere::makeSphere() {
    // Task 7: create a full sphere using the makeWedge() function you
    //         implemented in Task 6
    // Note: think about how param 2 comes into play here!

    const int wedges = std::max(3, m_param2);                         // ≥3 per spec
    const float dtheta = glm::two_pi<float>() / static_cast<float>(wedges);

    for (int k = 0; k < wedges; ++k) {
        float currentTheta = k * dtheta;
        float nextTheta    = (k + 1) * dtheta;
        makeWedge(currentTheta, nextTheta);
    }
}

void Sphere::setVertexData() {
    // Uncomment these lines to make a wedge for Task 6, then comment them out for Task 7:

    // float thetaStep = glm::radians(360.f / m_param2);
    // float currentTheta = 0 * thetaStep;
    // float nextTheta = 1 * thetaStep;
    // makeWedge(currentTheta, nextTheta);

    // Uncomment these lines to make sphere for Task 7:

    makeSphere();
}

// Inserts a glm::vec3 into a vector of floats.
// This will come in handy if you want to take advantage of vectors to build your shape!
void Sphere::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}
