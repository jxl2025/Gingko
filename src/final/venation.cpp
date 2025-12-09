#include "final/venation.h"
#include "final/KDTree3D.h"
#include "utils/sceneparser.h"

#include <glm/gtx/transform.hpp>
#include <glm/gtc/constants.hpp>
#include <algorithm>
#include <cmath>
#include <chrono>

namespace venation {

Venation::Venation(const Params &p)
    : m_params(p)
{
    // Seed with time; swap to fixed seed for reproducibility if you like.
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    m_rng.seed(seed);

    // Start with a small “trunk” near the bottom center of the domain
    glm::vec3 center = 0.5f * (m_params.domainMin + m_params.domainMax);
    glm::vec3 base   = glm::vec3(center.x, m_params.domainMin.y, center.z);

    m_nodes.push_back(Node{ base });
    glm::vec3 tip = base + glm::vec3(0.f, m_params.segmentLength, 0.f);
    m_nodes.push_back(Node{ tip });
    m_edges.push_back(Edge{0, 1});

    initSources();
}

glm::vec3 Venation::randomPointInBox()  {
    std::uniform_real_distribution<float> dx(m_params.domainMin.x, m_params.domainMax.x);
    std::uniform_real_distribution<float> dy(m_params.domainMin.y, m_params.domainMax.y);
    std::uniform_real_distribution<float> dz(m_params.domainMin.z, m_params.domainMax.z);
    return glm::vec3(dx(m_rng), dy(m_rng), dz(m_rng));
}

glm::vec3 Venation::jitterDirection(const glm::vec3 &dir) {
    if (m_params.directionJitter <= 0.f) return glm::normalize(dir);

    std::uniform_real_distribution<float> u(0.f, 1.f);
    float phi   = 2.f * glm::pi<float>() * u(m_rng);
    float theta = m_params.directionJitter * (u(m_rng) - 0.5f) * 2.f * 0.5f;

    glm::vec3 up(0.f, 1.f, 0.f);
    glm::vec3 n = glm::normalize(glm::cross(up, dir));
    if (glm::length(n) < 1e-6f) {
        n = glm::normalize(glm::cross(glm::vec3(1.f, 0.f, 0.f), dir));
        if (glm::length(n) < 1e-6f) {
            return glm::normalize(dir);
        }
    }

    glm::vec3 b = glm::normalize(glm::cross(dir, n));

    float sd = std::sin(theta);
    float cd = std::cos(theta);
    glm::vec3 d = cd * glm::normalize(dir) + sd * (std::cos(phi) * n + std::sin(phi) * b);
    return glm::normalize(d);
}

void Venation::initSources() {
    m_sources.clear();
    m_sources.reserve(m_params.initialSources);
    for (int i = 0; i < m_params.initialSources; ++i) {
        m_sources.push_back(randomPointInBox());
    }
}

void Venation::step() {
    if (m_sources.empty() || m_nodes.empty()) return;

    // Build KD-tree over vein nodes
    std::vector<glm::vec3> nodePositions;
    nodePositions.reserve(m_nodes.size());
    for (const auto &n : m_nodes) nodePositions.push_back(n.pos);

    KDTree3D tree(nodePositions);

    int nNodes = static_cast<int>(m_nodes.size());
    std::vector<glm::vec3> accumDir(nNodes, glm::vec3(0.f));
    std::vector<int>       accumCount(nNodes, 0);

    // Assign each source to its nearest vein node
    std::vector<int> sourceNearestIdx(m_sources.size(), -1);
    std::vector<float> sourceNearestDist2(m_sources.size(), 0.f);

    for (std::size_t i = 0; i < m_sources.size(); ++i) {
        float d2 = 0.f;
        int idx = tree.nearestNeighbor(m_sources[i], d2);
        sourceNearestIdx[i] = idx;
        sourceNearestDist2[i] = d2;

        if (idx >= 0 && idx < nNodes) {
            glm::vec3 v = m_nodes[idx].pos;
            glm::vec3 dir = glm::normalize(m_sources[i] - v);
            accumDir[idx]  += dir;
            accumCount[idx] += 1;
        }
    }

    // Create new vein tips
    std::vector<Node> newNodes;
    std::vector<Edge> newEdges;

    // Rebuild tree for min separation checks as we append new nodes
    auto canPlaceNode = [&](const glm::vec3 &p) -> bool {
        std::vector<glm::vec3> allPos;
        allPos.reserve(m_nodes.size() + newNodes.size());
        for (const auto &n : m_nodes)    allPos.push_back(n.pos);
        for (const auto &n : newNodes)   allPos.push_back(n.pos);
        KDTree3D tmpTree(allPos);
        float d2 = 0.f;
        int idx = tmpTree.nearestNeighbor(p, d2);
        if (idx < 0) return true;
        return d2 >= m_params.minNodeSeparation * m_params.minNodeSeparation;
    };

    for (int i = 0; i < nNodes; ++i) {
        if (accumCount[i] == 0) continue;

        glm::vec3 avgDir = accumDir[i] / static_cast<float>(accumCount[i]);
        if (glm::length(avgDir) < 1e-8f) continue;

        avgDir = glm::normalize(avgDir);
        avgDir = jitterDirection(avgDir);

        glm::vec3 p0 = m_nodes[i].pos;
        glm::vec3 p1 = p0 + m_params.segmentLength * avgDir;

        if (!canPlaceNode(p1)) continue;

        int newIndex = static_cast<int>(m_nodes.size() + newNodes.size());
        newNodes.push_back(Node{ p1 });
        newEdges.push_back(Edge{i, newIndex});
    }

    if (!newNodes.empty()) {
        int offset = static_cast<int>(m_nodes.size());
        m_nodes.insert(m_nodes.end(), newNodes.begin(), newNodes.end());
        m_edges.insert(m_edges.end(), newEdges.begin(), newEdges.end());
    }

    // Kill sources that got too close to any vein node
    std::vector<glm::vec3> kept;
    kept.reserve(m_sources.size());
    float kill2 = m_params.killDistance * m_params.killDistance;

    for (std::size_t i = 0; i < m_sources.size(); ++i) {
        if (sourceNearestDist2[i] >= kill2) {
            kept.push_back(m_sources[i]);
        }
    }
    m_sources.swap(kept);
}

void Venation::run() {
    for (int stepIdx = 0; stepIdx < m_params.maxSteps; ++stepIdx) {
        if (m_sources.empty()) break;
        step();
    }
}

// RenderData builder converts vein graph into oriented cylinders

RenderData buildProceduralTreeRenderData() {
    RenderData data;

    // Globals: simple Phong-ish setup
    data.globalData.ka = 0.2f;
    data.globalData.kd = 0.8f;
    data.globalData.ks = 0.4f;
    data.globalData.kt = 0.0f;

    // Camera: put it back a bit, looking toward origin
    SceneCameraData cam{};
    cam.pos   = glm::vec4(0.f, 2.f, 10.f, 1.f);
    cam.look  = glm::vec4(0.f, 0.f, -1.f, 0.f);
    cam.up    = glm::vec4(0.f, 1.f, 0.f, 0.f);
    cam.heightAngle = glm::radians(45.f);
    cam.aperture    = 0.f;
    cam.focalLength = 1.f;
    data.cameraData = cam;

    // One directional light
    SceneLightData L{};
    L.id      = 0;
    L.type    = LightType::LIGHT_DIRECTIONAL;
    L.color   = glm::vec4(1.f);
    L.function = glm::vec3(1.f, 0.f, 0.f);
    L.dir     = glm::vec4(glm::normalize(glm::vec3(-1.f, -1.f, -1.f)), 0.f);
    L.pos     = glm::vec4(0.f, 0.f, 0.f, 1.f);
    L.penumbra = 0.f;
    L.angle    = 0.f;
    L.width    = 0.f;
    L.height   = 0.f;
    data.lights.clear();
    data.lights.push_back(L);

    // Run venation simulation
    Params p;
    p.initialSources    = 200;
    p.maxSteps          = 100;
    p.segmentLength     = 0.18f;
    p.killDistance      = 0.22f;
    p.minNodeSeparation = 0.09f;
    p.knnPerSource      = 1;
    p.directionJitter   = 0.12f;

    Venation sim(p);
    sim.run();

    const auto &nodes = sim.nodes();
    const auto &edges = sim.edges();

    // Basic woody material
    SceneMaterial mat{};
    mat.cAmbient  = glm::vec4(0.10f, 0.06f, 0.03f, 1.f);
    mat.cDiffuse  = glm::vec4(0.55f, 0.35f, 0.20f, 1.f);
    mat.cSpecular = glm::vec4(0.08f, 0.08f, 0.08f, 1.f);
    mat.shininess = 10.f;

    data.shapes.clear();
    data.shapes.reserve(edges.size());

    float r = Venation::CylinderRadius;

    for (const Edge &e : edges) {
        if (e.from < 0 || e.to < 0 ||
            e.from >= static_cast<int>(nodes.size()) ||
            e.to   >= static_cast<int>(nodes.size())) {
            continue;
        }

        glm::vec3 p0 = nodes[e.from].pos;
        glm::vec3 p1 = nodes[e.to].pos;
        glm::vec3 d  = p1 - p0;
        float len    = glm::length(d);
        if (!(len > 1e-4f)) continue;

        glm::vec3 dir = d / len;

        glm::vec3 yAxis(0.f, 1.f, 0.f);
        glm::mat4 R(1.f);

        float cosTheta = glm::dot(yAxis, dir);
        cosTheta = std::clamp(cosTheta, -1.f, 1.f);
        glm::vec3 axis = glm::cross(yAxis, dir);

        if (glm::length(axis) < 1e-8f) {
            if (cosTheta < 0.f) {
                R = glm::rotate(glm::mat4(1.f), glm::pi<float>(), glm::vec3(1.f, 0.f, 0.f));
            }
        } else {
            axis = glm::normalize(axis);
            float angle = std::acos(cosTheta);
            R = glm::rotate(glm::mat4(1.f), angle, axis);
        }

        glm::vec3 center = 0.5f * (p0 + p1);

        glm::mat4 S = glm::scale(glm::mat4(1.f), glm::vec3(r, len, r));
        // Note: Cylinder primitive is [-0.5, 0.5] in Y, so scale by len/2 and rely
        // on model matrix to position the center correctly.
        glm::mat4 T = glm::translate(glm::mat4(1.f), center);

        ScenePrimitive prim{};
        prim.type     = PrimitiveType::PRIMITIVE_CYLINDER;
        prim.material = mat;

        RenderShapeData shape{};
        shape.primitive = prim;
        shape.ctm       = T * R * S;

        data.shapes.push_back(shape);
    }

    return data;
}

} // namespace venation
