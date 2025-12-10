#include "final/venation.h"
#include "final/KDTree3D.h"
#include "shapes/BezierRing.h"
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

// glm::vec3 Venation::randomPointInBox()  {
//     std::uniform_real_distribution<float> dx(m_params.domainMin.x, m_params.domainMax.x);
//     std::uniform_real_distribution<float> dy(m_params.domainMin.y, m_params.domainMax.y);
//     std::uniform_real_distribution<float> dz(m_params.domainMin.z, m_params.domainMax.z);
//     return glm::vec3(dx(m_rng), dy(m_rng), dz(m_rng));
// }

// to make the auxin distribute some way above the initial veins
// glm::vec3 Venation::randomPointInBox() {
//     // Compute global center and extents of the domain.
//     glm::vec3 center = 0.5f * (m_params.domainMin + m_params.domainMax);
//     float widthX  = m_params.domainMax.x - m_params.domainMin.x;
//     float widthZ  = m_params.domainMax.z - m_params.domainMin.z;
//     float heightY = m_params.domainMax.y - m_params.domainMin.y;

//     // ---- You can tweak these two knobs ----
//     // Fraction of the total height where auxin slab is centered.
//     const float slabCenterFrac  = 0.75f;  // 0 = at base, 1 = at top
//     // Fraction of total height used for slab thickness.
//     const float slabThicknessFrac = 0.30f; // 30% of height

//     float slabCenterY = m_params.domainMin.y + slabCenterFrac * heightY;
//     float slabHalfY   = 0.5f * slabThicknessFrac * heightY;

//     // Optionally tighten X/Z so auxin are near the trunk horizontally.
//     float halfX = 0.4f * widthX;
//     float halfZ = 0.4f * widthZ;

//     std::uniform_real_distribution<float> dx(center.x - halfX, center.x + halfX);
//     std::uniform_real_distribution<float> dz(center.z - halfZ, center.z + halfZ);
//     std::uniform_real_distribution<float> dy(slabCenterY - slabHalfY,
//                                              slabCenterY + slabHalfY);

//     return glm::vec3(dx(m_rng), dy(m_rng), dz(m_rng));
// }

// to use sphere instead
glm::vec3 Venation::randomPointInBox() {
    glm::vec3 center = 0.5f * (m_params.domainMin + m_params.domainMax);

    float heightY = m_params.domainMax.y - m_params.domainMin.y;
    // Center the sphere higher up in the box.
    float yCenter = m_params.domainMin.y + 0.75f * heightY;
    glm::vec3 sphereCenter(center.x, yCenter, center.z);

    // Radius: something modest relative to domain size.
    float radius = 0.6f * heightY;

    std::uniform_real_distribution<float> u(0.f, 1.f);
    std::uniform_real_distribution<float> angle01(0.f, 1.f);

    // Sample direction uniformly on sphere.
    float z   = 2.f * u(m_rng) - 1.f;
    float phi = 2.f * glm::pi<float>() * angle01(m_rng);
    float rXY = std::sqrt(std::max(0.f, 1.f - z*z));
    glm::vec3 dir(rXY * std::cos(phi), z, rXY * std::sin(phi));

    // Radius with cubic root so density is uniform inside the sphere.
    float r = radius * std::cbrt(u(m_rng));

    return sphereCenter + r * dir;
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

struct Walk {
    std::vector<int> nodeIndices;
};

// Build all maximal degree-2 walks in the vein tree.
// Endpoints have degree != 2 (root, leaves, branch junctions).
// Interior nodes in a walk have degree == 2.
static std::vector<Walk> buildWalks(const std::vector<Node> &nodes,
                                    const std::vector<Edge> &edges)
{
    const int n = static_cast<int>(nodes.size());
    std::vector<std::vector<int>> adj(n);
    std::vector<int> degree(n, 0);

    // Undirected adjacency from edges
    for (const Edge &e : edges) {
        if (e.from < 0 || e.to < 0 || e.from >= n || e.to >= n) {
            continue;
        }
        adj[e.from].push_back(e.to);
        adj[e.to].push_back(e.from);
        degree[e.from]++;
        degree[e.to]++;
    }

    auto isEndpoint = [&](int v) {
        return degree[v] != 2;
    };

    std::vector<bool> visited(n, false);
    std::vector<Walk> walks;

    // For every endpoint, launch walks along each neighbor.
    for (int v = 0; v < n; ++v) {
        if (!isEndpoint(v)) continue;

        for (int nbr : adj[v]) {
            // If the neighbor is already part of some walk interior, skip
            if (visited[nbr]) continue;

            Walk w;
            int prev = v;
            int cur  = nbr;

            w.nodeIndices.push_back(prev);

            while (true) {
                w.nodeIndices.push_back(cur);
                visited[cur] = true;

                if (isEndpoint(cur)) {
                    // Reached leaf or branching node: end of this walk
                    break;
                }

                // degree(cur) == 2: pick the neighbor that is not prev
                int next = (adj[cur][0] == prev) ? adj[cur][1] : adj[cur][0];
                prev = cur;
                cur  = next;
            }

            if (w.nodeIndices.size() >= 2) {
                walks.push_back(std::move(w));
            }
        }
    }

    return walks;
}

static void computeParallelTransportFrames(
    const std::vector<glm::vec3>& centers,
    std::vector<glm::vec3>& normals,
    std::vector<glm::vec3>& binormals)
{
    int n = centers.size();
    normals.resize(n);
    binormals.resize(n);

    glm::vec3 t0 = glm::normalize(centers[1] - centers[0]);
    glm::vec3 arbitrary = (std::abs(t0.y) < 0.99f) ? glm::vec3(0,1,0) : glm::vec3(1,0,0);
    normals[0]   = glm::normalize(glm::cross(t0, arbitrary));
    binormals[0]= glm::normalize(glm::cross(t0, normals[0]));

    for (int i = 1; i < n; ++i) {
        glm::vec3 ti = glm::normalize(centers[i] - centers[i-1]);
        glm::vec3 v = glm::cross(t0, ti);

        if (glm::length(v) < 1e-5f) {
            normals[i]   = normals[i-1];
            binormals[i]= binormals[i-1];
        } else {
            float a = std::acos(glm::clamp(glm::dot(t0, ti), -1.f, 1.f));
            v = glm::normalize(v);
            glm::mat3 R = glm::mat3(glm::rotate(glm::mat4(1.f), a, v));
            normals[i]   = R * normals[i-1];
            binormals[i]= R * binormals[i-1];
        }
        t0 = ti;
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
    p.maxSteps          = 80;
    p.segmentLength     = 0.18f;
    p.killDistance      = 0.22f;
    p.minNodeSeparation = 0.09f;
    p.knnPerSource      = 1;
    p.directionJitter   = 0.12f;

    Venation sim(p);
    sim.run();

    const auto &nodes = sim.nodes();
    const auto &edges = sim.edges();

    // ------------------------------------------------------------
    // Flood-fill style "flow" per node: how many leaf tips descend
    // through this node. Used later to drive branch thickness.
    // ------------------------------------------------------------
    const int nNodes = static_cast<int>(nodes.size());
    std::vector<int>   parent(nNodes, -1);
    std::vector<int>   outDegree(nNodes, 0);

    for (const Edge &e : edges) {
        if (e.from < 0 || e.to < 0 ||
            e.from >= nNodes || e.to >= nNodes) {
            continue;
        }
        parent[e.to]   = e.from;   // tree structure: one parent per node
        outDegree[e.from]++;       // count children
    }

    std::vector<float> flow(nNodes, 0.f);

    // Leaves (no children) start with flow = 1
    for (int i = 0; i < nNodes; ++i) {
        if (outDegree[i] == 0) {
            flow[i] = 1.f;
        }
    }

    // Accumulate flow upwards from leaves to root(s).
    // Thanks to construction, parents always have smaller indices,
    // so reverse index order gives a valid post-order.
    for (int i = nNodes - 1; i >= 0; --i) {
        int p = parent[i];
        if (p >= 0) {
            flow[p] += flow[i];
        }
    }

    // Find max flow for normalization
    float maxFlow = 0.f;
    for (int i = 0; i < nNodes; ++i) {
        maxFlow = std::max(maxFlow, flow[i]);
    }
    if (maxFlow <= 0.f) {
        maxFlow = 1.f;
    }

    // 1) Per-edge vs per-walk branches:
    //    - false: original behavior (one primitive per Edge)
    //    - true : new behavior   (one primitive per maximal degree-2 Walk)
    constexpr bool kUseWalkBranches      = false;

    // 2) Cylinder vs BezierRing for the actual tube primitive:
    constexpr bool kUseBezierRingBranches = false;

    // Basic woody material
    SceneMaterial mat{};
    mat.cAmbient  = glm::vec4(0.10f, 0.06f, 0.03f, 1.f);
    mat.cDiffuse  = glm::vec4(0.55f, 0.35f, 0.20f, 1.f);
    mat.cSpecular = glm::vec4(0.08f, 0.08f, 0.08f, 1.f);
    mat.shininess = 10.f;

    // Light grey floor material
    SceneMaterial floorMat{};
    floorMat.cAmbient  = glm::vec4(0.3f, 0.3f, 0.3f, 1.f);
    floorMat.cDiffuse  = glm::vec4(0.7f, 0.7f, 0.7f, 1.f);
    floorMat.cSpecular = glm::vec4(0.05f, 0.05f, 0.05f, 1.f);
    floorMat.shininess = 5.f;

    data.shapes.clear();

    //add floor
    {
        // Use the same domain extents as the venation sim
        float floorThickness = 0.1f; // tweak as you like

        float width  = p.domainMax.x - p.domainMin.x;
        float depth  = p.domainMax.z - p.domainMin.z;

        // Centered under the tree; top of the floor at domainMin.y
        glm::vec3 floorCenter;
        floorCenter.x = 0.5f * (p.domainMin.x + p.domainMax.x);
        floorCenter.z = 0.5f * (p.domainMin.z + p.domainMax.z);
        floorCenter.y = p.domainMin.y - 0.5f * floorThickness;

        ScenePrimitive floorPrim{};
        floorPrim.type     = PrimitiveType::PRIMITIVE_CUBE;
        floorPrim.material = floorMat;

        // Cube primitive is unit cube centered at origin; scale to our floor size
        glm::mat4 S_floor = glm::scale(glm::mat4(1.f),
                                       glm::vec3(width, floorThickness, depth));
        glm::mat4 T_floor = glm::translate(glm::mat4(1.f), floorCenter);

        RenderShapeData floorShape{};
        floorShape.primitive = floorPrim;
        floorShape.ctm       = T_floor * S_floor;

        data.shapes.push_back(floorShape);
    }

    // float r = Venation::CylinderRadius;
    float baseR = Venation::CylinderRadius;

    if (!kUseWalkBranches) {
        // ------------------------------------------------------------
        // ORIGINAL: one primitive per Edge
        // ----------------------------------------------------------
        int edgeCount = static_cast<int>(edges.size());
        data.shapes.reserve(edges.size());
        const float radiusScaleMin = 0.25f; // thinnest branches
        const float radiusScaleMax = 2.0f;  // thickest trunk

        // for (const Edge &e : edges) {
        for (int ei = 0; ei < edgeCount; ++ei) {
            const Edge &e = edges[ei];
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

            // float r = Venation::CylinderRadius;
            // float t = (edgeCount > 1) ? (float)ei / (float)(edgeCount - 1) : 0.f;
            // // version 1: linear scaling
            // // Map t in [0,1] to a radius scale.
            // // Example: oldest = 1.5 * baseR, newest = 0.5 * baseR
            // float r = baseR * (1.5f - t);  // linear: r in [0.5, 1.5]*baseR
            // // version 2: nonlinear
            // // Make the falloff non-linear: most edges get thin quickly.
            // // Try exponent between 2.0 and 3.0 for strong contrast.
            // float p = std::pow(t, 2.5f);
            // // Map to a radius scale:
            // //  - oldest edge: ~2.0 * baseR
            // //  - newest edge: ~0.2 * baseR
            // float radiusScale = 0.2f + 1.8f * (1.0f - p);
            // float r = baseR * radiusScale;
            // version 3: flood fill flow
            float f = flow[e.from]; // how many leaves flow through this branch
            // Normalize with log so large flows saturate visually:
            float t = std::log(1.f + f) / std::log(1.f + maxFlow); // in [0,1]
            float radiusScale =
                radiusScaleMin + (radiusScaleMax - radiusScaleMin) * t;
            float r = baseR * radiusScale;
            r = std::max(0.02f * baseR, r); // avoid degenerate tiny radii

            glm::vec3 yAxis(0.f, 1.f, 0.f);
            glm::mat4 R(1.f);

            float cosTheta = glm::dot(yAxis, dir);
            cosTheta = std::clamp(cosTheta, -1.f, 1.f);
            glm::vec3 axis = glm::cross(yAxis, dir);

            if (glm::length(axis) < 1e-8f) {
                if (cosTheta < 0.f) {
                    R = glm::rotate(glm::mat4(1.f),
                                    glm::pi<float>(),
                                    glm::vec3(1.f, 0.f, 0.f));
                }
            } else {
                axis = glm::normalize(axis);
                float angle = std::acos(cosTheta);
                R = glm::rotate(glm::mat4(1.f), angle, axis);
            }

            glm::vec3 center = 0.5f * (p0 + p1);

            glm::mat4 S = glm::scale(glm::mat4(1.f), glm::vec3(r, len, r));
            glm::mat4 T = glm::translate(glm::mat4(1.f), center);

            ScenePrimitive prim{};
            if (kUseBezierRingBranches) {
                prim.type = PrimitiveType::PRIMITIVE_BEZIER_RING;
            } else {
                prim.type = PrimitiveType::PRIMITIVE_CYLINDER;
            }
            prim.material = mat;

            RenderShapeData shape{};
            shape.primitive = prim;
            shape.ctm       = T * R * S;

            data.shapes.push_back(shape);
        }
    } else {
        // ------------------------------------------------------
        // NEW: one primitive per maximal degree-2 Walk
        // ------------------------------------------------------------
        auto walks = buildWalks(nodes, edges);
        data.shapes.reserve(walks.size());

        for (const Walk &w : walks) {
            std::vector<glm::vec3> centers;
            for (int idx : w.nodeIndices) {
                centers.push_back(nodes[idx].pos);
            }

            std::vector<glm::vec3> normals, binormals;
            computeParallelTransportFrames(centers, normals, binormals);

            std::vector<BezierRingFrame> frames;
            float r = 0.5f * Venation::CylinderRadius;

            for (int i = 0; i < (int)centers.size(); ++i) {
                frames.push_back({
                    centers[i],
                    normals[i],
                    binormals[i],
                    r
                });
            }

            //Register this walk globally
            int pathId = (int)g_bezierTubeFrames.size();
            g_bezierTubeFrames.push_back(frames);

            ScenePrimitive prim{};
            prim.type = PrimitiveType::PRIMITIVE_BEZIER_RING;
            prim.material = mat;

            // CONNECT THIS PRIMITIVE TO ITS FRAME PATH
            prim.tubePathId = pathId;

            RenderShapeData shape{};
            shape.primitive = prim;

            // World-space geometry. identity transform is correct
            shape.ctm = glm::mat4(1.0f);

            data.shapes.push_back(shape);
        }
    }

    return data;
}

} // namespace venation
