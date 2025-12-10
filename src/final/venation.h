#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <random>

// Forward declaration from sceneparser.h
struct RenderData;

namespace venation {

struct Node {
    glm::vec3 pos;
};

struct Edge {
    int from;
    int to;
};

struct Params {
    glm::vec3 domainMin = glm::vec3(-2.f, 0.f, -2.f);
    glm::vec3 domainMax = glm::vec3( 2.f, 4.f,  2.f);

    int   initialSources      = 400;
    int   maxSteps            = 300;
    float segmentLength       = 0.15f;
    float killDistance        = 0.25f;
    float minNodeSeparation   = 0.08f;
    float directionJitter     = 0.1f;  // radians
    int   knnPerSource        = 1;     // nearest vein node only
};

class Venation {
public:
    static constexpr float CylinderRadius = 0.05f; // fixed radius for all cylinders

    explicit Venation(const Params &p = Params());

    void run();

    const std::vector<Node> &nodes() const { return m_nodes; }
    const std::vector<Edge> &edges() const { return m_edges; }

private:
    Params m_params;
    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    std::vector<glm::vec3> m_sources;

    std::mt19937 m_rng;

    glm::vec3 randomPointInBox();
    glm::vec3 jitterDirection(const glm::vec3 &dir);

    void initSources();
    void step();
};

// Build a complete RenderData object containing a venation “tree”
// represented as cylinders, ready for Realtime to render.
RenderData buildProceduralTreeRenderData();

} // namespace venation
