#pragma once

#include <vector>
#include <glm/glm.hpp>

class KDTree3D {
public:
    KDTree3D() = default;
    explicit KDTree3D(const std::vector<glm::vec3> &points);

    void build(const std::vector<glm::vec3> &points);

    bool empty() const { return m_root < 0 || m_points == nullptr || m_points->empty(); }

    // Returns index into the original points vector, or -1 if tree is empty.
    int nearestNeighbor(const glm::vec3 &q, float &outDist2) const;

private:
    struct Node {
        int pointIndex = -1;
        int left = -1;
        int right = -1;
        int axis = 0;
    };

    const std::vector<glm::vec3> *m_points = nullptr;
    std::vector<Node> m_nodes;
    int m_root = -1;

    int buildRecursive(std::vector<int> &indices, int start, int end, int depth);
    void nearestRecursive(int nodeIdx,
                          const glm::vec3 &q,
                          int &bestIndex,
                          float &bestDist2) const;
};
