#include "KDTree3D.h"

#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>

KDTree3D::KDTree3D(const std::vector<glm::vec3> &points) {
    build(points);
}

void KDTree3D::build(const std::vector<glm::vec3> &points) {
    m_points = &points;
    m_nodes.clear();
    if (points.empty()) {
        m_root = -1;
        return;
    }

    std::vector<int> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);

    m_nodes.reserve(points.size());
    m_root = buildRecursive(indices, 0, static_cast<int>(indices.size()), 0);
}

int KDTree3D::buildRecursive(std::vector<int> &indices, int start, int end, int depth) {
    if (start >= end) return -1;

    int axis = depth % 3;
    int mid = (start + end) / 2;

    auto cmp = [&](int a, int b) {
        const glm::vec3 &pa = (*m_points)[a];
        const glm::vec3 &pb = (*m_points)[b];
        if (axis == 0) return pa.x < pb.x;
        if (axis == 1) return pa.y < pb.y;
        return pa.z < pb.z;
    };

    std::nth_element(indices.begin() + start,
                     indices.begin() + mid,
                     indices.begin() + end,
                     cmp);

    int nodeIdx = static_cast<int>(m_nodes.size());
    m_nodes.push_back({});
    Node &node = m_nodes.back();
    node.pointIndex = indices[mid];
    node.axis = axis;

    node.left  = buildRecursive(indices, start, mid, depth + 1);
    node.right = buildRecursive(indices, mid + 1, end, depth + 1);

    return nodeIdx;
}

int KDTree3D::nearestNeighbor(const glm::vec3 &q, float &outDist2) const {
    if (empty()) {
        outDist2 = std::numeric_limits<float>::infinity();
        return -1;
    }

    int bestIndex = -1;
    float bestDist2 = std::numeric_limits<float>::infinity();
    nearestRecursive(m_root, q, bestIndex, bestDist2);
    outDist2 = bestDist2;
    return bestIndex;
}

void KDTree3D::nearestRecursive(int nodeIdx,
                                const glm::vec3 &q,
                                int &bestIndex,
                                float &bestDist2) const {
    if (nodeIdx < 0) return;

    const Node &node = m_nodes[nodeIdx];
    const glm::vec3 &p = (*m_points)[node.pointIndex];

    float d2 = glm::length(q - p);
    if (d2 < bestDist2) {
        bestDist2 = d2;
        bestIndex = node.pointIndex;
    }

    float diff;
    if (node.axis == 0) diff = q.x - p.x;
    else if (node.axis == 1) diff = q.y - p.y;
    else diff = q.z - p.z;

    int first  = (diff <= 0.0f) ? node.left  : node.right;
    int second = (diff <= 0.0f) ? node.right : node.left;

    nearestRecursive(first, q, bestIndex, bestDist2);

    if (diff * diff < bestDist2) {
        nearestRecursive(second, q, bestIndex, bestDist2);
    }
}
