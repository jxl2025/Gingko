// BezierRing.h
#pragma once

#include <vector>
#include <glm/glm.hpp>


struct BezierRingFrame {
    glm::vec3 center;
    glm::vec3 normal;
    glm::vec3 binormal;
    float radius;
};

extern std::vector<std::vector<BezierRingFrame>> g_bezierTubeFrames;

class BezierRing
{
public:
    // param1: stacks along the axis (height direction, Y)
    // param2: slices around the ring (angular direction)
    void updateParams(int param1, int param2);
    std::vector<float> generateShape() { return m_vertexData; }
    BezierRing();

    void setFrames(const std::vector<BezierRingFrame>& frames);

private:
    void insertVec3(std::vector<float> &data, glm::vec3 v);
    void setVertexData();

    // Bézier evaluation (scalar radius scale field)
    float evalScale(float u, float v) const;
    float evalScaleDu(float u, float v) const;
    float evalScaleDv(float u, float v) const;

    void initControlNet();

    std::vector<float> m_vertexData;
    int m_param1;  // stacks (height)
    int m_param2;  // slices (angle)

    std::vector<BezierRingFrame> m_frames;
    bool m_useExternalFrames = false;

public:
    // YOU CAN CHANGE THESE:
    // Number of control points in u (height) and v (angle).
    // degree_u = CTRL_U - 1, degree_v = CTRL_V - 1.
    static constexpr int CTRL_U = 4; // e.g. set to 6 for more axial freedom
    static constexpr int CTRL_V = 6; // e.g. 6 control points around

private:
    // Control net: [i][j] → i in u (0..CTRL_U-1), j in v (0..CTRL_V-1)
    // We store a scalar "radius scale" in .x (y/z unused for now).
    glm::vec3 m_ctrl[CTRL_U][CTRL_V];
};
