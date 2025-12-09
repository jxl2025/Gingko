#pragma once

#include <vector>
#include <glm/glm.hpp>

class Cone
{
public:
    void updateParams(int param1, int param2);
    std::vector<float> generateShape() { return m_vertexData; }

private:
    void insertVec3(std::vector<float> &data, glm::vec3 v);
    void setVertexData();

    std::vector<float> m_vertexData;
    int m_param1;
    int m_param2;
    float m_radius = 0.5;

    // add declaration for new helpers
    void makeCapTile(float r0, float r1, float th0, float th1);
    void makeCapSlice(float th0, float th1);
    void makeSlopeTile(const glm::vec3 &TL, const glm::vec3 &TR,
                             const glm::vec3 &BL, const glm::vec3 &BR,
                       float thetaMid);
    void makeSlopeSlice(float th0, float th1);
    void makeWedge(float th0, float th1);
    void makeCone();
};
