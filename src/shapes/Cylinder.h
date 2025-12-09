#pragma once

#include <vector>
#include <glm/glm.hpp>

class Cylinder
{
public:
    void updateParams(int param1, int param2);
    std::vector<float> generateShape() { return m_vertexData; }

private:
    void insertVec3(std::vector<float> &data, glm::vec3 v);
    void setVertexData();

    // helpers
    void makeCapSlice(float y, float th0, float th1);
    void makeSideSlice(float th0, float th1);

    std::vector<float> m_vertexData;
    int m_param1;
    int m_param2;
    float m_radius = 0.5;
};
