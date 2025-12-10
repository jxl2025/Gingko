#pragma once

#include <vector>
#include <glm/glm.hpp>

class Cylinder
{
public:
    void updateParams(int param1, int param2);
    std::vector<float> generateShape() { return m_vertexData; }

private:
    void insertVec3(std::vector<float> &dst, glm::vec3 &pos, const glm::vec3 &n,
                              const glm::vec2 &uv,
                              const glm::vec3 &t,
                    const glm::vec3 &b);
    void pushVertex(std::vector<float> &dst, const glm::vec3 &pos, const glm::vec3 &n,
                    const glm::vec2 &uv,
                    const glm::vec3 &t,
                    const glm::vec3 &b);
    void setVertexData();

    // helpers
    void makeCapSlice(float y, float th0, float th1);
    void makeSideSlice(float th0, float th1);

    std::vector<float> m_vertexData;
    int m_param1;
    int m_param2;
    float m_radius = 0.5;
};
