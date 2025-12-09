#include "Triangle.h"

void Triangle::updateParams() {
    m_vertexData = std::vector<float>();
    setVertexData();
}

void Triangle::setVertexData() {
    // Task 1: update m_vertexData with the vertices and normals
    //         needed to tesselate a triangle
    // Note: you may find the insertVec3 function useful in adding your points into m_vertexData

    // if winding order wrong, tirangle wont show up and when it does, the shading is super dark
    glm::vec3 A(-0.5f,  0.5f, 0.0f);
    glm::vec3 B( 0.5f, -0.5f, 0.0f);
    glm::vec3 C(-0.5f, -0.5f, 0.0f);
    glm::vec3 n(0,0,1);

    insertVec3(m_vertexData, A); insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, C); insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, B); insertVec3(m_vertexData, n);
}

// Inserts a glm::vec3 into a vector of floats.
// This will come in handy if you want to take advantage of vectors to build your shape!
void Triangle::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}
