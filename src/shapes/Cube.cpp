#include "Cube.h"

void Cube::updateParams(int param1) {
    m_vertexData = std::vector<float>();
    m_param1 = param1;
    setVertexData();
}

void Cube::makeTile(glm::vec3 topLeft,
                    glm::vec3 topRight,
                    glm::vec3 bottomLeft,
                    glm::vec3 bottomRight) {
    // Task 2: create a tile (i.e. 2 triangles) based on 4 given points.
    // e1 = BL - TL, e2 = TR - TL (order matter a lot)
    glm::vec3 e1 = bottomLeft - topLeft;
    glm::vec3 e2 = topRight   - topLeft;
    glm::vec3 n  = glm::normalize(glm::cross(e1, e2));

    // Tri 1: TL, BL, TR
    insertVec3(m_vertexData, topLeft);    insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, bottomLeft); insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, topRight);   insertVec3(m_vertexData, n);

    // Tri 2: TR, BL, BR
    insertVec3(m_vertexData, topRight);    insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, bottomLeft);  insertVec3(m_vertexData, n);
    insertVec3(m_vertexData, bottomRight); insertVec3(m_vertexData, n);
}

static inline glm::vec3 bilerp(const glm::vec3 &TL, const glm::vec3 &TR,
                               const glm::vec3 &BL, const glm::vec3 &BR,
                               float u, float v) {
    // u: left to right in [0,1], v: top to bottom in [0,1]
    glm::vec3 top = glm::mix(TL, TR, u);
    glm::vec3 bot = glm::mix(BL, BR, u);
    return glm::mix(top, bot, v);
}

void Cube::makeFace(glm::vec3 topLeft,
                    glm::vec3 topRight,
                    glm::vec3 bottomLeft,
                    glm::vec3 bottomRight) {
    // Task 3: create a single side of the cube out of the 4
    //         given points and makeTile()
    // Note: think about how param 1 affects the number of triangles on
    //       the face of the cube
    int n = std::max(1, m_param1);
    float step = 1.0f / n;

    for (int i = 0; i < n; ++i) {        // rows (top to bottom)
        for (int j = 0; j < n; ++j) {    // cols (left to right)
            float u0 = j * step,     u1 = (j + 1) * step;
            float v0 = i * step,     v1 = (i + 1) * step;

            glm::vec3 cTL = bilerp(topLeft, topRight, bottomLeft, bottomRight, u0, v0);
            glm::vec3 cTR = bilerp(topLeft, topRight, bottomLeft, bottomRight, u1, v0);
            glm::vec3 cBL = bilerp(topLeft, topRight, bottomLeft, bottomRight, u0, v1);
            glm::vec3 cBR = bilerp(topLeft, topRight, bottomLeft, bottomRight, u1, v1);

            makeTile(cTL, cTR, cBL, cBR);
        }
    }


}

void Cube::setVertexData() {
    // Uncomment these lines for Task 2, then comment them out for Task 3:

    // makeTile(glm::vec3(-0.5f,  0.5f, 0.5f),
    //          glm::vec3( 0.5f,  0.5f, 0.5f),
    //          glm::vec3(-0.5f, -0.5f, 0.5f),
    //          glm::vec3( 0.5f, -0.5f, 0.5f));

    // Uncomment these lines for Task 3:

    // makeFace(glm::vec3(-0.5f,  0.5f, 0.5f),
    //          glm::vec3( 0.5f,  0.5f, 0.5f),
    //          glm::vec3(-0.5f, -0.5f, 0.5f),
    //          glm::vec3( 0.5f, -0.5f, 0.5f));

    // Task 4: Use the makeFace() function to make all 6 sides of the cube

    // cube need6 faces
    // +Z face (front): z = +0.5
    makeFace(glm::vec3(-0.5f,  0.5f,  0.5f),
             glm::vec3( 0.5f,  0.5f,  0.5f),
             glm::vec3(-0.5f, -0.5f,  0.5f),
             glm::vec3( 0.5f, -0.5f,  0.5f));

    // -Z face (back): z = -0.5.
    makeFace(glm::vec3( 0.5f,  0.5f, -0.5f),  // TL when looking from -z
             glm::vec3(-0.5f,  0.5f, -0.5f),  // TR
             glm::vec3( 0.5f, -0.5f, -0.5f),  // BL
             glm::vec3(-0.5f, -0.5f, -0.5f)); // BR

    // +X face (right): x = +0.5
    makeFace(glm::vec3( 0.5f,  0.5f,  0.5f),
             glm::vec3( 0.5f,  0.5f, -0.5f),
             glm::vec3( 0.5f, -0.5f,  0.5f),
             glm::vec3( 0.5f, -0.5f, -0.5f));

    // -X face (left): x = -0.5
    makeFace(glm::vec3(-0.5f,  0.5f, -0.5f),
             glm::vec3(-0.5f,  0.5f,  0.5f),
             glm::vec3(-0.5f, -0.5f, -0.5f),
             glm::vec3(-0.5f, -0.5f,  0.5f));

    // +Y face (top): y = +0.5
    makeFace(glm::vec3(-0.5f,  0.5f, -0.5f),
             glm::vec3( 0.5f,  0.5f, -0.5f),
             glm::vec3(-0.5f,  0.5f,  0.5f),
             glm::vec3( 0.5f,  0.5f,  0.5f));

    // -Y face (bottom): y = -0.5
    makeFace(glm::vec3(-0.5f, -0.5f,  0.5f),
             glm::vec3( 0.5f, -0.5f,  0.5f),
             glm::vec3(-0.5f, -0.5f, -0.5f),
             glm::vec3( 0.5f, -0.5f, -0.5f));
}

// Inserts a glm::vec3 into a vector of floats.
// This will come in handy if you want to take advantage of vectors to build your shape!
void Cube::insertVec3(std::vector<float> &data, glm::vec3 v) {
    data.push_back(v.x);
    data.push_back(v.y);
    data.push_back(v.z);
}
