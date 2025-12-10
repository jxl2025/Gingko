#include "ShapeMeshes.h"

#include <algorithm>
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "BezierRing.h"
// #include "Tet.h"
// #include "Triangle.h"

static std::vector<float> generatePrimitiveData(const ScenePrimitive &prim,
                                                int p1, int p2)
{
    int param1 = std::max(1, p1);
    int param2 = std::max(3, p2);

    switch (prim.type) {
    case PrimitiveType::PRIMITIVE_CUBE: {
        Cube s;
        s.updateParams(param1);
        return s.generateShape();
    }
    case PrimitiveType::PRIMITIVE_SPHERE: {
        Sphere s;
        s.updateParams(param1, param2);
        return s.generateShape();
    }
    case PrimitiveType::PRIMITIVE_CONE: {
        Cone s;
        s.updateParams(param1, param2);
        return s.generateShape();
    }
    case PrimitiveType::PRIMITIVE_CYLINDER: {
        Cylinder s;
        s.updateParams(param1, param2);
        return s.generateShape();
    }
    // case PrimitiveType::PRIMITIVE_BEZIER_RING: {
    //     BezierRing s;

    //     // Use world-space walk frames if available
    //     if (prim.tubePathId >= 0 &&
    //         prim.tubePathId < (int)g_bezierTubeFrames.size()) {

    //         s.setFrames(g_bezierTubeFrames[prim.tubePathId]);

    //     } else {
    //         // Fallback to straight procedural tube
    //         s.updateParams(param1, param2);
    //     }

    //     return s.generateShape();
    // }
    case PrimitiveType::PRIMITIVE_BEZIER_RING: {
        BezierRing s;

        // 1) Always apply tessellation params from UI
        s.updateParams(param1, param2);

        // 2) If this primitive is one of the venation tubes, attach its frames
        if (prim.tubePathId >= 0 &&
            prim.tubePathId < (int)g_bezierTubeFrames.size()) {

            s.setFrames(g_bezierTubeFrames[prim.tubePathId]);
            // setFrames() calls setVertexData(), which now uses the updated m_param2
        }

        return s.generateShape();
    }
    default:
        // fallback: tiny cube so something shows up instead of nothing
        Cube s;
        s.updateParams(param1);
        return s.generateShape();
    }
}

ShapeMesh createShapeMesh(const ScenePrimitive &prim,
                          int shapeParam1,
                          int shapeParam2)
{
    ShapeMesh mesh;

    std::vector<float> data = generatePrimitiveData(prim, shapeParam1, shapeParam2);
    if (data.empty()) return mesh;

    mesh.vertexCount = static_cast<int>(data.size() / 6); // pos(3) + normal(3)

    glGenVertexArrays(1, &mesh.vao);
    glGenBuffers(1, &mesh.vbo);

    glBindVertexArray(mesh.vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 data.size() * sizeof(float),
                 data.data(),
                 GL_STATIC_DRAW);

    // layout: location 0 -> position, location 1 -> normal
    GLsizei stride = 6 * sizeof(float);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                          stride,
                          reinterpret_cast<void*>(0));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE,
                          stride,
                          reinterpret_cast<void*>(3 * sizeof(float)));

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return mesh;
}

void deleteShapeMesh(ShapeMesh &mesh)
{
    if (mesh.vbo) {
        glDeleteBuffers(1, &mesh.vbo);
        mesh.vbo = 0;
    }
    if (mesh.vao) {
        glDeleteVertexArrays(1, &mesh.vao);
        mesh.vao = 0;
    }
    mesh.vertexCount = 0;
}
