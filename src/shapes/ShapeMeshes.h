#pragma once

#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>

#include "utils/scenedata.h"  // RenderShapeData, PrimitiveType, ScenePrimitive, etc.

// Forward-declare your shape classes
class Cube;
class Cylinder;
class Cone;
class Sphere;
// add others (Tet, Triangle) if you want them rendered too

struct ShapeMesh {
    GLuint vao = 0;
    GLuint vbo = 0;
    int    vertexCount = 0;   // number of vertices (not floats)
};

ShapeMesh createShapeMesh(const ScenePrimitive &prim,
                          int shapeParam1,
                          int shapeParam2);

void deleteShapeMesh(ShapeMesh &mesh);
