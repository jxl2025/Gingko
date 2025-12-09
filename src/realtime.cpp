#include "realtime.h"

#include <QCoreApplication>
#include <QMouseEvent>
#include <QKeyEvent>
#include <iostream>
#include "settings.h"
// added these
#include "utils/shaderloader.h"
#include <algorithm>
#include <cmath>
// added for final project
#include "final/venation.h"

// ================== Rendering the Scene!

Realtime::Realtime(QWidget *parent)
    : QOpenGLWidget(parent)
{
    m_prev_mouse_pos = glm::vec2(size().width()/2, size().height()/2);
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    m_keyMap[Qt::Key_W]       = false;
    m_keyMap[Qt::Key_A]       = false;
    m_keyMap[Qt::Key_S]       = false;
    m_keyMap[Qt::Key_D]       = false;
    m_keyMap[Qt::Key_Control] = false;
    m_keyMap[Qt::Key_Space]   = false;

    // If you must use this function, do not edit anything above this
    // add these to avoid spamming rebulids
    m_lastExtraCredit3  = settings.extraCredit3;
    m_lastSceneFilePath = settings.sceneFilePath;
}

// added
glm::mat4 Realtime::makeViewMatrix() const {
    glm::vec3 f = glm::normalize(m_camLook);      // forward
    glm::vec3 u = glm::normalize(m_camUp);
    glm::vec3 s = glm::normalize(glm::cross(f, u)); // right
    u = glm::cross(s, f);

    glm::mat4 R(1.f);
    R[0][0] =  s.x; R[1][0] =  s.y; R[2][0] =  s.z;
    R[0][1] =  u.x; R[1][1] =  u.y; R[2][1] =  u.z;
    R[0][2] = -f.x; R[1][2] = -f.y; R[2][2] = -f.z;

    glm::mat4 T(1.f);
    T[3][0] = -m_camEye.x;
    T[3][1] = -m_camEye.y;
    T[3][2] = -m_camEye.z;

    return R * T;
}

// added
glm::mat4 Realtime::makeProjectionMatrix(float fovy,
                                         float aspect,
                                         float nearP,
                                         float farP) const {
    float f = 1.f / std::tan(0.5f * fovy);

    glm::mat4 P(0.f);
    P[0][0] = f / aspect;
    P[1][1] = f;
    P[2][2] = (farP + nearP) / (nearP - farP);
    P[2][3] = -1.f;
    P[3][2] = (2.f * farP * nearP) / (nearP - farP);
    return P;
}

// added (base assignment)
// void Realtime::rebuildSceneMeshes() {
//     makeCurrent();

//     for (ShapeMesh &m : m_shapeMeshes) {
//         deleteShapeMesh(m);
//     }
//     m_shapeMeshes.clear();

//     int p1 = settings.shapeParameter1;
//     int p2 = settings.shapeParameter2;

//     for (const RenderShapeData &shapeData : m_renderData.shapes) {
//         const ScenePrimitive &prim = shapeData.primitive;
//         ShapeMesh mesh = createShapeMesh(prim, p1, p2);
//         if (mesh.vertexCount > 0) {
//             m_shapeMeshes.push_back(mesh);
//         }
//     }
//     doneCurrent();
// }

// added (extra credit for adaptive detail (both the 3 point and 7 point features)
void Realtime::rebuildSceneMeshes() {
    makeCurrent();

    for (ShapeMesh &m : m_shapeMeshes) {
        deleteShapeMesh(m);
    }
    m_shapeMeshes.clear();

    // Base tessellation from GUI sliders
    int baseP1 = std::max(1, settings.shapeParameter1);
    int baseP2 = std::max(1, settings.shapeParameter2);

    const std::size_t numShapes = m_renderData.shapes.size();

    // ---------- EC1: fewer tris when there are many objects ----------
    float countFactor = 1.f;
    if (settings.extraCredit1 && numShapes > 0) {
        // Up to ~10 objects: full detail; 40+ objects: down to 30% detail.
        const float maxObjects   = 40.f;
        const float minFactorEC1 = 0.3f;

        float t = std::min(static_cast<float>(numShapes) / maxObjects, 1.f);
        // Linear interpolation between 1.0 and minFactorEC1
        countFactor = (1.f - t) + minFactorEC1 * t;
    }

    // For distance-based EC2 we’ll need near/far; fall back to sane defaults
    float nearP = settings.nearPlane;
    float farP  = settings.farPlane;
    if (!(farP > nearP)) {
        nearP = 0.1f;
        farP  = 50.f;
    }

    // ---------- Build a mesh for each shape with its own tessellation ----------
    for (const RenderShapeData &shapeData : m_renderData.shapes) {
        // EC2: distance-based LOD (per-object)
        float distFactor = 1.f;
        if (settings.extraCredit2) {
            // Approximate "center" of the primitive as the CTM-applied origin
            glm::vec3 objPos = glm::vec3(shapeData.ctm * glm::vec4(0.f, 0.f, 0.f, 1.f));
            glm::vec3 eye    = m_camEye;

            float dist = glm::length(objPos - eye);

            // Clamp distance into [nearP, farP] so weird scenes don't explode
            float dClamped = std::max(nearP, std::min(dist, farP));
            float t = (dClamped - nearP) / (farP - nearP); // 0 = near, 1 = far

            // Near objects keep full detail; far ones drop to ~25% detail
            const float minFactorEC2 = 0.25f;
            distFactor = (1.f - t) + minFactorEC2 * t;
        }

        // Combine EC1 + EC2
        float totalFactor = countFactor * distFactor;
        // Never go *too* low or param sliders will feel dead
        const float minOverallFactor = 0.2f;
        if (totalFactor < minOverallFactor) {
            totalFactor = minOverallFactor;
        }

        // Scale base tessellation by the combined factor
        int p1 = std::max(1, static_cast<int>(std::round(baseP1 * totalFactor)));
        int p2 = std::max(3, static_cast<int>(std::round(baseP2 * totalFactor)));

        const ScenePrimitive &prim = shapeData.primitive;
        ShapeMesh mesh = createShapeMesh(prim, p1, p2);
        if (mesh.vertexCount > 0) {
            m_shapeMeshes.push_back(mesh);
        }
    }

    doneCurrent();
}

// added
void Realtime::uploadLights(GLuint program) const {
    const int MAX_LIGHTS = 8;
    int n = static_cast<int>(m_renderData.lights.size());
    if (n > MAX_LIGHTS) n = MAX_LIGHTS;

    GLint locNum = glGetUniformLocation(program, "u_numLights");
    glUniform1i(locNum, n);

    for (int i = 0; i < n; ++i) {
        const SceneLightData &L = m_renderData.lights[i];

        std::string base = "u_lights[" + std::to_string(i) + "].";

        auto loc1i = [program](const std::string &name) {
            return glGetUniformLocation(program, name.c_str());
        };
        auto loc3f = loc1i; // same type; we just reuse

        glUniform1i( loc1i(base + "type"),
                    static_cast<int>(L.type) );

        glUniform3f( loc3f(base + "color"),
                    L.color.r, L.color.g, L.color.b );

        glUniform3f( loc3f(base + "function"),
                    L.function.x, L.function.y, L.function.z );

        glUniform3f( loc3f(base + "pos"),
                    L.pos.x, L.pos.y, L.pos.z );

        glUniform3f( loc3f(base + "dir"),
                    L.dir.x, L.dir.y, L.dir.z );

        glUniform1f( loc1i(base + "angle"),
                    L.angle );

        glUniform1f( loc1i(base + "penumbra"),
                    L.penumbra );
    }

    // Globals
    const SceneGlobalData &g = m_renderData.globalData;
    glUniform1f(glGetUniformLocation(program, "u_globalKa"), g.ka);
    glUniform1f(glGetUniformLocation(program, "u_globalKd"), g.kd);
    glUniform1f(glGetUniformLocation(program, "u_globalKs"), g.ks);
}

// added
void Realtime::uploadMaterial(GLuint program, const SceneMaterial &mat) const {
    glm::vec3 kA(mat.cAmbient.r,  mat.cAmbient.g,  mat.cAmbient.b);
    glm::vec3 kD(mat.cDiffuse.r,  mat.cDiffuse.g,  mat.cDiffuse.b);
    glm::vec3 kS(mat.cSpecular.r, mat.cSpecular.g, mat.cSpecular.b);

    glUniform3f(glGetUniformLocation(program, "u_kA"), kA.x, kA.y, kA.z);
    glUniform3f(glGetUniformLocation(program, "u_kD"), kD.x, kD.y, kD.z);
    glUniform3f(glGetUniformLocation(program, "u_kS"), kS.x, kS.y, kS.z);
    glUniform1f(glGetUniformLocation(program, "u_shininess"), mat.shininess);
}
// added
glm::vec3 rotateVector(const glm::vec3 &v, const glm::vec3 &axis,float angle) {
    glm::vec3 a = glm::normalize(axis);
    float c = std::cos(angle);
    float s = std::sin(angle);
    return v * c +
           glm::cross(a, v) * s +
           a * glm::dot(a, v) * (1.0f - c);
}


void Realtime::finish() {
    killTimer(m_timer);
    this->makeCurrent();

    // Students: anything requiring OpenGL calls when the program exits should be done here

    this->doneCurrent();
}

void Realtime::initializeGL() {
    m_devicePixelRatio = this->devicePixelRatio();

    m_timer = startTimer(1000/60);
    m_elapsedTimer.start();

    // Initializing GL.
    // GLEW (GL Extension Wrangler) provides access to OpenGL functions.
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cerr << "Error while initializing GL: " << glewGetErrorString(err) << std::endl;
    }
    std::cout << "Initialized GL: Version " << glewGetString(GLEW_VERSION) << std::endl;

    // Allows OpenGL to draw objects appropriately on top of one another
    glEnable(GL_DEPTH_TEST);
    // Tells OpenGL to only draw the front face
    glEnable(GL_CULL_FACE);
    // Tells OpenGL how big the screen is
    glViewport(0, 0, size().width() * m_devicePixelRatio, size().height() * m_devicePixelRatio);

    // Students: anything requiring OpenGL calls when the program starts should be done here
    m_phongProgram = ShaderLoader::createShaderProgram(
        ":/resources/shaders/default.vert",
        ":/resources/shaders/default.frag"
        );
}

void Realtime::paintGL() {
    // Students: anything requiring OpenGL calls every frame should be done here
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (!m_phongProgram || m_shapeMeshes.empty() || m_renderData.shapes.empty()) {
        return;
    }

    glUseProgram(m_phongProgram);

    int w = width()  * m_devicePixelRatio;
    int h = height() * m_devicePixelRatio;
    float aspect = (h > 0) ? float(w) / float(h) : 1.f;

    float fovy = m_renderData.cameraData.heightAngle;
    if (!(fovy > 0.f)) fovy = glm::radians(45.f);

    float nearP = settings.nearPlane;
    float farP  = settings.farPlane;

    glm::mat4 V = makeViewMatrix();
    glm::mat4 P = makeProjectionMatrix(fovy, aspect, nearP, farP);

    GLint locV = glGetUniformLocation(m_phongProgram, "u_view");
    GLint locP = glGetUniformLocation(m_phongProgram, "u_proj");
    glUniformMatrix4fv(locV, 1, GL_FALSE, &V[0][0]);
    glUniformMatrix4fv(locP, 1, GL_FALSE, &P[0][0]);

    // Camera position
    glUniform3f(glGetUniformLocation(m_phongProgram, "u_camPos"),
                m_camEye.x, m_camEye.y, m_camEye.z);

    // Lights + global coefficients
    uploadLights(m_phongProgram);

    // Draw each shape
    std::size_t n = std::min(m_shapeMeshes.size(), m_renderData.shapes.size());
    for (std::size_t i = 0; i < n; ++i) {
        const RenderShapeData &shape = m_renderData.shapes[i];
        const ShapeMesh &mesh = m_shapeMeshes[i];

        // Model matrix from scene graph
        glm::mat4 M = shape.ctm;
        glUniformMatrix4fv(glGetUniformLocation(m_phongProgram, "u_model"),
                           1, GL_FALSE, &M[0][0]);

        // Material
        uploadMaterial(m_phongProgram, shape.primitive.material);

        // Geometry
        glBindVertexArray(mesh.vao);
        glDrawArrays(GL_TRIANGLES, 0, mesh.vertexCount);
    }

    glBindVertexArray(0);
}

void Realtime::resizeGL(int w, int h) {
    // Tells OpenGL how big the screen is
    glViewport(0, 0, size().width() * m_devicePixelRatio, size().height() * m_devicePixelRatio);

    // Students: anything requiring OpenGL calls when the program starts should be done here
}

void Realtime::sceneChanged() {
    RenderData newData;

    if (settings.extraCredit3) {
        // Procedural venation tree mode
        newData = venation::buildProceduralTreeRenderData();
    } else {
        // Original behavior: parse scene file
        if (settings.sceneFilePath.empty()) {
            std::cerr << "No scene file path set in settings.\n";
            return;
        }

        bool ok = SceneParser::parse(settings.sceneFilePath, newData);
        if (!ok) {
            std::cerr << "Failed to parse scene file: " << settings.sceneFilePath << "\n";
            return;
        }
    }

    m_renderData = std::move(newData);

    // Camera from renderData (works for both parsed and procedural scenes)
    const SceneCameraData &cam = m_renderData.cameraData;
    m_camEye = glm::vec3(cam.pos);

    glm::vec3 look = glm::vec3(cam.look);
    if (glm::length(look) < 1e-6f) {
        look = glm::vec3(0.f, 0.f, -1.f);
    }
    m_camLook = glm::normalize(look);

    glm::vec3 up = glm::vec3(cam.up);
    if (glm::length(up) < 1e-6f) {
        up = glm::vec3(0.f, 1.f, 0.f);
    }
    m_camUp = glm::normalize(up);

    rebuildSceneMeshes();
    update(); // asks for a PaintGL() call to occur

    // // old version from project
    // if (settings.sceneFilePath.empty()) {
    //     std::cerr << "No scene file path set in settings.\n";
    //     return;
    // }

    // RenderData newData;
    // bool ok = SceneParser::parse(settings.sceneFilePath, newData);
    // if (!ok) {
    //     std::cerr << "Failed to parse scene file: " << settings.sceneFilePath << "\n";
    //     return;
    // }

    // m_renderData = std::move(newData);

    // // Camera from scene
    // const SceneCameraData &cam = m_renderData.cameraData;
    // m_camEye = glm::vec3(cam.pos);
    // glm::vec3 look = glm::vec3(cam.look);
    // if (glm::length(look) < 1e-6f) {
    //     // If the JSON uses 'focus' not 'look', your parser should already have
    //     // converted it; this is just a fallback.
    //     look = glm::vec3(0.f, 0.f, -1.f);
    // }
    // m_camLook = glm::normalize(look);

    // glm::vec3 up = glm::vec3(cam.up);
    // if (glm::length(up) < 1e-6f) {
    //     up = glm::vec3(0.f, 1.f, 0.f);
    // }
    // m_camUp = glm::normalize(up);

    // rebuildSceneMeshes();
    // // added above
    // update(); // asks for a PaintGL() call to occur
}

void Realtime::settingsChanged() {
    bool modeChanged  = (settings.extraCredit3 != m_lastExtraCredit3);
    bool fileChanged  = (!settings.extraCredit3 && settings.sceneFilePath != m_lastSceneFilePath);

    m_lastExtraCredit3  = settings.extraCredit3;
    m_lastSceneFilePath = settings.sceneFilePath;

    if (modeChanged || fileChanged) {
        // Either switched between file/procedural, or chose a new scene file:
        // need to rebuild RenderData and meshes.
        sceneChanged();
    } else {
        // Geometry source is the same; only parameters that affect tessellation /
        // shading / near/far planes changed. Just rebuild VAOs/VBOs and redraw.
        rebuildSceneMeshes();
        update();
    }
    // sceneChanged();
    // // also from old project
    // rebuildSceneMeshes();
    // // from old project
    // update(); // asks for a PaintGL() call to occur
}

// ================== Camera Movement!

void Realtime::keyPressEvent(QKeyEvent *event) {
    m_keyMap[Qt::Key(event->key())] = true;
}

void Realtime::keyReleaseEvent(QKeyEvent *event) {
    m_keyMap[Qt::Key(event->key())] = false;
}

void Realtime::mousePressEvent(QMouseEvent *event) {
    if (event->buttons().testFlag(Qt::LeftButton)) {
        m_mouseDown = true;
        m_prev_mouse_pos = glm::vec2(event->position().x(), event->position().y());
    }
}

void Realtime::mouseReleaseEvent(QMouseEvent *event) {
    if (!event->buttons().testFlag(Qt::LeftButton)) {
        m_mouseDown = false;
    }
}

void Realtime::mouseMoveEvent(QMouseEvent *event) {
    if (m_mouseDown) {
        int posX = event->position().x();
        int posY = event->position().y();
        int deltaX = posX - m_prev_mouse_pos.x;
        int deltaY = posY - m_prev_mouse_pos.y;
        m_prev_mouse_pos = glm::vec2(posX, posY);

        // Use deltaX and deltaY here to rotate
        const float sensitivity = 0.005f;

        float yaw   = -deltaX * sensitivity; // rotate around world Y
        float pitch = -deltaY * sensitivity; // rotate around camera right

        glm::vec3 forward = glm::normalize(m_camLook);
        glm::vec3 up      = glm::normalize(m_camUp);
        glm::vec3 right   = glm::normalize(glm::cross(forward, up));

        // Yaw: world-space Y axis
        glm::vec3 worldY(0.0f, 1.0f, 0.0f);
        forward = rotateVector(forward, worldY, yaw);
        up      = rotateVector(up,      worldY, yaw);

        // Pitch: camera right axis (perpendicular to look and up)
        forward = rotateVector(forward, right, pitch);
        up      = rotateVector(up,      right, pitch);

        // Re-orthonormalize a bit
        forward = glm::normalize(forward);
        right   = glm::normalize(glm::cross(forward, up));
        up      = glm::normalize(glm::cross(right, forward));

        m_camLook = forward;
        m_camUp   = up;

        update(); // asks for a PaintGL() call to occur
    }
}

void Realtime::timerEvent(QTimerEvent *event) {
    int elapsedms   = m_elapsedTimer.elapsed();
    float deltaTime = elapsedms * 0.001f;
    m_elapsedTimer.restart();

    // Use deltaTime and m_keyMap here to move around

    // Move at 5 units/sec according to key spec.
    const float speed = 5.0f;

    glm::vec3 move(0.0f);
    glm::vec3 forward = glm::normalize(m_camLook);
    glm::vec3 upCam   = glm::normalize(m_camUp);
    glm::vec3 right   = glm::normalize(glm::cross(forward, upCam));
    glm::vec3 worldUp(0.0f, 1.0f, 0.0f);

    if (m_keyMap[Qt::Key_W])       move += forward;
    if (m_keyMap[Qt::Key_S])       move -= forward;
    if (m_keyMap[Qt::Key_D])       move += right;
    if (m_keyMap[Qt::Key_A])       move -= right;
    if (m_keyMap[Qt::Key_Space])   move += worldUp;
    if (m_keyMap[Qt::Key_Control]) move -= worldUp;

    if (glm::length(move) > 1e-6f) {
        move = glm::normalize(move);
        m_camEye += move * speed * deltaTime;
    }

    update(); // asks for a PaintGL() call to occur
}

// DO NOT EDIT
void Realtime::saveViewportImage(std::string filePath) {
    // Make sure we have the right context and everything has been drawn
    makeCurrent();

    int fixedWidth = 1024;
    int fixedHeight = 768;

    // Create Frame Buffer
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    // Create a color attachment texture
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, fixedWidth, fixedHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    // Optional: Create a depth buffer if your rendering uses depth testing
    GLuint rbo;
    glGenRenderbuffers(1, &rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, fixedWidth, fixedHeight);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cerr << "Error: Framebuffer is not complete!" << std::endl;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return;
    }

    // Render to the FBO
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, fixedWidth, fixedHeight);

    // Clear and render your scene here
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    paintGL();

    // Read pixels from framebuffer
    std::vector<unsigned char> pixels(fixedWidth * fixedHeight * 3);
    glReadPixels(0, 0, fixedWidth, fixedHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());

    // Unbind the framebuffer to return to default rendering to the screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // Convert to QImage
    QImage image(pixels.data(), fixedWidth, fixedHeight, QImage::Format_RGB888);
    QImage flippedImage = image.mirrored(); // Flip the image vertically

    // Save to file using Qt
    QString qFilePath = QString::fromStdString(filePath);
    if (!flippedImage.save(qFilePath)) {
        std::cerr << "Failed to save image to " << filePath << std::endl;
    }

    // Clean up
    glDeleteTextures(1, &texture);
    glDeleteRenderbuffers(1, &rbo);
    glDeleteFramebuffers(1, &fbo);
}
