#pragma once

// Defined before including GLEW to suppress deprecation messages on macOS
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GL/glew.h>
#include <glm/glm.hpp>

#include <unordered_map>
#include <QElapsedTimer>
#include <QOpenGLWidget>
#include <QTime>
#include <QTimer>

// added these
#include "utils/sceneparser.h"   // RenderData, SceneParser
#include "shapes/ShapeMeshes.h"   // ShapeMesh

class Realtime : public QOpenGLWidget
{
public:
    Realtime(QWidget *parent = nullptr);
    void finish();                                      // Called on program exit
    void sceneChanged();
    void settingsChanged();
    void saveViewportImage(std::string filePath);

public slots:
    void tick(QTimerEvent* event);                      // Called once per tick of m_timer

protected:
    void initializeGL() override;                       // Called once at the start of the program
    void paintGL() override;                            // Called whenever the OpenGL context changes or by an update() request
    void resizeGL(int width, int height) override;      // Called when window size changes

private:
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void timerEvent(QTimerEvent *event) override;

    // Tick Related Variables
    int m_timer;                                        // Stores timer which attempts to run ~60 times per second
    QElapsedTimer m_elapsedTimer;                       // Stores timer which keeps track of actual time between frames

    // Input Related Variables
    bool m_mouseDown = false;                           // Stores state of left mouse button
    glm::vec2 m_prev_mouse_pos;                         // Stores mouse position
    std::unordered_map<Qt::Key, bool> m_keyMap;         // Stores whether keys are pressed or not

    // Device Correction Variables
    double m_devicePixelRatio;

    // added
    RenderData m_renderData;
    std::vector<ShapeMesh> m_shapeMeshes;

    // Camera state (to be) derived from m_renderData.cameraData
    glm::vec3 m_camEye  {0.f, 0.f, 5.f};
    glm::vec3 m_camLook {0.f, 0.f, -1.f};  // direction
    glm::vec3 m_camUp   {0.f, 1.f, 0.f};

    // Helper to rebuild all meshes from m_renderData.shapes
    void rebuildSceneMeshes();

    // Helpers for lighting
    GLuint m_phongProgram = 0;
    void uploadLights(GLuint program) const;
    void uploadMaterial(GLuint program, const SceneMaterial &mat) const;

    // Helper to make view/projection matrices without glm::lookAt/perspective
    glm::mat4 makeViewMatrix() const;
    glm::mat4 makeProjectionMatrix(float fovy, float aspect, float nearP, float farP) const;

    // reduce unnecessary rebuilding of vao/vbos!
    bool        m_lastExtraCredit3 = false;
    std::string m_lastSceneFilePath;
};
