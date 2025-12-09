#pragma once

#include "utils/scenedata.h"
#include <glm/glm.hpp>

class Camera {
public:
    Camera(const SceneCameraData& sc, float aspect);

    glm::mat4 getViewMatrix() const;
    glm::mat4 getProjectionMatrix(float nearPlane, float farPlane) const;

    float getAspectRatio() const;
    float getHeightAngle() const;

    glm::vec3 getEye() const { return m_eye; }
    glm::vec3 getLook() const { return m_look; }   // direction
    glm::vec3 getUp() const   { return m_up; }

    void setAspectRatio(float aspect);
    void translate(const glm::vec3 &delta);
    void yaw(float radians);    // rotate around world Y
    void pitch(float radians);  // rotate around camera right

    float getFocalLength() const;
    float getAperture() const;

private:
    static glm::vec3 rotateVector(const glm::vec3 &v,
                                  const glm::vec3 &axis,
                                  float angle);
    void orthonormalize();

    glm::vec3 m_eye{0.f};
    glm::vec3 m_look{0.f, 0.f, -1.f}; // direction
    glm::vec3 m_up{0.f, 1.f, 0.f};
    float     m_heightAngle{0.78539816339f};
    float     m_aspectRatio{1.f};

    float     m_focalLength{0.f};
    float     m_aperture{0.f};
};
