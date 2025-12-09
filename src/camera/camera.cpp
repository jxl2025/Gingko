#include "camera.h"
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

Camera::Camera(const SceneCameraData& sc, float aspect)
    : m_eye(glm::vec3(sc.pos))
    , m_look(glm::vec3(sc.look))
    , m_up(glm::normalize(glm::vec3(sc.up)))
    , m_heightAngle(sc.heightAngle)
    , m_aspectRatio(aspect)
    , m_focalLength(sc.focalLength)
    , m_aperture(sc.aperture)
{
    if (glm::length(m_up) < 1e-6f) m_up = glm::vec3(0.f, 1.f, 0.f);
    if (glm::length(m_look) < 1e-6f) m_look = glm::vec3(0.f, 0.f, -1.f);

    const float maxF = glm::radians(175.f);
    if (!(m_heightAngle > 0.f && m_heightAngle < maxF))
        m_heightAngle = glm::radians(45.f);
    if (!(m_aspectRatio > 0.f))
        m_aspectRatio = 1.f;

    orthonormalize();
}

glm::mat4 Camera::getViewMatrix() const {
    glm::vec3 w = glm::normalize(-m_look);
    glm::vec3 u = glm::normalize(glm::cross(m_up, w));
    glm::vec3 v = glm::cross(w, u);

    glm::mat4 view(1.f);
    view[0] = glm::vec4(u, -glm::dot(u, m_eye));
    view[1] = glm::vec4(v, -glm::dot(v, m_eye));
    view[2] = glm::vec4(w, -glm::dot(w, m_eye));
    view[3] = glm::vec4(0.f, 0.f, 0.f, 1.f);
    return view;
}

glm::mat4 Camera::getProjectionMatrix(float nearPlane, float farPlane) const {
    float n = std::max(nearPlane, 0.001f);
    float f = std::max(farPlane, n + 0.001f);

    float t = std::tan(0.5f * m_heightAngle);
    float a = m_aspectRatio;

    glm::mat4 P(0.f);
    P[0][0] = 1.f / (a * t);
    P[1][1] = 1.f / t;
    P[2][2] = -(f + n) / (f - n);
    P[2][3] = -1.f;
    P[3][2] = -(2.f * f * n) / (f - n);
    return P;
}

float Camera::getAspectRatio() const {
    return m_aspectRatio;
}

float Camera::getHeightAngle() const {
    return m_heightAngle;
}

void Camera::setAspectRatio(float aspect) {
    if (aspect > 0.f) m_aspectRatio = aspect;
}

float Camera::getFocalLength() const {
    if (m_focalLength > 0.f) return m_focalLength;
    return 1.f / std::tan(0.5f * m_heightAngle);
}

float Camera::getAperture() const {
    return m_aperture;
}

void Camera::translate(const glm::vec3 &delta) {
    m_eye += delta;
}

glm::vec3 Camera::rotateVector(const glm::vec3 &v,
                               const glm::vec3 &axis,
                               float angle) {
    glm::vec3 a = glm::normalize(axis);
    float c = std::cos(angle);
    float s = std::sin(angle);
    return v * c +
           glm::cross(a, v) * s +
           a * glm::dot(a, v) * (1.f - c);
}

void Camera::yaw(float radians) {
    glm::vec3 axis(0.f, 1.f, 0.f);
    m_look = rotateVector(m_look, axis, radians);
    m_up   = rotateVector(m_up,   axis, radians);
    orthonormalize();
}

void Camera::pitch(float radians) {
    glm::vec3 right = glm::normalize(glm::cross(m_look, m_up));
    m_look = rotateVector(m_look, right, radians);
    m_up   = rotateVector(m_up,   right, radians);
    orthonormalize();
}

void Camera::orthonormalize() {
    glm::vec3 w = glm::normalize(-m_look);
    glm::vec3 u = glm::normalize(glm::cross(m_up, w));
    glm::vec3 v = glm::cross(w, u);

    m_look = -w;
    m_up   = v;
}
