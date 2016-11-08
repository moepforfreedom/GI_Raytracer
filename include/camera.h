#pragma once

#include <glm/glm.hpp>

/// Represents the camera with information about the 'sensor' size.
struct Camera {
    explicit Camera(glm::dvec3 pos) : Camera(pos, {0, 0, 0}) {}
    Camera(glm::dvec3 pos, glm::dvec3 lookAt) : pos(pos), up({0, 0, 1.0}), forward(lookAt - pos) {
        forward = glm::normalize(forward);
    }

    glm::dvec3 pos;
    glm::dvec3 up;
    glm::dvec3 forward;              // normalized vector of the view direction
    const double sensorDiag = 0.035; // diagonal of the sensor
    const double focalDist = 0.04;   // focal distance
};
