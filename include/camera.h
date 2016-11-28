#pragma once

#include <glm/glm.hpp>
#define FOCAL_DIST 540

/// Represents the camera with information about the 'sensor' size.
struct Camera {
    explicit Camera(glm::dvec3 pos) : Camera(pos, {0, 0, 0}) {}
    Camera(glm::dvec3 pos, glm::dvec3 lookAt) : pos(pos), up({0, 0, 1.0}), forward(lookAt - pos) 
	{
        forward = glm::normalize(forward);
		up = glm::dvec3(0, 1.0, 0);
		right = glm::normalize(glm::cross(up, forward));
		up = glm::cross(forward, right);
    }

    glm::dvec3 pos;
    glm::dvec3 up;
    glm::dvec3 forward;              // normalized vector of the view direction
	glm::dvec3 right;
    const double sensorDiag = 0.035*FOCAL_DIST; // diagonal of the sensor
    const double focalDist = 0.04*FOCAL_DIST;   // focal distance
};
