#pragma once

#include <glm/glm.hpp>

struct Photon
{
	Photon(glm::dvec3 origin, glm::dvec3 dir) : origin(std::move(origin)), dir(glm::normalize(dir))
	{
		invDir = glm::dvec3(1.0 / dir);
	}

    glm::dvec3 origin;
    glm::dvec3 dir; // normalized directional vector
	glm::dvec3 invDir;
	glm::dvec3 col;
};
