#pragma once

#include <glm/glm.hpp>

struct Photon
{
	Photon(glm::dvec3 origin, glm::dvec3 dir, glm::dvec3 color) : origin(std::move(origin)), dir(glm::normalize(dir)), col(color)
	{

	}

    glm::dvec3 origin;
    glm::dvec3 dir; // normalized directional vector
	glm::dvec3 col;
};
