#pragma once
#include "util.h"
#include <glm/glm.hpp>

struct Ray
{
	Ray(glm::dvec3 origin, glm::dvec3 d) : origin(origin), dir(d)
	{
		setDir(d);
	}

	void setDir(const glm::dvec3& d)
	{
		dir = glm::normalize(d);

		invDir = glm::dvec3(1.0 / dir);

		/*for (int i = 0; i < 24; i++)
		{
			r[i] = origin[i % 3];
			invD[i] = invDir[i % 3];
			invlz[i] = (invDir[i % 3] < 0.0);
		}*/
	}

	/*float r[24];
	float invD[24];
	float invlz[24];*/
    glm::dvec3 origin;
    glm::dvec3 dir; // normalized directional vector
	glm::dvec3 invDir;
};
