#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

//Struct for all scene lights
struct Light
{
	Light(glm::dvec3 position, glm::dvec3 color, double radius) : pos(position), col(color), rad(radius)
	{

	}

	//returns a point on the surface if the light sourceS
	glm::dvec3 getPoint() const
	{
		return pos + rad*randomVec();
	}

	glm::dvec3 pos = { 0, 0, 0 };
	glm::dvec3 col = { 0, 0, 0 };
	double rad = 0;
};