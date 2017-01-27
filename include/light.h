#pragma once

#include <array>
#include <memory>
#include <vector>
#include "util.h"
#include <glm/glm.hpp>

//Struct for all scene lights
struct Light
{

	std::vector<glm::dvec3> points;
	int count = 0;

	Light(glm::dvec3 position, glm::dvec3 color, double radius) : pos(position), col(color), rad(radius)
	{
		points.reserve(250);
		subrandUnitVec(points, 250);
	}

	//returns a point on the surface of the light source
	glm::dvec3 getPoint()
	{
		return pos + rad*points[rand() % points.size()];
	}

	glm::dvec3 getPoint(int i)
	{
		return pos + rad*points[rand() + 1097*i % points.size()];
	}

	glm::dvec3 pos = { 0, 0, 0 };
	glm::dvec3 col = { 0, 0, 0 };
	double rad = 0;
};
