#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

//Struct for all scene lights
struct Light
{

	std::vector<glm::dvec3> points;
	int count = 0;

	Light(glm::dvec3 position, glm::dvec3 color, double radius) : pos(position), col(color), rad(radius)
	{
		points.reserve(2500);
		subrandUnitVec(points, 2500);
	}

	//returns a point on the surface if the light sourceS
	glm::dvec3 getPoint()
	{
		//std::cout << count << ": " << points[count % points.size()].x << ", " << points[count % points.size()].y << ", " << points[count % points.size()].z << "\n";
		count = (count+1) % points.size();
		return pos + rad*randomVec();//points[count];
	}

	glm::dvec3 pos = { 0, 0, 0 };
	glm::dvec3 col = { 0, 0, 0 };
	double rad = 0;
};
