#pragma once

#include <array>
#include <memory>
#include <vector>
#include "util.h"
#include "bbox.h"
#include <glm/glm.hpp>

//Struct for atmosphere entities (fog, clouds etcs)
struct AtmosphereEntity
{
		
	AtmosphereEntity(glm::dvec3 position, glm::dvec3 size, glm::dvec3 color, double scatter) : pos(position), col(color), sc(scatter), bbox(position - .5*size, position + .5*size)
	{		
	}

	//returns the density at the specified position
	virtual double density(glm::dvec3& pos)
	{
		return 0;
	}

	glm::dvec3 pos = { 0, 0, 0 };
	glm::dvec3 col = { 0, 0, 0 };
	BoundingBox bbox;
	double sc = 0;
};

struct HeightFog : AtmosphereEntity
{

	HeightFog(glm::dvec3 position, glm::dvec3 size, glm::dvec3 color, double scatter) :AtmosphereEntity(position, size, color, scatter)
	{
	}

	//returns the density at the specified position
	virtual double density(glm::dvec3& pos)
	{
		return clamp(pos.y, 0, 1);
	}

};
