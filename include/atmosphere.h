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

	//returns the density at the specified position, defined as intersection probability per unit length
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
	double d;
	HeightFog(glm::dvec3 position, glm::dvec3 size, glm::dvec3 color, double density, double scatter, int noiseScale) : AtmosphereEntity(position, size, color, scatter), d(density)
	{
	}

	//returns the density at the specified position, defined as intersection probability per unit length
	virtual double density(glm::dvec3& p)
	{
		double ymax = pos.y + .5*bbox.dy();
		return d*fastPow((ymax - p.y)/bbox.dy(), 2);
	}

};
