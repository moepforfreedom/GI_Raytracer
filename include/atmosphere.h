#pragma once

#include <array>
#include <memory>
#include <vector>
#include "util.h"
#include "bbox.h"
#include <glm/glm.hpp>

//Struct for atmosphere entities (fog, clouds etc)
struct AtmosphereEntity
{
		
	AtmosphereEntity(glm::dvec3 position, glm::dvec3 size, glm::dvec3 color, double scatter) : pos(position), col(color), sc(scatter), bbox(position - .5*size, position + .5*size)
	{		
	}

	//returns the density at the specified position, defined as intersection probability per unit length
	virtual double density(const glm::dvec3& pos)
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
	int nscale;
	std::vector<double> noiseGrid;
	glm::dvec3 s;

	HeightFog(glm::dvec3 position, glm::dvec3 size, glm::dvec3 color, double density, double scatter, int noiseScale) : AtmosphereEntity(position, size, color, scatter), d(density), nscale(noiseScale), s(size)
	{
		noiseGrid.reserve((size.x + 1)*(size.y + 1)*(size.z + 1)*std::pow(noiseScale, 3));

		for (int i = 0; i < (size.x + 1)*(size.y + 1)*(size.z + 1)*std::pow(noiseScale, 3); i++)
		{
			noiseGrid.push_back(drand());
		}

		nscale = 1;
	}

	//returns the density at the specified position, defined as intersection probability per unit length
	virtual double density(const glm::dvec3& p)
	{
		double ymax = pos.y + .5*s.y;

		glm::dvec3 rel = (double)nscale*(p - bbox.min);

		//trilinear interpolation from grid values
		double dx = nscale*(rel.x - (int)rel.x);
		double dy = nscale*(rel.y - (int)rel.y);
		double dz = nscale*(rel.z - (int)rel.z);

		double c00 = (1 - dx)*noiseGrid[((int)rel.x * nscale*s.x + (int)rel.y) * nscale*s.z + (int)rel.z]
						 + dx*noiseGrid[(((int)rel.x + 1) * nscale*s.x + (int)rel.y) * nscale*s.z + (int)rel.z];

		double c01 = (1 - dx)*noiseGrid[((int)rel.x * nscale*s.x + (int)rel.y) * nscale*s.z + (int)rel.z + 1]
						 + dx*noiseGrid[(((int)rel.x + 1) * nscale*s.x + (int)rel.y) * nscale*s.z + (int)rel.z + 1];

		double c10 = (1 - dx)*noiseGrid[((int)rel.x * nscale*s.x + ((int)rel.y + 1)) * nscale*s.z + (int)rel.z]
						 + dx*noiseGrid[(((int)rel.x + 1) * nscale*s.x + ((int)rel.y + 1)) * nscale*s.z + (int)rel.z];

		double c11 = (1 - dx)*noiseGrid[((int)rel.x * nscale*s.x + ((int)rel.y + 1)) * nscale*s.z + (int)rel.z + 1]
						 + dx*noiseGrid[(((int)rel.x + 1) * nscale*s.x + ((int)rel.y + 1)) * nscale*s.z + (int)rel.z + 1];

		double c0 = c00*(1 - dy) + c10*dy;

		double c1 = c01*(1 - dy) + c11*dy;

		
		double noise = fastPow((1 - dz)*c0 + dz*c1, 7);

		return d*noise*fastPow((ymax - p.y) / s.y, 2);
	}

};
