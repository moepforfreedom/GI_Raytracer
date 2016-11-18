#pragma once

#include <array>
#include <cassert>

#include <glm/glm.hpp>
#include "ray.h"

/// Represents an axis-aligned bounding box.
struct BoundingBox 
{
    BoundingBox(glm::dvec3 min, glm::dvec3 max) : min(min), max(max) 
	{
        assert(min.x < max.x);
        assert(min.y < max.y);
        assert(min.z < max.z);
    }

    double dx() const { return max.x - min.x; }
    double dy() const { return max.y - min.y; }
    double dz() const { return max.z - min.z; }

    const glm::dvec3 min;
    const glm::dvec3 max;

    /// Check if another bounding box intersects the current bounding box.
    bool intersect(const BoundingBox& other) const 
	{
		/*if (max.x < other.min.x) return false; // a is left of b
		if (min.x > other.max.x) return false; // a is right of b
		if (max.y < other.min.y) return false; // a is above b
		if (min.y > other.max.y) return false; // a is below b
		if (max.z < other.min.z) return false; // a is in front of b
		if (min.z > other.max.z) return false; // a is behind b
		return true; // boxes overlap*/

		return (min.x <= other.max.x && max.x >= other.min.x) &&
			(min.y <= other.max.y && max.y >= other.min.y) &&
			(min.z <= other.max.z && max.z >= other.min.z);
    }

    /// Check if a point lies within the bounding box.
    bool contains(glm::dvec3 point) const 
	{
		if ((point - min).x > 0 && (point - min).y > 0 && (point - min).z > 0 && (point - max).x > 0 && (point - max).y > 0 && (point - max).z > 0)
			return true;
		else
			return false;
    }

	//checks if a Ray intersects the bounding box, implementation based on Peter Shirleys algorithm
	inline bool intersect(const Ray& ray) const 
	{
		double tmin = 0;
		double tmax = INFINITY;
		//store vectors in arrays to avoid code duplication
		double minPos[3], maxPos[3], rayOrigin[3], rayDir[3];

		minPos[0] = min.x;
		minPos[1] = min.y;
		minPos[2] = min.z;

		maxPos[0] = max.x;
		maxPos[1] = max.y;
		maxPos[2] = max.z;

		rayOrigin[0] = ray.origin.x;
		rayOrigin[1] = ray.origin.y;
		rayOrigin[2] = ray.origin.z;

		rayDir[0] = ray.dir.x;
		rayDir[1] = ray.dir.y;
		rayDir[2] = ray.dir.z;

		for (int i = 0; i < 3; i++) 
		{
			double invD = 1.0f / rayDir[i];
			double t0 = (minPos[i] - rayOrigin[i]) * invD;
			double t1 = (maxPos[i] - rayOrigin[i]) * invD;

			if (invD < 0.0f) 
			{
				double tmp = t0;
				t0 = t1;
				t1 = tmp;
			}
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
				return false;
		}
		return true;
	}
};
