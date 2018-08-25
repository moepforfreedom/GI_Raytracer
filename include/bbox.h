#pragma once

#include <array>
#include <cassert>
#include <iostream>

#include <glm/glm.hpp>
#include "ray.h"
#include "util.h"

/// Represents an axis-aligned bounding box.
struct BoundingBox
{
    BoundingBox(glm::dvec3 min, glm::dvec3 max) : min(min), max(max)
	{
        assert(min.x <= max.x);
        assert(min.y <= max.y);
        assert(min.z <= max.z);
    }

	double dx() const { return max.x - min.x; }
    double dy() const { return max.y - min.y; }
    double dz() const { return max.z - min.z; }

	glm::dvec3 center() const { return min + 0.5*(max-min); }

	glm::dvec3 size() const { return max - min; }

    glm::dvec3 min;
    glm::dvec3 max;

    //Check if another bounding box intersects the current bounding box.
    bool intersect(const BoundingBox& other) const
	{
		return (min.x <= other.max.x && max.x >= other.min.x) &&
			(min.y <= other.max.y && max.y >= other.min.y) &&
			(min.z <= other.max.z && max.z >= other.min.z);
    }

    //Check if a point lies within the bounding box.
    inline bool contains(glm::dvec3 point) const
	{
		return point.x >= min.x && point.y >= min.y && point.z >= min.z && point.x < max.x && point.y < max.y && point.z < max.z;
    }

	//checks if a Ray intersects the bounding box, implementation based on Peter Shirleys algorithm
	inline bool intersect(const Ray& ray, double tmin, double tmax, double& toutmin, double& toutmax) const
	{
		for (int i = 0; i < 3; i++)
		{
            double t0 = ((&min.x)[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];
			double t1 = ((&max.x)[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];

			if ((&ray.invDir.x)[i] < 0.0)
			{
				double tmp = t0;
				t0 = t1;
				t1 = tmp;
			}

			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
			{
				toutmin = INFINITY;
				toutmax = -INFINITY;
				return false;
			}
		}
		toutmin = tmin;
		toutmax = tmax;
		return true;
	}

	inline bool intersectMulti(const Ray& ray, double tmin, double tmax, double& toutmin, double& toutmax, float* t0, float*t1) const
	{
		for (int i = 0; i < 3; i++)
		{

			/*if ((&ray.invDir.x)[i] < 0.0)
			{
				double tmp = t0[i];
				t0[i] = t1[i];
				t1[i] = tmp;
			}*/

			tmin = t0[i] > tmin ? t0[i] : tmin;
			tmax = t1[i] < tmax ? t1[i] : tmax;
			if (tmax < tmin)
			{
				toutmin = INFINITY;
				toutmax = -INFINITY;
				return false;
			}
		}

		//if (tmax > tmin)
		//{
			toutmin = tmin;
			toutmax = tmax;
			return true;
		//}
		//return false;
	}

	inline bool intersectSimpleMulti(const Ray& ray, double tmin, double tmax, float* t0, float*t1) const
	{
		for (int i = 0; i < 3; i++)
		{
			tmin = t0[i] > tmin ? t0[i] : tmin;// std::max(t0[i], tmin);
			tmax = t1[i] < tmax ? t1[i] : tmax;
		}
		return tmax > tmin;
	}

	//checks if a Ray intersects the bounding box, implementation based on Peter Shirleys algorithm
	inline bool intersectSimple(const Ray& ray, double tmin, double tmax) const
	{
		for (int i = 0; i < 3; i++)
		{
			double t0 = ((&min.x)[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];
			double t1 = ((&max.x)[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];

			if ((&ray.invDir.x)[i] < 0.0)
			{
				double tmp = t0;
				t0 = t1;
				t1 = tmp;
			}
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
			{
				return false;
			}
		}
		return true;
	}
};
