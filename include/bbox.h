#pragma once

#include <array>
#include <cassert>

#include <glm/glm.hpp>

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
};
