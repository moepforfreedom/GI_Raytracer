#pragma once

#include <array>
#include <cassert>

#include <glm/glm.hpp>

/// Represents an axis-aligned bounding box.
struct BoundingBox {
    BoundingBox(glm::dvec3 min, glm::dvec3 max) : min(min), max(max) {
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
    bool intersect(const BoundingBox& other) const {
        // TODO Implement this
        return false;
    }

    /// Check if a point lies within the bounding box.
    bool contains(glm::dvec3 point) const {
        // TODO Implement this
        return false;
    }
};
