#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "bbox.h"
#include "entities.h"
#include "octree.h"


Octree::Octree(glm::dvec3 min, glm::dvec3 max) : _root(Node({min, max})) {}


/// Store an entity in the correct position of the octree.
void Octree::push_back(Entity* object) {
        // TODO Implement this



	_root._entities.push_back(object);
}

/// Returns list of entities that have the possibility to be intersected by the ray.
std::vector<Entity*> Octree::intersect(const Ray& ray) const {
        // TODO Implement this
        return _root._entities;
}


Octree::Node::Node(const BoundingBox& bbox) : _bbox(bbox) {}

/// Subdivides the current node into 8 children.
void Octree::Node::partition()
{
	// TODO Implement this
};

bool Octree::Node::is_leaf() const { return _children[0] == nullptr; }

