#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "bbox.h"


struct Entity;

class Octree {
  public:

	Octree(glm::dvec3 min, glm::dvec3 max);

	std::vector<Entity*> tmp_entities;

    /// Store an entity in the correct position of the octree.
	void push_back(Entity* object);

    /// Returns list of entities that have the possibility to be intersected by the ray.
	std::vector<Entity*> intersect(const Ray& ray) const;

  private:
    struct Node {

		explicit Node(const BoundingBox& bbox);

		void partition();

		bool is_leaf() const;

		int depth;
        BoundingBox _bbox;
        std::vector<Entity*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

    Node _root;
};
