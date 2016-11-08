#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "bbox.h"
#include "entities.h"

class Octree {
  public:
    Octree(glm::dvec3 min, glm::dvec3 max) : _root(Node({min, max})) {}

    /// Store an entity in the correct position of the octree.
    void push_back(Entity* object) {
        // TODO Implement this
        _root._entities.push_back(object);
    }

    /// Returns list of entities that have the possibility to be intersected by the ray.
    std::vector<Entity*> intersect(const Ray& ray) const {
        // TODO Implement this
        return _root._entities;
    }

  private:
    struct Node {
        explicit Node(const BoundingBox& bbox) : _bbox(bbox) {}

        /// Subdivides the current node into 8 children.
        void partition(){
            // TODO Implement this
        };

        bool is_leaf() const { return _children[0] == nullptr; }

        BoundingBox _bbox;
        std::vector<Entity*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

    Node _root;
};
