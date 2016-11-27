#pragma once

#include <array>
#include <memory>
#include <vector>
#include <set>

#include <glm/glm.hpp>

#include "bbox.h"
#include "light.h"


struct Entity;

class Octree {
  public:
    struct Node {

		explicit Node(const BoundingBox& bbox);

		void partition();

		void debugVis(Node* root, Node* current);

		bool is_leaf() const;

		void intersect(const Ray& ray, std::vector<Entity*>& res, double tmin, double tmax) const;

		void Octree::Node::intersectSorted(const Ray& ray, std::vector<const Node*>& res, double tmin, double tmax) const;

        BoundingBox _bbox;
		double* mint = new double(0);
		double* maxt = new double(0);
        std::vector<Entity*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

public:

	Octree(glm::dvec3 min, glm::dvec3 max);

	std::vector<Light*> lights;

	/// Store an entity in the correct position of the octree.
	void push_back(Entity* object);
	void push_back(Light* object);

	void rebuild();

	/// Returns list of entities that have the possibility to be intersected by the ray.
	std::vector<Entity*> intersect(const Ray& ray, double tmin, double tmax) const;

	std::vector<const Octree::Node*> Octree::intersectSorted(const Ray& ray, double tmin, double tmax) const;

	bool valid;
    Node _root;
};
