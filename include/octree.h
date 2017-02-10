#pragma once

#include <array>
#include <memory>
#include <vector>
#include <set>

#include <glm/glm.hpp>

#include "bbox.h"
#include "light.h"
#include "atmosphere.h"


struct Entity;

class Octree {
  public:
    struct Node {

		explicit Node(const BoundingBox& bbox);

		void partition();

		void debugVis(Node* root, Node* current);

		bool is_leaf() const;

		void intersect(const Ray& ray, std::vector<Entity*>& res, double tmin, double tmax, float* tval0, float* tval1, int n) const;

		void intersectSorted(const Ray& ray, std::vector<const Node*>& res, double tmin, double tmax, float* tval0, float* tval1, int n) const;

        BoundingBox _bbox;
		double* mint = new double(0);
		double* maxt = new double(0);
		//float boxes[48];
        std::vector<Entity*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

public:

	Octree(glm::dvec3 min, glm::dvec3 max);

	std::vector<Light*> lights;
	std::vector<AtmosphereEntity*> at;

	/// Store an entity in the correct position of the octree.
	void push_back(Entity* object);
	void push_back(Light* object);
	void push_back(AtmosphereEntity* entity);

	void rebuild();

	/// Returns list of entities that have the possibility to be intersected by the ray.
	std::vector<Entity*> intersect(const Ray& ray, double tmin, double tmax) const;

	std::vector<const Octree::Node*> intersectSorted(const Ray& ray, double tmin, double tmax) const;

	// Returns density, color and scattering coefficient of the atmosphere at the specified position
	double atmosphereDensity(glm::dvec3& pos, glm::dvec3&col, double& scatter);

	bool atmosphereBounds(Ray r, double& mint, double& maxt);

	bool valid;
    Node _root;
};
