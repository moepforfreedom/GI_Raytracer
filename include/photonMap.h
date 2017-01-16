#pragma once

#include <array>
#include <memory>
#include <vector>
#include <set>

#include <glm/glm.hpp>

#include "bbox.h"


struct Photon;

class PhotonMap {
  public:
    struct Node {

		explicit Node(const BoundingBox& bbox);

		void partition();

		bool is_leaf() const;

		void get(glm::dvec3& pos, std::vector<Photon*> res, double dist) const;

        BoundingBox _bbox;
        std::vector<Photon*> _entities;
        std::array<std::unique_ptr<Node>, 8> _children;
    };

public:

	PhotonMap(glm::dvec3 min, glm::dvec3 max);

	/// Store an entity in the correct position of the octree.
	void push_back(Photon* object);

	void rebuild();

	/// Returns list of entities that have the possibility to be intersected by the ray.
	std::vector<Photon*> getInRange(glm::dvec3& pos, double dist) const;

	bool valid;
    Node _root;
};