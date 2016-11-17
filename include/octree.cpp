#pragma once

#include <array>
#include <memory>
#include <vector>

#include <glm/glm.hpp>

#include "bbox.h"
#include "entities.h"
#include "octree.h"


Octree::Octree(glm::dvec3 min, glm::dvec3 max) : _root(Node({ min, max }))
{

}


/// Store an entity in the correct position of the octree.
void Octree::push_back(Entity* object) {
        // TODO Implement this


	_root._entities.push_back(object);

	if (_root._entities.size() > 870)
	{
		std::cout << "subdividing root...\n";
		/*_root.partition();

		Node* current = &_root;

		while (!current->is_leaf())
		{

			for (int i = 0; i < 8; i++)
			{
				glm::dvec3 minPos = current->_children[i]->_bbox.min;
				glm::dvec3 maxPos = current->_children[i]->_bbox.max;

				_root._entities.push_back(new sphere(maxPos, .5, Material(glm::dvec3(1*(1.0/i), 0, 0), glm::dvec3(1, 0, 0))));
			}

			current = current->_children[0].get();
		}*/
	}
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
	glm::dvec3 mid = glm::mix(_bbox.min, _bbox.max, .5);

	std::cout << "creating child nodes...\n";

	_children[0] = std::unique_ptr<Node>(new Node(BoundingBox(_bbox.min, mid)));
	_children[1] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z))));
	_children[2] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x , _bbox.min.y, _bbox.min.z+ .5*_bbox.dz()), glm::dvec3(mid.x, mid.y, mid.z + .5*_bbox.dz()))));
	_children[3] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z + .5*_bbox.dz()))));
	_children[4] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z))));
	_children[5] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y + .5*_bbox.dy(), mid.z))));
	_children[6] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z + .5*_bbox.dz()))));
	_children[7] = std::unique_ptr<Node>(new Node(BoundingBox(mid, _bbox.max)));

	std::vector<Entity*>::iterator it = _entities.begin();

  std::cout << "adding entities...\n";
	while (it != _entities.end())
	{
		Entity* current = *it;

		for (int i = 0; i < 8; i++)
		{
			if (_children[i]->_bbox.intersect(current->boundingBox()) && current->boundingBox().dx() > 0.0001)
			{
				_children[i]->_entities.push_back(current);
			}
		}
		it++;
	}

	for (int i = 0; i < 8; i++)
	{
		if (_children[i]->_entities.size() > 10 && _bbox.dx() > .5)
		{
			std::cout << "subdividing node, size: " << _children[i]->_bbox.dx() << ", entities: " << _children[i]->_entities.size() << "\n";
			_children[i]->partition();

		}
	}
};

bool Octree::Node::is_leaf() const { return _children[0] == nullptr; }
