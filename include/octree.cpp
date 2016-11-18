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

}

void Octree::rebuild()
{
	if (_root._entities.size() > MAX_ENTITIES_PER_LEAF)
	{
		std::cout << "subdividing root...\n";
		_root.partition();

		_root.debugVis(&_root, &_root);
	}
				
}

void Octree::Node::debugVis(Node* root, Node* current)
{
	std::cout << "inserting debug vis object...\n";
	if (current->is_leaf())
	{
		//root->_entities.push_back(new triangle(new vertex(current->_bbox.min), new vertex(current->_bbox.min + glm::dvec3(current->_bbox.dx(), 0, 0)), new vertex(current->_bbox.max), Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));

		root->_entities.push_back(new sphere(current->_bbox.min + .0*current->_bbox.max, .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), 0, 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, current->_bbox.dy(), 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, 0, current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), current->_bbox.dy(), 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, current->_bbox.dy(), current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), 0, current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
		root->_entities.push_back(new sphere(current->_bbox.max, .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
	}
	else
	{
		for (int i = 0; i < 8; i++)
		{
			debugVis(root, current->_children[i].get());
		}
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
		if (_children[i]->_entities.size() > MAX_ENTITIES_PER_LEAF && _bbox.dx() > MIN_LEAF_SIZE)
		{
			std::cout << "subdividing node, size: " << _children[i]->_bbox.dx() << ", entities: " << _children[i]->_entities.size() << "\n";
			_children[i]->partition();

		}
	}
};

bool Octree::Node::is_leaf() const { return _children[0] == nullptr; }
