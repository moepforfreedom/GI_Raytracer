#pragma once

#include <array>
#include <memory>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include <glm/glm.hpp>

#include "bbox.h"
#include "photonMap.h"
#include "photon.h"


PhotonMap::PhotonMap(glm::dvec3 min, glm::dvec3 max) : _root(Node({ min, max }))
{

}

/// Store an entity in the root node
void PhotonMap::push_back(Photon* object)
{
	_root._entities.push_back(object);
	valid = false;
}

void PhotonMap::rebuild()
{
	std::cout << "total entities: " << _root._entities.size() << "\n";
	std::cout << "rebuilding photon map...\n";
	if (_root._entities.size() > MAX_ENTITIES_PER_LEAF)
	{
		std::cout << "subdividing root...\n";
		_root.partition();
	}

	valid = true;

}

/// Returns list of photons in range of the given hit point
std::vector<Photon*> PhotonMap::getInRange(glm::dvec3& pos, double dist) const
{
	std::vector<Photon*> res;
	res.reserve(256);

	//return _root._entities;

	_root.get(pos, res, dist);

	return res;
}

PhotonMap::Node::Node(const BoundingBox& bbox) : _bbox(bbox) {}

//Returns list of photons in range of the given hit point
void PhotonMap::Node::get(glm::dvec3& pos, std::vector<Photon*> res, double dist) const
{
	if (is_leaf())
	{
		if (_entities.size() > 0)
		{
			res.insert(res.end(), _entities.begin(), _entities.end());
		}
	}
	else
	{
		int i = 0;

		while (!_children[i]->_bbox.contains(pos))
			i++;
		
		_children[i]->get(pos, res, dist);
	}

}

/// Subdivides the current node into 8 children.
void PhotonMap::Node::partition()
{
	glm::dvec3 mid = glm::mix(_bbox.min, _bbox.max, .5);
	double avg_entities = 0;

	_children[0] = std::unique_ptr<Node>(new Node(BoundingBox(_bbox.min, mid)));
	_children[1] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z))));
	_children[2] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x , _bbox.min.y, _bbox.min.z+ .5*_bbox.dz()), glm::dvec3(mid.x, mid.y, mid.z + .5*_bbox.dz()))));
	_children[3] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z + .5*_bbox.dz()))));
	_children[4] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z))));
	_children[5] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y + .5*_bbox.dy(), mid.z))));
	_children[6] = std::unique_ptr<Node>(new Node(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z + .5*_bbox.dz()))));
	_children[7] = std::unique_ptr<Node>(new Node(BoundingBox(mid, _bbox.max)));

	std::vector<Photon*>::iterator it = _entities.begin();

	while (it != _entities.end())
	{
		Photon* current = *it;

		for (int i = 0; i < 8; i++)
		{
			if (_children[i]->_bbox.contains(current->origin))
			{
				_children[i]->_entities.push_back(current);
			}
		}
		it++;
	}

	for (int i = 0; i < 8; i++)
	{
		avg_entities += (double)_children[i]->_entities.size();
	}

	avg_entities /= 8;

	//stop subdividing if the last subdivision didn't improve the entity counts
	if (avg_entities > MAX_SUBDIV_RATIO * _entities.size())
	{
		_entities.clear();
		return;
	}

	_entities.clear();

	for (int i = 0; i < 8; i++)
	{
		if (_children[i]->_entities.size() > MAX_ENTITIES_PER_LEAF && _children[i]->_bbox.dx() > MIN_LEAF_SIZE)
		{
			//std::cout << "subdividing node, size: " << _children[i]->_bbox.dx() << ", entities: " << _children[i]->_entities.size() << "\n";
			_children[i]->partition();
		}
	}
};

bool PhotonMap::Node::is_leaf() const { return _children[0] == nullptr; }
