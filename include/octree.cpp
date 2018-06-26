#include <array>
#include <memory>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include <glm/glm.hpp>

#include "bbox.h"
#include "entities.h"
#include "atmosphere.h"
#include "octree.h"
#include "light.h"


Octree::Octree(glm::dvec3 min, glm::dvec3 max) : _root(Node({ min, max }))
{

}

int nodes = 0;
int skipped_subdiv = 0;

/// Store an entity in the root node
void Octree::push_back(Entity* object)
{
	_root._entities.push_back(object);
	valid = false;
}

//Lights are stored in a separate vector
void Octree::push_back(Light* light)
{
	lights.push_back(light);
	//_root._entities.push_back(new sphere(light->pos, light->rad - SHADOW_BIAS, Material(new texture(glm::dvec3(0, 0, 0)), new texture(5.0*light->col), 1, 1)));
	//valid = false;
}

void Octree::push_back(AtmosphereEntity* entity)
{
	at.push_back(entity);
}

void Octree::rebuild()
{
	std::cout << "total entities: " << _root._entities.size() << "\n";
	std::cout << "precomputing photon bounds...\n";

	glm::dvec3 avgPos(0, 0, 0);
	double count = 0;

	for(Entity* current : _root._entities)
	{
		if(current->material.roughness < 0.1)
		{
			BoundingBox bbox = current->boundingBox();
			avgPos += bbox.center();

			count++;
		}
	}

	if(count > 0)
		avgPos /= count;

	std::cout << "average pos: " << avgPos.x  << ", " << avgPos.y << ", " << avgPos.z << "\n";

	for(Light* l : lights)
	{
		double maxAngle = 0;

		l->dir = glm::normalize(avgPos - l->pos);

		for(Entity* current : _root._entities)
		{
			if(current->material.roughness < 0.1)
			{
				BoundingBox bbox = current->boundingBox();

				avgPos += bbox.center();

				double angle = 1.0 - acos(glm::dot(l->dir, glm::normalize(l->pos - bbox.min)))/M_PI;

				maxAngle = std::max(maxAngle, angle);
				count++;
			}
		}

		std::cout << "max photon angle: " << maxAngle << "\n";
		l->angle = maxAngle;
	}

	std::cout << "rebuilding scene octree...\n";

	if (_root._entities.size() > MAX_ENTITIES_PER_LEAF)
	{
		std::cout << "subdividing root...\n";
		_root.partition();

		//_root.debugVis(&_root, &_root);
	}

	std::cout << "done, total nodes: " << nodes << "\n";
	std::cout << "skipped: " << skipped_subdiv << "\n";

	valid = true;

}

//Creates a debug view where the boundong boxes of leaf nodes are visualized using spheres
void Octree::Node::debugVis(Node* root, Node* current)
{
	std::cout << "inserting debug vis object...\n";
	if (current->is_leaf())
	{
		if (current->_entities.size() > 0)
		{

			/*root->_entities.push_back(new sphere(current->_bbox.min + .0*current->_bbox.max, .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), 0, 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, current->_bbox.dy(), 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, 0, current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), current->_bbox.dy(), 0), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(0, current->_bbox.dy(), current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.min + glm::dvec3(current->_bbox.dx(), 0, current->_bbox.dz()), .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));
			root->_entities.push_back(new sphere(current->_bbox.max, .125, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 1, 0))));*/
		}
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
std::vector<Entity*> Octree::intersect(const Ray& ray, double tmin, double tmax) const
{
	std::vector<Entity*> res;
	res.reserve(256);

	//return _root._entities;

	float t0[24];
	float t1[24];

	for (int i = 0; i < 3; i++)
	{
		t0[i] = ((&(_root._bbox.min.x))[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];
		t1[i] = ((&(_root._bbox.max.x))[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];

		if ((&ray.invDir.x)[i] < 0.0)
		{
			double tmp = t0[i];
			t0[i] = t1[i];
			t1[i] = tmp;
		}
	}

	_root.intersect(ray, res, tmin, tmax, t0, t1, 0);


	/*//remove duplicates if lots of objects are returned
	if (res.size() > DUPLICATE_THRESHOLD)
	{
		std::sort(res.begin(), res.end());

		res.erase(std::unique(res.begin(), res.end()), res.end());
	}*/

	return res;
}

//Returns a sorted list of nodes that intersect the given ray
std::vector<std::pair<const Octree::Node*, double>> Octree::intersectSorted(const Ray& ray, double tmin, double tmax) const
{
	std::vector<std::pair<const Node*, double>> res;

	/*float t0[24];
	float t1[24];

	for (int i = 0; i < 3; i++)
	{
		t0[i] = ((&(_root._bbox.min.x))[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];
		t1[i] = ((&(_root._bbox.max.x))[i] - (&ray.origin.x)[i]) * (&ray.invDir.x)[i];

		if ((&ray.invDir.x)[i] < 0.0)
		{
			double tmp = t0[i];
			t0[i] = t1[i];
			t1[i] = tmp;
		}
	}*/

	_root.intersectSorted(ray, res, tmin, tmax);

	return res;
}

//returns the atmosphere density, color and scattering coefficient at the specified position
double Octree::atmosphereDensity(const glm::dvec3& pos, glm::dvec3&col, double& scatter)
{
	double d = 0;
	for (AtmosphereEntity* current : at)
	{
		if (current->bbox.contains(pos))
		{
			col = current->col;
			d += RAYMARCH_STEPSIZE*current->density(pos);
		}
	}
	return d;
}

//returns true if the ray intersects an atmosphere entity and determines the distance bounds
bool Octree::atmosphereBounds(const Ray& r, double& mint,  double& maxt)
{
	double min = 0;
	double max = 0;
	bool intersected = false;

	for (AtmosphereEntity* current : at)
	{
		double tmpmin = 0;
		double tmpmax = 0;

		if (current->bbox.intersect(r, mint, maxt, tmpmin, tmpmax))
		{
			min = std::min(min, tmpmin);
			max = std::max(max, tmpmax);
			intersected = true;
		}
	}

	mint = std::max(mint, min);
	maxt = std::min(maxt, max);
	return intersected;
}

Octree::Node::Node(const BoundingBox& bbox) : _bbox(bbox) {}

//returns a list of objects within a node that can potentially intersect the given ray
void Octree::Node::intersect(const Ray& ray, std::vector<Entity*>& res, double tmin, double tmax, float* tval0, float* tval1, int n) const
{
	float ta[24];
	float tb[24];

	if (_bbox.intersectSimple(ray, tmin, tmax)/*_bbox.intersectSimpleMulti(ray, tmin, tmax, &(tval0[3 * n]), &(tval1[3 * n]))*/)
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
			//intersectSIMD(ta, tb, boxes, ray.r, ray.invD, ray.invlz);
			for (int i = 0; i < 8; i++)
			{
				_children[i]->intersect(ray, res, tmin, tmax, ta, tb, i);

			}
		}
	}

}

//inserts the node into a sorted list if it intersects the given ray
void Octree::Node::intersectSorted(const Ray& ray, std::vector<std::pair<const Node*, double>>& res, double tmin, double tmax) const
{
	double t0, t1;
	/*float ta[24];
	float tb[24];*/
	if (_bbox.intersect(ray, tmin, tmax, t0, t1)/*_bbox.intersectMulti(ray, tmin, tmax, t0, t1, &(tval0[3*n]), &(tval1[3*n]))*/)
	{
		if (is_leaf())
		{
			//std::cout << t0 << "\n";
			if (_entities.size() > 0)
			{
				std::vector<std::pair<const Node*, double>>::iterator it = std::partition_point(res.begin(), res.end(),
																							    [&t0](std::pair<const Node*, double> n){return t0 >= n.second;});

				res.insert(it, { this , t0});
			}
		}
		else
		{
			//intersectSIMD(ta, tb, boxes, ray.r, ray.invD, ray.invlz);
			for (int i = 0; i < 8; i++)
			{
				_children[i]->intersectSorted(ray, res, tmin, tmax);

			}
		}
	}
}

/// Subdivides the current node into 8 children.
void Octree::Node::partition()
{
	glm::dvec3 mid = glm::mix(_bbox.min, _bbox.max, .5);
	double avg_entities = 0;

	_children[0] = std::make_unique<Node>(BoundingBox(_bbox.min, mid));
	_children[1] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z)));
	_children[2] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x , _bbox.min.y, _bbox.min.z+ .5*_bbox.dz()), glm::dvec3(mid.x, mid.y, mid.z + .5*_bbox.dz())));
	_children[3] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y, _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y, mid.z + .5*_bbox.dz())));
	_children[4] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z)));
	_children[5] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x + .5*_bbox.dx(), _bbox.min.y + .5*_bbox.dy(), _bbox.min.z), glm::dvec3(mid.x + .5*_bbox.dx(), mid.y + .5*_bbox.dy(), mid.z)));
	_children[6] = std::make_unique<Node>(BoundingBox(glm::dvec3(_bbox.min.x, _bbox.min.y + .5*_bbox.dy(), _bbox.min.z + .5*_bbox.dz()), glm::dvec3(mid.x, mid.y + .5*_bbox.dy(), mid.z + .5*_bbox.dz())));
	_children[7] = std::make_unique<Node>(BoundingBox(mid, _bbox.max));

	std::vector<Entity*>::iterator it = _entities.begin();

	while (it != _entities.end())
	{
		Entity* current = *it;

		for (int i = 0; i < 8; i++)
		{
			if (_children[i]->_bbox.intersect(current->boundingBox()) && current->intersect(_children[i]->_bbox) && current->boundingBox().dx() > EPSILON)
			{
				_children[i]->_entities.push_back(current);
			}
		}
		++it;
	}

	for (int i = 0; i < 8; i++)
	{
		avg_entities += (double)_children[i]->_entities.size();
	}

	avg_entities /= 8;

	//stop subdividing if the last subdivision didn't improve the entity counts
	if (avg_entities > MAX_SUBDIV_RATIO * _entities.size())
	{
		skipped_subdiv++;
		_entities.clear();
		return;
	}

	_entities.clear();
	_entities.shrink_to_fit();

	for (int i = 0; i < 8; i++)
	{
		if (_children[i]->_entities.size() > MAX_ENTITIES_PER_LEAF && _children[i]->_bbox.dx() > MIN_LEAF_SIZE)// && avg_entities < MAX_SUBDIV_RATIO * _entities.size())
		{
			//std::cout << "subdividing node, size: " << _children[i]->_bbox.dx() << ", entities: " << _children[i]->_entities.size() << "\n";
			_children[i]->partition();
		}
	}

	nodes += 8;
};

bool Octree::Node::is_leaf() const { return _children[0] == nullptr; }
