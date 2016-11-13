#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "bbox.h"
#include "material.h"
#include "ray.h"
#include "util.h"
#include <iostream>

/// A base class for all entities in the scene.
struct Entity 
{

    constexpr Entity() : material(Material(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0))) {}
    constexpr Entity(const Material& material) : material(material) {}

    /// Check if a ray intersects the object. The arguments intersect and normal will contain the
    /// point of intersection and its normals.
    virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const = 0;

    /// Returns an axis-aligned bounding box of the entity.
    virtual BoundingBox boundingBox() const = 0;

    glm::dvec3 pos = {0, 0, 0};
	glm::dvec3 rot = {0, 0, 0};
    Material material;
};

struct sphere: Entity
{
	double rad;

	sphere(glm::dvec3 position, double radius, const Material& material) : Entity(material), rad(radius)
	{
		pos = position;		
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{

		double dot = glm::dot(ray.dir, (ray.origin - pos));

		double r = (pow(dot, 2) - pow(glm::length(ray.origin - pos), 2) + pow(rad, 2));

		if (r < 0)
			return false;

		double t_1 = -1*dot - sqrt(r);
		double t_2 = -1*dot + sqrt(r);

		//std::cout << t_1 << ", " << t_2 << ": ";

		if (t_1 < 0 && t_2 < 0)
			return false;
		else
		{
			if ((t_1 < t_2 && t_1 > 0) || t_2 < 0)
				intersect = ray.origin + ray.dir*t_1;
			else
				intersect = ray.origin + ray.dir*t_2;

			normal = glm::normalize(intersect - pos);
			//std::cout << "intersection\n";
			return true;
		}
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(pos + rad*glm::dvec3(-1, -1, -1), (pos + rad*glm::dvec3(1, 1, 1)));
	}
};

struct cone : Entity
{
	double rad, height, posX;

	glm::dmat3x3 rot;

	cone(glm::dvec3 position, glm::dvec3 rotation, double radius, double height, const Material& material) : Entity(material), rad(radius), height(height)
	{

		pos = position;		

		rot = glm::inverse(glm::eulerAngleXYZ(rotation.x, rotation.y, rotation.z)); 
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		//instead of transforming the cylinder, apply the inverse transform to the ray
		glm::dvec3 tmpOrigin = ray.origin*rot;
		glm::dvec3 tmpDir = ray.dir*rot;

		glm::dvec3 origin = (ray.origin - pos)*rot;
		glm::dvec3 dir = ray.dir*rot;

		double t_1, t_2, thit, phi;
		double phiMax = 2*M_PI; //can be used for radial clipping
		glm::dvec3 phit;

		//coefficients for quadratic equation
		double k = pow(rad / height, 2);

		double A = dir.x * dir.x + dir.y * dir.y -
				   k * dir.z * dir.z;
		double B = 2 * (dir.x * origin.x + dir.y * origin.y -
		          k * dir.z * (origin.z - height));
		double C = origin.x * origin.x + origin.y * origin.y -
				  k * (origin.z - height) * (origin.z - height);

		double discrim = B * B - 4.f * A * C;

		if (discrim < 0)
			return false;

		double rootDiscrim = sqrt(discrim);

		// Compute quadratic t values
		double q;

		if (B < 0)
			q = -.5f * (B - rootDiscrim);
		else
			q = -.5f * (B + rootDiscrim);

		t_1 = q / A;
		t_2 = C / q;

		//both intersection points are behind ray origin
		if (t_1 < 0 && t_2 < 0)
			return false;

		if (t_1 > t_2)
		{
			double tmp = t_1;
			t_1 = t_2;
			t_2 = tmp;
		}
		
		thit = t_1;

		//one intersection point is behind ray origin
		if (t_1 < 0)
			thit = t_2;
		else if (t_2 < 0)
			thit = t_1;	

		phit = origin + dir*thit;
		phi = atan2f(phit.y, phit.x);

		if (phi < 0.)
			phi += 2.f*M_PI;

		// Test cone intersection against clipping plane
		if (phit.z < 0 || phit.z > height || phi > phiMax) 
		 {
			 if (thit == t_2) 
				 return false;

			 thit = t_2;

			// Compute cone inverse mapping
			phit = origin + dir*thit;
			phi = atan2f(phit.y, phit.x);

			if (phi < 0.)
				phi += 2.f*M_PI;

			if (phit.z < 0 || phit.z > height || phi > phiMax)
				return false;
		}

		 //std::cout << t_1 << ", " << t_2 << "\n";

		 intersect = ray.origin + thit*ray.dir;

		 // Find parametric representation of cone hit
		 double u = phi / phiMax;
		 double v = phit.z / height;

		 // Compute cone $\dpdu$ and $\dpdv$
		 glm::dvec3 dpdu(-phiMax * phit.y, phiMax * phit.x, 0);
		 glm::dvec3 dpdv(-phit.x / (1.f - v),- phit.y / (1.f - v), height);

		 // Compute cone $\dndu$ and $\dndv$
		 glm::dvec3 d2Pduu = -phiMax * phiMax * glm::dvec3(phit.x, phit.y, 0.);
		 glm::dvec3 d2Pduv = phiMax / (1.f - v) * glm::dvec3(-phit.y, -phit.x, 0.);
		 glm::dvec3 d2Pdvv(0, 0, 0);

		 normal = glm::normalize(glm::cross(dpdu, dpdv));

		 return true;
	}

	virtual BoundingBox boundingBox() const
	{
		glm::dvec3 verts[8];
		glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY); //TODO: change this somehow

		verts[0] = rad*glm::dvec3(-1, -1, 0);
		verts[1] = rad*glm::dvec3(-1, 1, 0);
		verts[2] = rad*glm::dvec3(1, -1, 0);
		verts[3] = rad*glm::dvec3(1, 1, 0);

		for (int j = 0; j < 4; j++)
			verts[4+j] = verts[j] + glm::dvec3(0, 0, height);

		for (int i = 0; i < 8; i++)
		{
			verts[i] = verts[i]*glm::inverse(rot) + pos;

			if (verts[i].x < min.x)
				min.x = verts[i].x;
			if (verts[i].x > max.x)
				max.x = verts[i].x;

			if (verts[i].y < min.y)
				min.y = verts[i].y;
			if (verts[i].y > max.y)
				max.y = verts[i].y;

			if (verts[i].z < min.z)
				min.z = verts[i].z;
			if (verts[i].z > max.z)
				max.z = verts[i].z;
		}

		return BoundingBox(min, max);

	}
};

// TODO Implement implicit sphere
// TODO Implement implicit cone

// TODO Implement triangle
// TODO Implement explicit sphere (triangles)
// TODO Implement explicit quad (triangles)
// TODO Implement explicit cube (triangles)
