#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>

#include "bbox.h"
#include "material.h"
#include "ray.h"
#include "util.h"
#include <iostream>
#include "octree.h"

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
		phi = atan2(phit.y, phit.x);

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
			phi = atan2(phit.y, phit.x);

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

		 // Compute cone tangent vectors
		 glm::dvec3 dpdu(-phiMax * phit.y, phiMax * phit.x, 0);
		 glm::dvec3 dpdv(-phit.x / (1.f - v),- phit.y / (1.f - v), height);

		 normal = glm::normalize(glm::cross(dpdu, dpdv));

		 return true;
	}

	virtual BoundingBox boundingBox() const
	{
		glm::dvec3 verts[8];
		glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY); //TODO: change this somehow


		//vertices of bounding pyramid in object space
		verts[0] = rad*glm::dvec3(-1, -1, 0);
		verts[1] = rad*glm::dvec3(-1, 1, 0);
		verts[2] = rad*glm::dvec3(1, -1, 0);
		verts[3] = rad*glm::dvec3(1, 1, 0);
		verts[4] = glm::dvec3(0, 0, height);

		/*for (int j = 0; j < 4; j++)
			verts[4+j] = verts[j] + glm::dvec3(0, 0, height);*/

		//transform into world space and compute axis aligned bounding box TODO: find a more space efficient solution
		for (int i = 0; i < 5; i++)
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

struct vertex
{
	glm::dvec3 pos;
	glm::dvec3 norm = {0, 0, 0};
	glm::dvec2 texCoord;

	vertex(glm::dvec3 position, glm::dvec3 normal, glm::dvec2 uv)
	{
		pos = position;
		norm = glm::normalize(normal);
		texCoord = uv;
	}

	vertex(glm::dvec3 position, glm::dvec3 normal)
	{
		pos = position;
		norm = glm::normalize(normal);
	}

	vertex(glm::dvec3 position)
	{
		pos = position;
	}
};

struct triangle: Entity
{
	glm::dmat3x3 rot;
	vertex* vertices[3];
	glm::dvec3 norm;

	triangle(vertex* v1, vertex* v2, vertex* v3, const Material& material) : Entity(material)
	{
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;

		norm = glm::normalize(glm::cross((v2->pos - v1->pos), (v3->pos - v1->pos)));
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		double dot = glm::dot(ray.dir, norm);

		if (dot <= EPSILON && dot >= -EPSILON)
			return false;

		double t = glm::dot((vertices[0]->pos - ray.origin), norm) / dot;

		if (t < 0)
			return false;

		intersect = ray.origin + t*ray.dir;

		if (glm::dot(glm::cross((vertices[1]->pos - vertices[0]->pos), (intersect - vertices[0]->pos)), norm) < 0)
			return false;

		if (glm::dot(glm::cross((vertices[2]->pos - vertices[1]->pos), (intersect - vertices[1]->pos)), norm) < 0)
			return false;

		if (glm::dot(glm::cross((vertices[0]->pos - vertices[2]->pos), (intersect - vertices[2]->pos)), norm) < 0)
			return false;

		//Interpolate normal if vertex normals are set
		if (vecLengthSquared(vertices[0]->norm) > 0 && vecLengthSquared(vertices[1]->norm) > 0 && vecLengthSquared(vertices[2]->norm) > 0)
		{

			glm::dvec3 d1 = intersect - vertices[0]->pos;
			glm::dvec3 d2 = intersect - vertices[1]->pos;
			glm::dvec3 d3 = intersect - vertices[2]->pos;


			double area = glm::length(glm::cross(vertices[0]->pos - vertices[1]->pos, vertices[0]->pos - vertices[2]->pos));

			double a1 = glm::length(glm::cross(d2, d3)) / area;
			double a2 = glm::length(glm::cross(d3, d1)) / area;
			double a3 = glm::length(glm::cross(d1, d2)) / area;


			normal = a1*vertices[0]->norm + a2*vertices[1]->norm + a3*vertices[2]->norm;
		}
		else
			normal = norm;

		return true;

	}

	virtual BoundingBox boundingBox() const
	{
		glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY);
		for (int i = 0; i < 3; i++)
		{
			//vertices[i] = verts[i] * glm::inverse(rot) + pos;

			glm::dvec3 vpos = vertices[i]->pos;

			if (vpos.x < min.x)	min.x = vpos.x;
			if (vpos.x > max.x)	max.x = vpos.x;
			if (vpos.y < min.y)	min.y = vpos.y;
			if (vpos.y > max.y)	max.y = vpos.y;
			if (vpos.z < min.z)	min.z = vpos.z;
			if (vpos.z > max.z)	max.z = vpos.z;

      //special cases for axis aligned triangles
      if(max.x <= min.x) max.x += EPSILON;
      if(max.y <= min.y) max.y += EPSILON;
      if(max.z <= min.z) max.z += EPSILON;

      /*std::cout << "bbox min: " << min.x << ", " << min.y << ", " << min.z << "\n";
      std::cout << "bbox max: " << max.x << ", " << max.y << ", " << max.z << "\n";*/

		}

		return BoundingBox(min, max);
	}
};


//Mesh objects only generate triangles and add them to the scene, intersection/bbox tests are performed on the triangle objects
struct sphereMesh : Entity
{
	double rad;
	Octree* scene;
	int count = 0;

	sphereMesh(Octree* o, glm::dvec3 position, double radius, int subdivs, const Material& material) : Entity(material), rad(radius)
	{
		pos = position;
		scene = o;

		std::vector<std::array<glm::dvec3, 3>> tris;
		std::vector<std::array<glm::dvec3, 3>> tmp;

		//Idea: create octahedron, then subdivide and normalize to approximate a sphere
		tris.push_back({ glm::dvec3(-1, 0, 0), glm::dvec3(0, -1, 0) , glm::dvec3(0, 0, -1) });
		tris.push_back({ glm::dvec3(0, -1, 0), glm::dvec3(1, 0, 0) , glm::dvec3(0, 0, -1) });
		tris.push_back({ glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0) , glm::dvec3(0, 0, -1) });
		tris.push_back({ glm::dvec3(0, 1, 0), glm::dvec3(-1, 0, 0) , glm::dvec3(0, 0, -1) });

		tris.push_back({ glm::dvec3(-1, 0, 0), glm::dvec3(0, -1, 0) , glm::dvec3(0, 0, 1) });
		tris.push_back({ glm::dvec3(0, -1, 0), glm::dvec3(1, 0, 0) , glm::dvec3(0, 0, 1) });
		tris.push_back({ glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0) , glm::dvec3(0, 0, 1) });
		tris.push_back({ glm::dvec3(0, 1, 0), glm::dvec3(-1, 0, 0) , glm::dvec3(0, 0, 1) });

		for (int j = 0; j < subdivs; j++)
		{

			for (int i = 0; i < tris.size(); i++)
			{
				
				//subdivide each triangle into 4 triangles
				/**
				*      v2
				*     /   \
				*    a --- b
				*   /  \ /  \
				*  v1---c---v3
				*/

				glm::dvec3 v1 = glm::normalize(tris[i][0]);
				glm::dvec3 v2 = glm::normalize(tris[i][1]);
				glm::dvec3 v3 = glm::normalize(tris[i][2]);

				glm::dvec3 a = glm::normalize(glm::mix(v1, v2, .5));
				glm::dvec3 b = glm::normalize(glm::mix(v2, v3, .5));
				glm::dvec3 c = glm::normalize(glm::mix(v1, v3, .5));

				tmp.push_back({ v1, a, c });
				tmp.push_back({ a, v2, b });
				tmp.push_back({ a, b,  c });
				tmp.push_back({ c, b, v3 });				
			}

			tris = std::move(tmp);
		}

		for (int i = 0; i < tris.size(); i++)
		{
			scene->push_back(new triangle(new vertex(rad*tris[i][0] + pos, tris[i][0]), new vertex(rad*tris[i][1] + pos, tris[i][1]), new vertex(rad*tris[i][2] + pos, tris[i][2]), material));
			count++;
		}

		std::cout << "sphere tris: " << count << "\n";

	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(glm::dvec3(0, 0, 0), glm::dvec3(1, 1, 1));
	}
};

struct coneMesh : Entity
{
	double rad, height;
	Octree* scene;
	glm::dmat3x3 rot;
	int count = 0;

	coneMesh(Octree* o, glm::dvec3 position, glm::dvec3 rotation, double radius, double height, int tris, const Material& material) : Entity(material), rad(radius), height(height)
	{
		pos = position;
		scene = o;
		rot = glm::eulerAngleXYZ(rotation.x, rotation.y, rotation.z);

		glm::dmat3x3 baseRot = glm::eulerAngleXYZ(0.0, 0.0, 2 * M_PI / tris);

		vertex top(glm::dvec3(0, 0, height), glm::dvec3(0, 1, 0));
		vertex base(glm::dvec3(0, 0, 0), glm::dvec3(0, -1, 0));

		glm::dvec3 last(rad, 0, 0);

		for (int i = 0; i < tris; i++)
		{
			glm::dvec3 next = baseRot*last;

			scene->push_back(new triangle(new vertex(rot*last + pos, rot*last), new vertex(rot*next + pos, rot*next), new vertex(rot*glm::dvec3(0, 0, height) + pos, rot*last), material));
			scene->push_back(new triangle(new vertex(rot*last + pos, rot*glm::dvec3(0, 0, -1)), new vertex(rot*next + pos, rot*glm::dvec3(0, 0, -1)), new vertex(glm::dvec3(0, 0, 0) + pos, rot*glm::dvec3(0, 0, -1)), material));

			last = next;
			count += 2;
		}

	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		glm::dvec3 verts[8];
		glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY); //TODO: change this somehow


		//vertices of bounding pyramid in object space
		verts[0] = rad*glm::dvec3(-1, -1, 0);
		verts[1] = rad*glm::dvec3(-1, 1, 0);
		verts[2] = rad*glm::dvec3(1, -1, 0);
		verts[3] = rad*glm::dvec3(1, 1, 0);
		verts[4] = glm::dvec3(0, 0, height);

		/*for (int j = 0; j < 4; j++)
		verts[4+j] = verts[j] + glm::dvec3(0, 0, height);*/

		//transform into world space and compute axis aligned bounding box TODO: find a more space efficient solution
		for (int i = 0; i < 5; i++)
		{
			verts[i] = verts[i] * glm::inverse(rot) + pos;

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

struct quadMesh : Entity
{
	quadMesh(Octree* o, glm::dvec3 v1, glm::dvec3 v2, glm::dvec3 v3, glm::dvec3 v4, const Material& material)
	{
		o->push_back(new triangle(new vertex(v1), new vertex(v2), new vertex(v3), material));
		o->push_back(new triangle(new vertex(v3), new vertex(v2), new vertex(v4), material));
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(glm::dvec3(0, 0, 0), glm::dvec3(1, 1, 1));
	}
};

struct boxMesh : Entity
{
	glm::dmat3x3 rot;

	boxMesh(Octree* o, glm::dvec3 position, glm::dvec3 size, glm::dvec3 rotation, const Material& material)
	{
		pos = position;
		rot = glm::eulerAngleXYZ(rotation.x, rotation.y, rotation.z);

		std::vector<std::array<glm::dvec3, 3>> tris;

		tris.push_back({ glm::dvec3(-1, -1, -1), glm::dvec3(-1, 1, -1), glm::dvec3(1, -1, -1)  });
		tris.push_back({ glm::dvec3(-1, 1, -1) , glm::dvec3(1, 1, -1), glm::dvec3(1, -1, -1) });

		tris.push_back({ glm::dvec3(-1, -1, -1), glm::dvec3(-1, -1, 1) , glm::dvec3(-1, 1, -1) });
		tris.push_back({ glm::dvec3(-1, -1, 1), glm::dvec3(-1, 1, 1) , glm::dvec3(-1, 1, -1) });

		tris.push_back({ glm::dvec3(-1, -1, -1) ,glm::dvec3(1, -1, -1), glm::dvec3(-1, -1, 1) });
		tris.push_back({ glm::dvec3(-1, -1, 1) ,glm::dvec3(1, -1, -1), glm::dvec3(1, -1, 1) });

		tris.push_back({ glm::dvec3(-1, -1, 1), glm::dvec3(1, -1, 1), glm::dvec3(-1, 1, 1) });
		tris.push_back({ glm::dvec3(1, -1, 1), glm::dvec3(1, 1, 1), glm::dvec3(-1, 1, 1) });

		tris.push_back({ glm::dvec3(-1, 1, 1), glm::dvec3(1, 1, 1), glm::dvec3(1, 1, -1) });
		tris.push_back({ glm::dvec3(-1, 1, 1), glm::dvec3(1, 1, -1), glm::dvec3(-1, 1, -1) });

		tris.push_back({ glm::dvec3(1, -1, -1), glm::dvec3(1, 1, -1), glm::dvec3(1, -1, 1) });
		tris.push_back({ glm::dvec3(1, -1, 1), glm::dvec3(1, 1, -1), glm::dvec3(1, 1, 1) });


		for (int i = 0; i < 12; i++)
		{
			o->push_back(new triangle(new vertex(rot*(glm::normalize(tris[i][0])*size) + pos), new vertex(rot*(glm::normalize(tris[i][1])*size) + pos), new vertex(rot*(glm::normalize(tris[i][2])*size) + pos), material));
		}
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal) const
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(glm::dvec3(0, 0, 0), glm::dvec3(1, 1, 1));
	}
};

// TODO Implement implicit sphere - done
// TODO Implement implicit cone - done

// TODO Implement triangle - done
// TODO Implement explicit sphere (triangles) - done
// TODO Implement explicit quad (triangles) - done
// TODO Implement explicit cube (triangles) - done
