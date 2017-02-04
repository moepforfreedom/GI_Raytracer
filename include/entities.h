#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/component_wise.hpp>

#include "bbox.h"
#include "material.h"
#include "ray.h"
#include <iostream>
#include "octree.h"
//#include "util.h"

/// A base class for all entities in the scene.
struct Entity
{
    Entity() : material(Material(new texture(glm::dvec3(1, 0, 0)), new texture(glm::dvec3(0, 0, 0)), .75, 1)) {}
    Entity(const Material& material) : material(material) {}

    /// Check if a ray intersects the object. The arguments intersect and normal will contain the
    /// point of intersection and its normals.
    virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv) = 0;

	virtual bool intersect(const Ray& ray, double& t)
	{
		glm::dvec3 hit, norm;
        glm::dvec2 uv;

		bool i = intersect(ray, hit, norm, uv);

		t = glm::length(hit - ray.origin);

		return i;
	}

	virtual inline bool intersect(BoundingBox bbox)
	{
		return true;
	}

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

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
	{

		double dot = glm::dot(ray.dir, (ray.origin - pos));

		double r = (pow(dot, 2) - vecLengthSquared(ray.origin - pos) + pow(rad, 2));

		if (r < 0)
			return false;

		double sr = sqrt(r);

		double t_1 = -1*dot - sr;
		double t_2 = -1*dot + sr;


		if (t_1 < 0 && t_2 < 0)
			return false;
		else
		{
			if ((t_1 < t_2 && t_1 > 0) || t_2 < 0)
				intersect = ray.origin + ray.dir*t_1;
			else
				intersect = ray.origin + ray.dir*t_2;

			glm::dvec3 hitNorm = glm::normalize(intersect - pos);

			/*if (glm::dot(hitNorm, ray.dir) > 0)
				hitNorm = -1.0*hitNorm;*/

			normal = hitNorm;

			glm::dvec3 d = (pos - intersect)/rad;

			double v = .5 + asin(d.y) / M_PI;
            double u = .5 + atan2(d.z, d.x) / (2 * M_PI);

            uv = glm::dvec2(u, v);
			//std::cout << "intersection\n";
			return true;
		}
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(pos + rad*glm::dvec3(-1, -1, -1), (pos + rad*glm::dvec3(1, 1, 1)));
	}

	virtual inline bool intersect(BoundingBox bbox)
	{
		auto check = [&](
			const double pn,
			const double bmin,
			const double bmax) -> double
		{
			double out = 0;
			double v = pn;

			if (v < bmin)
			{
				double val = (bmin - v);
				out += val * val;
			}

			if (v > bmax)
			{
				double val = (v - bmax);
				out += val * val;
			}

			return out;
		};

		// Squared distance
		double sq = 0.0;

		sq += check(pos.x, bbox.min.x, bbox.max.x);
		sq += check(pos.y, bbox.min.y, bbox.max.y);
		sq += check(pos.z, bbox.min.z, bbox.max.z);

		return sq <= (rad*rad);
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

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
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
	vertex* vertices[3];
	glm::dvec3 edges[3];
	glm::dvec3 norm;
	double inv_area;

	triangle(vertex* v1, vertex* v2, vertex* v3, const Material& material) : Entity(material)
	{
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;

		edges[0] = vertices[1]->pos - vertices[0]->pos;
		edges[1] = vertices[2]->pos - vertices[1]->pos;
		edges[2] = vertices[0]->pos - vertices[2]->pos;

		norm = glm::normalize(glm::cross((v2->pos - v1->pos), (v3->pos - v1->pos)));

		inv_area = 1.0 / glm::length(glm::cross(vertices[0]->pos - vertices[1]->pos, vertices[0]->pos - vertices[2]->pos));
	}

	//Möller–Trumbore intersection test, based on https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	inline bool intersectMT(const glm::dvec3&   V1,  // Triangle vertices
					 const glm::dvec3&   V2,
					 const glm::dvec3&   V3,
					 const glm::dvec3&    O,  //Ray origin
					 const glm::dvec3&    D,  //Ray direction
					 double& tout)
	{
		glm::dvec3 e1, e2;  //Edge1, Edge2
		glm::dvec3 P, Q, T;
		float det, inv_det, u, v;
		float t;

		//Find vectors for two edges sharing V1
		//e1 = V2 - V1;
		e2 = V3 - V1;
		//Begin calculating determinant - also used to calculate u parameter
		P = glm::cross(D, e2);
		//if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
		det = glm::dot(edges[0], P);
		//NOT CULLING
		if (det > -EPSILON && det < EPSILON) return 0;
		inv_det = 1.f / det;

		//calculate distance from V1 to ray origin
		T = O - V1;

		//Calculate u parameter and test bound
		u = glm::dot(T, P) * inv_det;
		//The intersection lies outside of the triangle
		if (u < 0.f || u > 1.f) return 0;

		//Prepare to test v parameter
		Q = glm::cross(T, edges[0]);

		//Calculate V parameter and test bound
		v = glm::dot(D, Q) * inv_det;
		//The intersection lies outside of the triangle
		if (v < 0.f || u + v  > 1.f) return 0;

		t = glm::dot(e2, Q) * inv_det;

		if (t > EPSILON)
		{ 
			//ray intersection
			tout = t;
			return 1;
		}

		// No hit, no win
		return 0;
	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
	{		
		glm::dvec3 d1, d2, d3;

		double dot = glm::dot(ray.dir, norm);

		if (dot <= EPSILON && dot >= -EPSILON)
			return false;

		double t = glm::dot((vertices[0]->pos - ray.origin), norm) / dot;

		if (t < 0)
			return false;

		intersect = ray.origin + t*ray.dir;

		d1 = intersect - vertices[0]->pos;
		if (glm::dot(glm::cross(edges[0], d1), norm) < 0) return false;

		d2 = (intersect - vertices[1]->pos);
		if (glm::dot(glm::cross(edges[1], d2), norm) < 0) return false;

		d3 = (intersect - vertices[2]->pos);
		if (glm::dot(glm::cross(edges[2], d3), norm) < 0) return false;

		/*if (dot > 0)
			hitNorm = -1.0*norm;
		else
			hitNorm = norm;*/

		/*double t;
		if (!intersectMT(vertices[0]->pos, vertices[1]->pos, vertices[2]->pos, ray.origin, ray.dir, t))
			return false;

		intersect = ray.origin + t*ray.dir;*/

		//Interpolate normal if vertex normals are set
		if (vecLengthSquared(vertices[0]->norm) > 0 && vecLengthSquared(vertices[1]->norm) > 0 && vecLengthSquared(vertices[2]->norm) > 0)
		{			

			double a1 = glm::length(glm::cross(d2, d3)) * inv_area;
			double a2 = glm::length(glm::cross(d3, d1)) * inv_area;
			double a3 = glm::length(glm::cross(d1, d2)) * inv_area;


			normal = a1*vertices[0]->norm + a2*vertices[1]->norm + a3*vertices[2]->norm;

			/*if (glm::dot(hitNorm, normal) < 0)
				normal = -1.0*normal;*/

			//normal = tmpNorm;

			uv = a1*vertices[0]->texCoord + a2*vertices[1]->texCoord + a3*vertices[2]->texCoord;
		}
		else
			normal = norm;

		return true;

	}

	virtual bool intersect(const Ray& ray, double& t)
	{
		glm::dvec3 d1, d2, d3;

		double dot = glm::dot(ray.dir, norm);

		if (dot <= EPSILON && dot >= -EPSILON)
			return false;

		double t0 = glm::dot((vertices[0]->pos - ray.origin), norm) / dot;

		if (t < 0)
			return false;

		glm::dvec3 intersect = ray.origin + t0*ray.dir;

		d1 = intersect - vertices[0]->pos;
		if (glm::dot(glm::cross(edges[0], d1), norm) < 0) return false;

		d2 = (intersect - vertices[1]->pos);
		if (glm::dot(glm::cross(edges[1], d2), norm) < 0) return false;

		d3 = (intersect - vertices[2]->pos);
		if (glm::dot(glm::cross(edges[2], d3), norm) < 0) return false;

		t = t0;

		return true;
	}

	virtual bool intersect(BoundingBox bbox)
	{
		BoundingBox tmpBox = BoundingBox(bbox.min - EPSILON, bbox.max + EPSILON);
		glm::dvec3 verts[3] = { vertices[0]->pos, vertices[1]->pos , vertices[2]->pos };

		return triBoxOverlap(tmpBox.center(), glm::dvec3(tmpBox.dx() / 2, tmpBox.dy() / 2, tmpBox.dz() / 2), verts);
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
      /*if(max.x <= min.x)*/ max.x += EPSILON;
      /*if(max.y <= min.y)*/ max.y += EPSILON;
     /* if(max.z <= min.z)*/ max.z += EPSILON;

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
			/*double u = .5*acos((intersect.y - pos.y) / rad) / (M_PI)+.5;
			double v = .5*atan((intersect.z - pos.z) / (intersect.x - pos.x)) / (2 * M_PI) + .5;*/
			scene->push_back(new triangle(new vertex(rad*tris[i][0] + pos, tris[i][0], glm::dvec2(.5*acos((tris[i][0]).y) / (M_PI)+.5, .5*atan((tris[i][0]).z / (tris[i][0]).x) / (2 * M_PI) + .5)),
										  new vertex(rad*tris[i][1] + pos, tris[i][1], glm::dvec2(.5*acos((tris[i][1]).y) / (M_PI)+.5, .5*atan((tris[i][1]).z / (tris[i][1]).x) / (2 * M_PI) + .5)),
										  new vertex(rad*tris[i][2] + pos, tris[i][2], glm::dvec2(.5*acos((tris[i][2]).y) / (M_PI)+.5, .5*atan((tris[i][2] ).z / (tris[i][2]).x) / (2 * M_PI) + .5)), material));
			count++;
		}

		std::cout << "sphere tris: " << count << "\n";

	}

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
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

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		glm::dvec3 verts[8];
		glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY); //TODO: change this somehow
		glm::dmat3x3 irot = glm::inverse(rot);

		//vertices of bounding pyramid in object space
		verts[0] = rad*glm::dvec3(-1, -1, 0);
		verts[1] = rad*glm::dvec3(-1, 1, 0);
		verts[2] = rad*glm::dvec3(1, -1, 0);
		verts[3] = rad*glm::dvec3(1, 1, 0);
		verts[4] = glm::dvec3(0, 0, height);

		//transform into world space and compute axis aligned bounding box
		for (int i = 0; i < 5; i++)
		{
			verts[i] = verts[i] * irot  + pos;

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

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
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

	virtual bool intersect(const Ray& ray, glm::dvec3& intersect, glm::dvec3& normal, glm::dvec2& uv)
	{
		return false;
	}

	virtual BoundingBox boundingBox() const
	{
		return BoundingBox(glm::dvec3(0, 0, 0), glm::dvec3(1, 1, 1));
	}
};
