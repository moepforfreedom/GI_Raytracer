#include <glm/glm.hpp>
#include <iostream>
#include "entities.h"
#include "octree.h"
#include "material.h"
#include "meshLoader.h"

//imports a mesh from an OBJ file and adds it to the scene, based on http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/
void loadOBJ(Octree* scene, const char* fname, glm::dvec3 pos, glm::dvec3 rotation, const Material& material)
{
	std::vector< glm::vec3 > verts;
	std::vector< glm::vec2 > uvs;
	std::vector< glm::vec3 > normals;

	glm::dmat3x3 rot = glm::eulerAngleXYZ(rotation.x, rotation.y, rotation.z);

	FILE* f;

	fopen_s(&f, fname, "r");

	if (f == NULL)
	{
		std::cout << "error while opening file: " << fname << "\n";
		return;
	}

	while (1)
	{
		char lineHeader[128];
		// read the first word of the line
		int res = fscanf_s(f, "%s", lineHeader, 128);

		if (res == EOF)
			break;

		if (strcmp(lineHeader, "v") == 0)
		{
			glm::vec3 vertex;
			fscanf(f, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
			verts.push_back(rot*vertex + pos);
		}
		else if (strcmp(lineHeader, "vt") == 0)
		{
			glm::vec2 uv;
			fscanf(f, "%f %f\n", &uv.x, &uv.y);
			uvs.push_back(uv);
		}
		else if (strcmp(lineHeader, "vn") == 0)
		{
			glm::vec3 normal;
			fscanf(f, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
			normals.push_back(rot*normal);
		}
		else if (strcmp(lineHeader, "f") == 0)
		{

			unsigned int vertIndex[3], uvIndex[3], normIndex[3];
			int matches = fscanf(f, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertIndex[0], &uvIndex[0], &normIndex[0], &vertIndex[1], &uvIndex[1], &normIndex[1], &vertIndex[2], &uvIndex[2], &normIndex[2]);
			
			if (matches != 9) 
			{
				std::cout << "error while reading faces \n";
				return;
			}

			scene->push_back(new triangle(new vertex(verts[vertIndex[0] - 1], normals[normIndex[0] - 1], uvs[uvIndex[0] - 1]),
										  new vertex(verts[vertIndex[1] - 1], normals[normIndex[1] - 1], uvs[uvIndex[1] - 1]),
										  new vertex(verts[vertIndex[2] - 1], normals[normIndex[2] - 1], uvs[uvIndex[2] - 1]), material));
		}
	}
}

