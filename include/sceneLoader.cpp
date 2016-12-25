#include <glm/glm.hpp>
#include <iostream>
#include "entities.h"
#include "octree.h"
#include "material.h"
#include "sceneLoader.h"
#include "meshLoader.h"
#include <cstring>

//imports a scene frm a file
void loadScene(Octree* o, const char* fname)
{
	//std::cout << "loading scene: " << fname << "\n";
	std::vector<texture*> tex;
	std::vector<Material*> mats;

	#ifndef _MSC_VER
		std::locale::global(std::locale("POSIX"));
	#endif

	std::string path = fname;
	std::size_t botDirPos = path.find_last_of("/");
	std::string dir = path.substr(0, botDirPos);

	FILE* f = fopen(fname, "r");

	if (f == NULL)
	{
		std::cout << "error while opening file: " << fname << "\n";
		return;
	}
	int faces = 0;

	std::cout << "loading scene: " << fname << "\n";

	while (1)
	{
		char lineHeader[128];
		// read the first word of the line
		int res = fscanf(f, "%s", lineHeader);

		//std::cout << lineHeader << "\n";

		if (res == EOF)
			break;

		if (std::strcmp(lineHeader, "imTex") == 0)
		{
			char fn[100];
			std::string filePath = dir;
			int utile, vtile;
			fscanf(f, "%s %d %d\n", fn, &utile, &vtile);
			filePath += fn;

			tex.push_back(new imageTexture(filePath.c_str(), glm::dvec2(utile, vtile)));
		}
		else if (std::strcmp(lineHeader, "checkerboardTex") == 0)
		{
			glm::dvec3 a, b;
			int tiles;
			fscanf(f, "%lf %lf %lf %lf %lf %lf %d\n", &a.x, &a.y, &a.z, &b.x, &b.y, &b.z, &tiles);

			tex.push_back(new checkerboard(tiles, a, b));
		}
		else if (std::strcmp(lineHeader, "colorTex") == 0)
		{
			glm::dvec3 col;
			fscanf(f, "%lf %lf %lf\n", &col.x, &col.y, &col.z);

			std::cout << "creating texture: " << col.x << ", " << col.y << ", " << col.z << "\n";

			tex.push_back(new texture(col));
		}
		else if (std::strcmp(lineHeader, "mat") == 0)
		{			
			int dif, em;
			double r, o;
			fscanf(f, "%d %d %lf %lf\n", &dif, &em, &r, &o);
			std::cout << "creating mat: " << dif << ", " << em << "\n";

			mats.push_back(new Material(tex[dif], tex[em], r, o));
		}
		else if (std::strcmp(lineHeader, "mesh") == 0)
		{
			char fn[100];
			std::string filePath = dir;
			glm::dvec3 pos;
			glm::dvec3 rot;
			int mat;
			fscanf(f, "%s %lf %lf %lf %lf %lf %lf %d\n", &fn, &pos.x, &pos.y, &pos.z, &rot.x, &rot.y, &rot.z, &mat);
			filePath += fn;

			std::cout << "loading mesh: " << fn << ", " << mat << "\n";

			loadOBJ(o, fn, pos, rot, *mats[mat]);
		}
		else if (std::strcmp(lineHeader, "sphere") == 0)
		{
			glm::dvec3 pos;
			double rad;
			int mat;
			fscanf(f, "%lf %lf %lf %lf %d\n", &pos.x, &pos.y, &pos.z, &rad, &mat);

			o->push_back(new sphere(pos, rad, *mats[mat]));
		}
		else if (std::strcmp(lineHeader, "box") == 0)
		{
			glm::dvec3 pos;
			glm::dvec3 size;
			glm::dvec3 rot;
			int mat;
			fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &pos.x, &pos.y, &pos.z, &size.x, &size.y, &size.z, &rot.x, &rot.y, &rot.z, &mat);

			boxMesh(o, pos, size, rot, *mats[mat]);
		}
		else if (std::strcmp(lineHeader, "light") == 0)
		{
			glm::dvec3 pos;
			glm::dvec3 col;
			double rad;
			fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &pos.x, &pos.y, &pos.z, &col.x, &col.y, &rad);

			o->push_back(new Light(pos, col, rad));
		}
	}

	std::cout << "read faces: " << faces << "\n";

	fclose(f);
}