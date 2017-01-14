#include <QApplication>

#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>

#include "camera.h"
#include "gui.h"
#include "meshLoader.h"
#include "sceneLoader.h"
#include "util.h"
#include "atmosphere.h"
#include <ctime>
#include <iomanip>
#include <limits>

int main(int argc, char** argv)
{
	//std::srand(std::time(0));

	srand(std::time(0));

    QApplication app(argc, argv);

	Camera camera({ 10, 5, 0 }, {0, 0, 0});
    glm::dvec3 light{10, 10, 10};

    RayTracer raytracer(camera, light);


    // Set up scene
    Octree* scene = new Octree({-70, -70, -70}, {70, 70, 70});

	float test0[24], test1[24];
	float boxes[48];

	Ray r(glm::dvec3(0, 0, 0), glm::dvec3(0, 1, 0));


	intersectSIMD(test0, test1, boxes, r.r, r.invD, 0, INFINITY);


	/*scene->push_back(new sphere(glm::dvec3(0, 0, 0), 1, Material(new checkerboard(8, glm::dvec3(1, 0, 0), glm::dvec3(1, 1, 1)), new texture(glm::dvec3(0, 0, 0)), .1, 1)));
	scene->push_back(new sphere(glm::dvec3(2, -4.5, 0), .75, Material(new texture(glm::dvec3(1, 1, 1)), new texture(glm::dvec3(0, 0, 0)), .0001, .1)));
	scene->push_back(new sphere(glm::dvec3(6, .5, .75), .25, Material(new texture(glm::dvec3(0, 0, 1)), new texture(glm::dvec3(.125, 0, 0)), .75, 1)));
	glm::dvec3 pos = glm::dvec3(4.5, .5, .75);*/

	/*std::vector<glm::dvec3> points;
	points.reserve(25000);
	subrandUnitVec(points, 25000);*/

	/*Halton_sampler halton_sampler;

    halton_sampler.init_faure();
	const Halton_enum halton_enum(width, height);

	for (int i = 0; i < 5500; i++)
	{
		glm::dvec2 p = hammersley2d(rand() % 5500, 5500);

		glm::dvec3 norm = glm::normalize(glm::dvec3(0, 1, -1));

		//const unsigned index = halton_enum.get_index(i, x, y);

		const float sx = halton_sampler.sample(0, i);
        const float sy = halton_sampler.sample(1, i);

		double z = abs(norm.z);

		glm::dmat3x3 rot(z + (1.0 / (1 + z))*-norm.y*-norm.y, (1.0 / (1 + z))*(norm.x*-norm.y), -norm.x,
						 (1.0 / (1 + z))*(norm.x*-norm.y), z + (1.0 / (1 + z))*-norm.x*-norm.x, -norm.y,
						 norm.x, norm.y, z);

		glm::dvec3 dir = rot*hemisphereSample_cos(sx, sy, 1);

		scene->push_back(new sphere(2.0*dir,  .01, Material(new texture(glm::dvec3(0, 0, 1)), new texture(glm::dvec3(.125, 0, 0)), .75, 1)));
	}*/
	//scene->push_back(new triangle(new vertex(glm::dvec3(0, 0, 0)), new vertex(glm::dvec3(0, 3, 0)), new vertex(glm::dvec3(0, 0, 3)), Material(new texture(glm::dvec3(0, 1, 0)), new texture(glm::dvec3(0, 0, 0)), .75, 1)));
	//scene->push_back(new cone(glm::dvec3(2.5, 0, .75), glm::dvec3(0, 0, 0), .35, 1, Material(new texture(glm::dvec3(1, 0, 0)), new texture(glm::dvec3(.125, 0, 0)), .75, 1)));
	//scene->push_back(new sphere(glm::dvec3(0, 0, 0), 75, Material(new texture(glm::dvec3(0, 0, 0)), new imageTexture("HDR_029_Sky_Cloudy_Bg", glm::dvec2(1, 1)), 1, 1)));

	/*scene->push_back(new coneMesh(scene, glm::dvec3(3, 0, 0), glm::dvec3(0.5, -.5, 0.0), 1.5, 3.0, 32, Material(new texture(glm::dvec3(1, 0, 0)), new texture(glm::dvec3(0, 0, 0)), 1, 1)));

	sphereMesh(scene, glm::dvec3(0, 0, -5), 3, 3, Material(new texture(glm::dvec3(1, 0, 1)), new texture(glm::dvec3(0, 0, 0)), 1, 1));

	//scene->push_back(new sphere(glm::dvec3(0, 0, -5), 3, Material(new texture(glm::dvec3(1, 0, 1)), new texture(glm::dvec3(0, 0, 0)), 1, 1)));

	//scene->push_back(new quadMesh(scene, glm::dvec3(1, 0, 0), glm::dvec3(1, 1, 0), glm::dvec3(0, 0, 1), glm::dvec3(0, 1, .5), Material(glm::dvec3(0, 0, 1), glm::dvec3(0, 0, 0))));

	scene->push_back(new boxMesh(scene, glm::dvec3(-1, -1, -1), glm::dvec3(2, 2, 2), glm::dvec3(.5, 1, 0), Material(new texture(glm::dvec3(1, 0, 0)), new texture(glm::dvec3(0, 0, 0)), .75, 1)));

	loadOBJ(scene, "teapot.obj", glm::dvec3(-6, 0, 4), glm::dvec3(0, 0, 1), Material(new texture(glm::dvec3(0.25, 1, 1)), new texture(glm::dvec3(0, 0, 0)), .0002, .1));

	*///loadOBJ(scene, "terrain.obj", glm::dvec3(-0, -9, -10), glm::dvec3(0, 0, 0), Material(new imageTexture("Stone_01_Diffuse", glm::dvec2(16, 16)), new texture(glm::dvec3(0, 0, 0)), 1, 1));

	loadScene(scene, "scenes/cornell/test.scn");

	//scene->push_back(new HeightFog(glm::dvec3(0, .5, 0), glm::dvec3(10, 1, 10), glm::dvec3(1, 1, 1), 220, 0.5, 4));

	//loadOBJ(scene, "test.obj", glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0), Material(new texture(glm::dvec3(1, 1, 1)), new texture(glm::dvec3(0, 0, 0)), 1, 1));
	
	scene->push_back(new Light(glm::dvec3(0, 5, 0), 4.0*glm::dvec3(1, 1, 1), .5));


    raytracer.setScene(scene);

    Gui window(1000, 1000, raytracer);
    window.show();
    return app.exec();
}
