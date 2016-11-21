#include <QApplication>

#include <iostream>

#include "camera.h"
#include "gui.h"
#include "meshLoader.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

	Camera camera({ 20, 10, 0 }, {0, 0, 0});
    glm::dvec3 light{10, 10, 10};

    RayTracer raytracer(camera, light);

    // Set up scene
    Octree* scene = new Octree({-30, -30, -30}, {30, 30, 30});
    // TODO Add objects to the scene
    // scene.push_back(...);


	scene->push_back(new sphere(glm::dvec3(0, 0, 0), 1, Material(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0))));
	scene->push_back(new sphere(glm::dvec3(5, .5, 0), .5, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 0, 0))));
	scene->push_back(new sphere(glm::dvec3(6, .5, .75), .25, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0))));
	glm::dvec3 pos = glm::dvec3(4.5, .5, .75);
	for (int i = 0; i < 16; i++)
	{
		pos.x -= 4;
		pos.y -= 1;
		cone* a = new cone(pos, glm::dvec3(0.0 - .2*i, -0.0, 0), .25, 2, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0)));

		scene->push_back(a);
	}
	scene->push_back(new triangle(new vertex(glm::dvec3(0, 0, 0)), new vertex(glm::dvec3(0, 3, 0)), new vertex(glm::dvec3(0, 0, 3)), Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 0, 0))));
	scene->push_back(new cone(glm::dvec3(2.5, 0, .75), glm::dvec3(0, 0, 0), .35, 1, Material(glm::dvec3(1, 0, 0), glm::dvec3(.125, 0, 0))));
	scene->push_back(new sphere(glm::dvec3(0, 0, 0), 128, Material(glm::dvec3(0.25, 0.25, .25), glm::dvec3(.125, 0.125, .125))));

	scene->push_back(new coneMesh(scene, glm::dvec3(3, 0, 0), glm::dvec3(0.5, -.5, 0.0), 1.5, 3.0, 32, Material(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0))));

	scene->push_back(new sphereMesh(scene, glm::dvec3(0, 0, -5), 3, 3, Material(glm::dvec3(1, 0, 1), glm::dvec3(0, 0, 0))));

	scene->push_back(new quadMesh(scene, glm::dvec3(1, 0, 0), glm::dvec3(1, 1, 0), glm::dvec3(0, 0, 1), glm::dvec3(0, 1, .5), Material(glm::dvec3(0, 0, 1), glm::dvec3(0, 0, 0))));

	scene->push_back(new boxMesh(scene, glm::dvec3(-1, -1, -1), glm::dvec3(2, 2, 2), glm::dvec3(.5, 1, 0), Material(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0))));

	loadOBJ(scene, "teapot.obj", glm::dvec3(-6, 0, 4), glm::dvec3(0, 0, 1), Material(glm::dvec3(0.25, 1, 1), glm::dvec3(0, 0, 0)));

	loadOBJ(scene, "terrain.obj", glm::dvec3(-0, -6, -10), glm::dvec3(0, 0, 0), Material(glm::dvec3(1, .65, .25), glm::dvec3(0, 0, 0)));

    raytracer.setScene(scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
}
