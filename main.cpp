#include <QApplication>

#include <iostream>

#include "camera.h"
#include "gui.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

	Camera camera({ 10, 0, 0 }, {0, 0, 0});
    glm::dvec3 light{10, 10, 10};

	camera.up = glm::dvec3(0, 1.0, 0);

    RayTracer raytracer(camera, light);

	

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
    // TODO Add objects to the scene
    // scene.push_back(...);

	scene.push_back(&sphere(glm::dvec3(0, 0, 0), 1, Material(glm::dvec3(1, 0, 0), glm::dvec3(0, 0, 0))));
	scene.push_back(&sphere(glm::dvec3(5, .5, 0), .5, Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 0, 0))));
	scene.push_back(&sphere(glm::dvec3(6, .5, .75), .25, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0))));
	glm::dvec3 pos = glm::dvec3(4.5, .5, .75);
	for (int i = 0; i < 8; i++)
	{		
		pos.x -= 4;
		pos.y -= 1;
		cone* a = new cone(pos, glm::dvec3(0.0 - .2*i, -0.0, 0), .25, 2, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0)));
		BoundingBox b = a->boundingBox();
		scene.push_back(new sphere(b.min, .05, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0))));
		scene.push_back(new sphere(b.max, .05, Material(glm::dvec3(0, 0, 1), glm::dvec3(.125, 0, 0))));
		scene.push_back(a);
	}
	scene.push_back(&triangle(&vertex(glm::dvec3(0, 0, 0), glm::dvec3(0, 0, 0)), &vertex(glm::dvec3(0, 3, 0), glm::dvec3(0, 0, 0)), &vertex(glm::dvec3(0, 0, 3), glm::dvec3(0, 0, 0)), Material(glm::dvec3(0, 1, 0), glm::dvec3(0, 0, 0))));
	scene.push_back(&cone(glm::dvec3(2.5, 0, .75), glm::dvec3(0, 0, 0), .35, 1, Material(glm::dvec3(1, 0, 0), glm::dvec3(.125, 0, 0))));
	scene.push_back(&sphere(glm::dvec3(0, 0, 0), 32, Material(glm::dvec3(0.25, 0.25, .25), glm::dvec3(.125, 0.125, .125))));

	scene.push_back(&coneMesh(&scene, glm::dvec3(0, 5, 0), glm::dvec3(0.0, -0.0, 0), 2.25, 5.0, 16, Material(glm::dvec3(1, 0, 1), glm::dvec3(0, 0, 0))));

    raytracer.setScene(&scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
}
