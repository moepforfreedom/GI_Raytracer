#include <QApplication>

#include <iostream>

#include "camera.h"
#include "gui.h"

int main(int argc, char** argv) {
    QApplication app(argc, argv);

	Camera camera({ 10, 0, 0 }, {0, 0, 0});
    glm::dvec3 light{10, 10, 10};

    RayTracer raytracer(camera, light);

    // Set up scene
    Octree scene({-20, -20, -20}, {20, 20, 20});
    // TODO Add objects to the scene
    // scene.push_back(...);

	scene.push_back(&sphere(glm::dvec3(0, 0, 0), 1, Material(glm::dvec3(1, 0, 0))));
	scene.push_back(&sphere(glm::dvec3(5, .5, 0), .5, Material(glm::dvec3(0, 1, 0))));
	scene.push_back(&sphere(glm::dvec3(6, .5, .75), .25, Material(glm::dvec3(0, 0, 1))));

    raytracer.setScene(&scene);

    Gui window(500, 500, raytracer);
    window.show();
    return app.exec();
}
