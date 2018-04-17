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
	srand(std::time(0));

    QApplication app(argc, argv);

	Camera camera({ 10, 5, 0 }, {0, 0, 0});

    RayTracer raytracer(camera);


    // Set up scene, TODO: determine octree bounds from scene geometry
    Octree* scene = new Octree({-170, -170, -170}, {170, 170, 170});

	if(argc > 1)
		loadScene(scene, raytracer, argv[1]);
	else
		loadScene(scene, raytracer, "scenes/foliage/foliage.scn");

    raytracer.setScene(scene);

    Gui window(1000, 1000, raytracer);
    window.show();
    return app.exec();
}
