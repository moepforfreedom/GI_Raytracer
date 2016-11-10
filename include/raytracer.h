#pragma once

#include <algorithm>
#include <memory>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"
#include "omp.h"

class RayTracer {
  public:
    RayTracer() = delete;
    RayTracer(const Camera& camera, glm::dvec3 light)
        : _camera(camera), _light(light), _image(std::make_shared<Image>(0,0)){};

    void setScene(const Octree* scene) { _scene = scene; }

    void run(int w, int h) {
        // TODO Implement this
        _image = std::make_shared<Image>(w, h);

        double sensorHalfWidth = (_camera.sensorDiag*w)/(sqrt((double)w*w + h*h));
        double sensorHalfHeight = sensorHalfWidth * ((double)h/w);

		glm::dvec3 screenCenter = _camera.pos + _camera.focalDist*_camera.forward;
		glm::dvec3 cameraRight = glm::cross(_camera.forward, _camera.up);

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for schedule(dynamic, 10)
        for (int y = 0; y < h; ++y) {
          if(_running)
          {
            for (int x = 0; x < w; ++x) {
                // TODO Implement this

				glm::dvec3 pixelPos = screenCenter + (sensorHalfWidth*((double)x/w - .5))*_camera.up + (sensorHalfHeight*((double)y / h - .5))*cameraRight;
				glm::dvec3 color(0, 0, 0);

				glm::dvec3 hit, minHit, minNorm;
				glm::dvec3 norm;

				bool intersected = false;

				Ray ray(_camera.pos, glm::normalize(pixelPos - _camera.pos));

                //simultate an expensive operation for performance testing
				/*for (int j = 0; j < 10000; j++)
				{
					int a = j;
					a = a*a;
				}*/

				std::vector<Entity*> objects = _scene->intersect(ray);

				std::vector<Entity*>::iterator it = objects.begin();

				Entity* current;

				while (it != objects.end())
				{		
					Entity* tmp = *it;
					if (tmp->intersect(ray, hit, norm))
					{
						if (it == objects.begin() || glm::length(hit - _camera.pos) < glm::length(minHit - _camera.pos))
						{
							current = tmp;
							minHit = hit;
							minNorm = norm;
							intersected = true;
						}
					}
					++it;
				}
				
				if (intersected)
				{
					double l = glm::dot(minNorm, glm::normalize(_light - _camera.pos));

					if (l < 0)
						l = 0;

					color = current->material.color * l;
				}
				else
					color = glm::dvec3(0, 0, 0);

                #pragma omp critical
                {
					_image->setPixel(x, y, color);// {1 * ((double)y / h), 1 * ((double)x / w), 0});
                }
            }
          }
        }
    }

    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    std::shared_ptr<Image> getImage() const { return _image; }

  private:
    bool _running = false;
    const Octree* _scene;
    Camera _camera;
    glm::dvec3 _light;
    std::shared_ptr<Image> _image;
};
