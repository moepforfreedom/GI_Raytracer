#pragma once

#include <algorithm>
#include <memory>
#include <set>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"
#include "omp.h"

class RayTracer
{
  public:
    RayTracer() = delete;
    RayTracer(const Camera& camera, glm::dvec3 light)
        : _camera(camera), _light(light), _image(std::make_shared<Image>(0,0)){};

    void setScene(Octree* scene) { _scene = scene; }

    void run(int w, int h)
	{
        _image = std::make_shared<Image>(w, h);

		if (!_scene->valid)
		{
			_scene->rebuild();
		}

		glm::dmat3x3 crot = glm::eulerAngleXYZ(0.0, 0.02, 0.0);
		_camera.pos =  _camera.pos*crot;
		_camera.forward = glm::normalize(glm::dvec3(0, 0, 0) - _camera.pos);

		_camera.up = glm::dvec3(0, 1.0, 0);
		_camera.right = glm::normalize(glm::cross(_camera.up, _camera.forward));
		_camera.up = glm::cross(_camera.forward, _camera.right);

        double sensorHalfWidth = (_camera.sensorDiag*w)/(sqrt((double)w*w + h*h));
        double sensorHalfHeight = sensorHalfWidth * ((double)h/w);

		glm::dvec3 screenCenter = _camera.pos + _camera.focalDist*_camera.forward;
		glm::dvec3 cameraRight =  glm::normalize(glm::cross(_camera.forward, _camera.up));

		std::cout << cameraRight.x << ", " << cameraRight.y << ", " << cameraRight.z << ", " << "\n";


		double avgTests = 0;

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for schedule(dynamic, 10) //OpenMP
        for (int y = 0; y < h; ++y) 
		{
          if(_running)
          {
            for (int x = 0; x < w; ++x)
			{

				glm::dvec3 pixelPos = screenCenter + (sensorHalfWidth*((double)x/w - .5))*cameraRight - (sensorHalfHeight*((double)y / h - .5))*_camera.up;
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

				std::vector<Entity*> objects = _scene->intersect(ray, 0, INFINITY);

				avgTests += objects.size();

				std::vector<Entity*>::iterator it = objects.begin();

				Entity* current;

				while (it != objects.end())
				{
					Entity* tmp = *it;
					if (tmp->intersect(ray, hit, norm))
					{
						if (!intersected ||vecLengthSquared(hit - _camera.pos) < vecLengthSquared(minHit - _camera.pos))
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
					glm::dvec3 i(0, 0, 0);

					for (Light* light : _scene->lights)
					{
						bool shadow = false;						
						glm::dvec3 lightDir = light->getPoint() - minHit;
						double maxt = glm::length(lightDir);

						Ray shadow_ray(minHit + SHADOW_BIAS*minNorm, lightDir);


						std::vector<Entity*> shadow_objects = _scene->intersect(shadow_ray, 0, maxt);

						std::vector<Entity*>::iterator shadow_it = shadow_objects.begin();

						while (!shadow && shadow_it != shadow_objects.end())
						{
							Entity* t = *shadow_it;
							glm::dvec3 shadow_hit, shadow_norm;

							if (t->intersect(shadow_ray, shadow_hit, shadow_norm))
							{
								shadow = (vecLengthSquared(shadow_hit - minHit) < maxt*maxt);
							}
							shadow_it++;
						}

						if (!shadow)
						{
							double l = glm::dot(minNorm, glm::normalize(lightDir));

							if (l < 0)
								l = 0;

							i += light->col*l;
						}
					}					

					color = glm::clamp(current->material.color*i+ current->material.emissive, 0.0, 1.0);
				}
				else
					color = glm::dvec3(0, 0, 0);

                #pragma omp critical
                {
                  _image->setPixel(x, y, color);
                }
            }
          }
        }

		avgTests /= w*h;

		std::cout << "average intersection tests: " << avgTests << "\n";
    }

    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    std::shared_ptr<Image> getImage() const { return _image; }

  private:
    bool _running = false;
    Octree* _scene;
    Camera _camera;
    glm::dvec3 _light;
    std::shared_ptr<Image> _image;
};
