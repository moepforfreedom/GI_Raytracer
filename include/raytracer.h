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

		std::vector<double> xrand = subrand(2048);
		std::vector<double> yrand = subrand(2048);

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for schedule(dynamic, 10) //OpenMP
        for (int y = 0; y < h; ++y)
		{
          if(_running)
          {
            for (int x = 0; x < w; ++x)
			{
				glm::dvec3 color(0, 0, 0);

				for (int s = 0; s < SAMPLES; s++)
				{
					double xr = xrand[((x + w*y)*SAMPLES + s) % xrand.size()];
					double yr = yrand[((x + w*y)*SAMPLES + s) % yrand.size()];

					double dx = (double)x + AA_JITTER*xr;
					double dy = (double)y + AA_JITTER*yr;

					glm::dvec3 pixelPos = screenCenter + (sensorHalfWidth*(dx / w - .5))*cameraRight - (sensorHalfHeight*(dy / h - .5))*_camera.up;

					glm::dvec3 eyePos = _camera.pos + FOCAL_BLUR*xr*cameraRight + FOCAL_BLUR*yr *_camera.up;

					Ray ray(eyePos, glm::normalize(pixelPos - eyePos));

					color += (1.0/ SAMPLES)*radiance(ray, 0);

				}

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

	glm::dvec3 radiance(Ray& ray, int depth)
	{
		if (depth > MAX_DEPTH)
			return glm::dvec3(0, 0, 0);


		glm::dvec3 hit, minHit, norm, minNorm;
		glm::dvec2 uv, minUV;

		bool intersected = false;

		std::vector<const Octree::Node*> nodes = _scene->intersectSorted(ray, 0, INFINITY);

		std::vector<const Octree::Node*>::iterator nd = nodes.begin();

		//avgTests += objects.size();

		Entity* current;
		bool lastIntersect = false;

		while (nd != nodes.end() && !lastIntersect)
		{
			lastIntersect = intersected;

			const Octree::Node* curNode = *nd;

			std::vector<Entity*>::const_iterator it = curNode->_entities.begin();			

			while (it != curNode->_entities.end())
			{
				Entity* tmp = *it;
				if (tmp->intersect(ray, hit, norm, uv))
				{
					if ((!intersected || vecLengthSquared(hit - _camera.pos) < vecLengthSquared(minHit - _camera.pos)) && (vecLengthSquared(hit - _camera.pos) < pow(*curNode->maxt + .5, 2))	)
					{
						current = tmp;
						minHit = hit;
						minNorm = norm;
						minUV = uv;
						intersected = true;
					}
				}
				++it;
			}
			nd++;
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


				std::vector<Entity*> shadow_objects = _scene->intersect(shadow_ray, SHADOW_BIAS, maxt);

				std::vector<Entity*>::iterator shadow_it = shadow_objects.begin();

				while (!shadow && shadow_it != shadow_objects.end())
				{
					Entity* t = *shadow_it;
					glm::dvec3 shadow_hit, shadow_norm;
					glm::dvec2 shadow_uv;
					double t_shadow;

					if (t->intersect(shadow_ray, t_shadow))
					{
						shadow = (t_shadow < maxt);// (vecLengthSquared(shadow_hit - minHit) < maxt*maxt);
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

			return glm::clamp(current->material.diffuse->get(minUV)*i + current->material.emissive->get(minUV), 0.0, 1.0);// +radiance(Ray(minHit + SHADOW_BIAS*minNorm, glm::mix(glm::reflect(ray.dir, minNorm), randomUnitVec(), current->material.roughness)), ++depth), 0.0, 1.0);
		}
		else
			return glm::dvec3(1, 1, 0);
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
    std::shared_ptr<Image> _image;};


