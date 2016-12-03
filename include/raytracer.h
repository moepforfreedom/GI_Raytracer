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
#include <ctime>

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

		srand(std::time(0));

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

		std::vector<double> vars;
		vars.reserve(w*h);

		std::vector<double> xrand;
		subrand(xrand, 25000);
		std::vector<double> yrand;
		subrand(yrand, 25000);

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for schedule(dynamic, 10) //OpenMP
        for (int y = 0; y < h; ++y)
		{
		  srand(std::time(0) + rand());
          if(_running)
          {
            for (int x = 0; x < w; ++x)
			{
				glm::dvec3 color(0.5, 0.5, 0.5);
				glm::dvec3 lastCol(0, 0, 0);
				double var = 0;
				int samps = 0;
				int s = 0;

				while (s < SAMPLES && samps < MIN_SAMPLES)
				{
					lastCol = color;

					double xr = 1*xrand[((x + w*y)*SAMPLES + s) % xrand.size()] + 0*drand();
					double yr = 1*yrand[((x + w*y)*SAMPLES + s) % yrand.size()] + 0*drand();

					//std::cout << (xr - yr) << "\n";

					double dx = (double)x + AA_JITTER*xrand[((x + w*y)*SAMPLES + s) % xrand.size() + 1];
					double dy = (double)y + AA_JITTER*yrand[((x + w*y)*SAMPLES + s) % yrand.size() + 1];

					glm::dvec3 pixelPos = screenCenter + (sensorHalfWidth*(dx / w - .5))*cameraRight - (sensorHalfHeight*(dy / h - .5))*_camera.up;

					glm::dvec3 eyePos = _camera.pos + FOCAL_BLUR*(xr-.5)*cameraRight + FOCAL_BLUR*(yr-.5) *_camera.up;

					Ray ray(eyePos, glm::normalize(pixelPos - eyePos));

					if (s == 0)
						color = radiance(ray, 0);
					else
					color = (1.0*s*color + radiance(ray, 0))*(1.0 / (s + 1));// (1.0 / SAMPLES)*radiance(ray, 0);

					if (s > 0)
					{
						var = (1.0*5*var + glm::length(color - lastCol))*(1.0 / (5 + 1));

						//var = .5*var + .5*vars[clamp(0, w, x - 1) + w*clamp(0, h, y)];
					}

					/*if (s > MIN_SAMPLES && var < NOISE_THRESH)
						s+=4;*/

					if (s > 0 && var > NOISE_THRESH)
						samps-=2;

					s++;
					samps++;
				}

				vars[x + w*y] = var;

                #pragma omp critical
                {
                  _image->setPixel(x, y, glm::clamp(color/*1.0*(double)s / SAMPLES/*SAMPLES*4*var**//**glm::dvec3(1, 1, 1)*/, 0.0, 1.0));
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
		bool term = false;

		while (nd != nodes.end() && !term)
		{
			//lastIntersect = intersected;

			const Octree::Node* curNode = *nd;

			std::vector<Entity*>::const_iterator it = curNode->_entities.begin();

			while (it != curNode->_entities.end())
			{
				Entity* tmp = *it;
				if (tmp->intersect(ray, hit, norm, uv))
				{
					if ((!intersected || vecLengthSquared(hit - _camera.pos) < vecLengthSquared(minHit - _camera.pos)))
					{
						current = tmp;
						minHit = hit;
						minNorm = norm;
						minUV = uv;
						intersected = true;

						if((vecLengthSquared(hit - _camera.pos) < pow(*curNode->maxt, 2)))
							term = true;
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

			glm::dvec2 p = hammersley2d(rand() % 250, 250);		

			double z = abs(minNorm.z);

			glm::dmat3x3 rot(z + (1.0 / (1 + z))*-minNorm.y*-minNorm.y, (1.0 / (1 + z))*(minNorm.x*-minNorm.y), -minNorm.x,
				(1.0 / (1 + z))*(minNorm.x*-minNorm.y), z + (1.0 / (1 + z))*-minNorm.x*-minNorm.x, -minNorm.y,
				minNorm.x, minNorm.y, z);

			glm::dvec3 refDir = rot*hemisphereSample_uniform(p.x, p.y);			

			if (minNorm.z < 0)
				refDir.z = -1.0*refDir.z;

			/*if (glm::dot(minNorm, refDir) < 0)
				minNorm = -1.0*minNorm;*/

			double w = PowerCosHemispherePdfW(minNorm, refDir, 1);

			refDir = glm::mix(glm::reflect(hit - _camera.pos, minNorm), refDir, 1.0);

			i += 2*glm::dot(minNorm, refDir)*radiance(Ray(minHit + SHADOW_BIAS*minNorm, refDir), ++depth);

			return glm::clamp(current->material.diffuse->get(minUV)*i + current->material.emissive->get(minUV), 0.0, 1.0);
		}
		else
			return glm::dvec3(1, 1, 1);
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
