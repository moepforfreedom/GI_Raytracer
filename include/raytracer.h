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
#include "util.h"
#include "halton_enum.h"
#include "halton_sampler.h"

int p = 0;
int width, height;

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

        width = w;
        height = h;

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

        Halton_sampler sampler;
        sampler.init_faure();

        Halton_enum halton_enum(w, h);

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

                    int idx = halton_enum.get_index(s, x, y);

					if (s == 0)
						color = radiance(ray, 0, sampler, halton_enum, /*(x + w*y)*SAMPLES + s*/ idx, glm::dvec3(0, 0, 0));
					else
					color = (1.0*s*color + radiance(ray, 0, sampler, halton_enum, /*(x + w*y)*SAMPLES + s*/ idx, glm::dvec3(0, 0, 0)))*(1.0 / (s + 1));// (1.0 / SAMPLES)*radiance(ray, 0);

					if (s > 0)
					{
						var = (1.0*5*var + glm::length(color - lastCol))*(1.0 / (5 + 1));

						//var = .5*var + .5*vars[clamp(0, w, x - 1) + w*clamp(0, h, y)];
					}

					if (s > 0 && var > NOISE_THRESH)
						samps-=2;

					s++;
					samps++;
				}

				vars[x + w*y] = var;

                #pragma omp critical
                {
                  _image->setPixel(x, y, glm::clamp(color/*1.0*(double)s / SAMPLES/*SAMPLES*4*var*//**glm::dvec3(1, 1, 1)*/, 0.0, 1.0));
                }
            }
          }
        }

		avgTests /= w*h;

		std::cout << "average intersection tests: " << avgTests << "\n";
    }

	glm::dvec3 radiance(Ray ray, int depth, Halton_sampler& halton_sampler, Halton_enum& halton_enum, int sample, glm::dvec3 contrib)
	{
		if (depth > MAX_DEPTH)
			return glm::dvec3(0, 0, 0);


        float sx = halton_sampler.sample(0 + 2*depth, sample);
        float sy = halton_sampler.sample(1 + 2*depth, sample);

		double offset = SHADOW_BIAS;

        sx = fmod(halton_enum.scale_x(sx), 1.0);
        sy = fmod(halton_enum.scale_y(sy), 1.0);

        //std::cout << sx << ", " << sy << "\n";

		glm::dvec3 minHit, minNorm;
		glm::dvec2 minUV;

		Entity* current;
        bool backface = false;

		bool intersected = trace(ray, minHit, minNorm, minUV, current);

		if (intersected)
		{
			glm::dvec3 i(0, 0, 0);

			if (glm::dot(minNorm, ray.dir) > 0)
            {
				minNorm *= -1.0;
                backface = true;
            }

            /*if (vecLengthSquared(minNorm) > 1)
            {
				minNorm *= .5;
                backface = true;
            }*/

			for (Light* light : _scene->lights)
			{
				bool shadow = false;
				glm::dvec3 lightDir = light->getPoint() - minHit;
				double maxt = glm::length(lightDir);

				Ray shadow_ray(minHit + SHADOW_BIAS*minNorm, lightDir);

				shadow = !visible(shadow_ray, maxt);

				if (!shadow)
				{
					double d = glm::dot(minNorm, glm::normalize(lightDir));

					if (d < 0)
						d = 0;

					double l = pow(d, 1 / current->material.roughness);

					i += light->col*l;
				}
			}

			glm::dvec2 p = hammersley2d(rand() % 150, 150);

			double z = abs(minNorm.z);

			glm::dmat3x3 rot(z + (1.0 / (1 + z))*-minNorm.y*-minNorm.y, (1.0 / (1 + z))*(minNorm.x*-minNorm.y), -minNorm.x,
				(1.0 / (1 + z))*(minNorm.x*-minNorm.y), z + (1.0 / (1 + z))*-minNorm.x*-minNorm.x, -minNorm.y,
				minNorm.x, minNorm.y, z);


			glm::dvec2 g = importance_sample_ggx(p.x, p.y, .5);
			//glm::dvec3 tmpNorm = rot*hemisphereSample_cos(p.x, p.y, 2 / current->material.roughness);

			glm::dvec3 refDir;
			//refDir = glm::refract(ray.dir, minNorm, .75);

			if (minNorm.z < 0)
			{
				refDir.z = -1.0*refDir.z;
				//tmpNorm.z = -1.0*tmpNorm.z;
			}

            double IOR = current->material.IOR;

            double r0 = pow((1-IOR)/(1+IOR), 2);

            //Schlicks approximation of the fresnel term
            double fs = r0 + (1-r0)*pow(1-glm::dot(glm::reflect(ray.dir, minNorm), minNorm), 5);

            //type of secondary ray, 0 for reflection, 1 for refraction, 2 for glossy
            int type = 2;

			if(current->material.roughness < .001)
            {
                type = 0;
            }

            if (current->material.opacity < 1)
            {
                if(drand() < fs)
                    type = 0;
                else
                    type = 1;
            }


            if(type == 1)
			{
                if(backface)
                {
                    refDir = refr(ray.dir, minNorm, current->material.IOR);
                }
                else
                {
                    refDir = refr(ray.dir, minNorm, 1.0 / current->material.IOR);
                }
				offset *= -1;
			}
            else if (type == 0)
				refDir = glm::reflect(ray.dir, minNorm);
            else
                refDir = rot*hemisphereSample_cos(sx, sy, 2);



            if(vecLengthSquared(contrib) > 0.001)
            {
                i += 1.0/*glm::dot(ray.dir, refDir)*/*radiance(Ray(minHit + offset*minNorm, refDir), ++depth, halton_sampler, halton_enum, sample, contrib);
            }

			return current->material.diffuse->get(minUV)*i + current->material.emissive->get(minUV);
		}
		else
			return glm::dvec3(0, 0, 0);
	}


	//checks if the given ray is visible, meaning that nothing overlaps it
	bool visible(Ray& ray, double mt)
	{
		bool hit = false;
		std::vector<Entity*> shadow_objects = _scene->intersect(ray, SHADOW_BIAS, mt);

		std::vector<Entity*>::iterator shadow_it = shadow_objects.begin();

		while (!hit && shadow_it != shadow_objects.end())
		{
			Entity* t = *shadow_it;
			double t_shadow;

			if (t->intersect(ray, t_shadow))
			{
				hit = (t_shadow < mt);
			}
			++shadow_it;
		}

		return !hit;
	}

	bool trace(Ray& ray, glm::dvec3& minHit, glm::dvec3&minNorm, glm::dvec2& minUV, Entity*& obj)
	{
		glm::dvec3 hit, norm;
		glm::dvec2 uv;

		bool intersected = false;

		std::vector<const Octree::Node*> nodes = _scene->intersectSorted(ray, 0, INFINITY);

		std::vector<const Octree::Node*>::iterator nd = nodes.begin();

		//avgTests += objects.size();

		Entity* current;
		bool term = false;

		while (nd != nodes.end() && !term)
		{
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

						if ((vecLengthSquared(hit - _camera.pos) < pow(*curNode->maxt, 2)))
							term = true;
					}
				}
				++it;
			}
			++nd;
		}

		obj = current;
		return intersected;
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
