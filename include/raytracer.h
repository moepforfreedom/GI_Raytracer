#pragma once

#include <algorithm>
#include <memory>
#include <set>

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
#include "photonMap.h"

int p = 0;
int width, height;

class RayTracer
{
  public:
    RayTracer() = delete;
	RayTracer(const Camera& camera, glm::dvec3 light)
		: _camera(camera), _light(light), _image(std::make_shared<Image>(0, 0))
	{

	};

    void setScene(Octree* scene)
	{
		_scene = scene;
		_photon_map = new PhotonMap(_scene->_root._bbox.min, _scene->_root._bbox.max);
	}

    void run(int w, int h)
	{
		srand(std::time(0));

		Halton_sampler sampler;
		sampler.init_faure();

		Halton_enum halton_enum(w, h);

        _image = std::make_shared<Image>(w, h);

        width = w;
        height = h;

		if (!_scene->valid)
		{
			_scene->rebuild();
		}

		if (!_photon_map->valid)
		{
            long start = clock();
			std::cout << "emitting photons...\n";
			tracePhotons(5, 750000, sampler, halton_enum);

            std::cout << "photon time: " << (clock() - start)/CLOCKS_PER_SEC << "\n";

			_photon_map->rebuild();
			//_scene->rebuild();
		}

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

					double dx = (double)x + AA_JITTER*xrand[((x + w*y)*SAMPLES + s) % xrand.size()];
					double dy = (double)y + AA_JITTER*yrand[((x + w*y)*SAMPLES + s) % yrand.size()];

					glm::dvec3 pixelPos = screenCenter + (sensorHalfWidth*(dx / w - .5))*cameraRight - (sensorHalfHeight*(dy / h - .5))*_camera.up;

					glm::dvec3 eyePos = _camera.pos +FOCAL_BLUR*(xr - .5)*cameraRight + FOCAL_BLUR*(yr - .5) *_camera.up;

					Ray ray(eyePos, glm::normalize(pixelPos - eyePos));

                    int idx = halton_enum.get_index(s, x, y);

					if (s == 0)
						color = radiance(ray, 0, sampler, halton_enum, /*(x + w*y)*SAMPLES + s*/ idx, glm::dvec3(1, 1, 1));
					else
					color = (1.0*s*color + radiance(ray, 0, sampler, halton_enum, /*(x + w*y)*SAMPLES + s*/ idx, glm::dvec3(1, 1, 1)))*(1.0 / (s + 1));// (1.0 / SAMPLES)*radiance(ray, 0);

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

				color = gamma(color, 2.2);

				//vars[x + w*y] = var;

                #pragma omp critical (im_update)
                {
                  _image->setPixel(x, y, glm::clamp(color/*1.0*(double)s / SAMPLES/*SAMPLES*4*var s**glm::dvec3(1, 1, 1)*/, 0.0, 1.0));
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


		glm::dvec3 minHit, minNorm;
		glm::dvec2 minUV;

		Entity* current;
        bool backface = false;

		bool intersected = trace(ray, minHit, minNorm, minUV, current);

		if (intersected)
		{
			glm::dvec3 i(0, 0, 0);

			glm::dvec3 refDir;
			glm::dvec3 color = current->material.diffuse->get(minUV);
			double roughness = current->material.roughness;

            //type of secondary ray, 0 for reflection, 1 for refraction, 2 for glossy
            // int type = rayType(current, ray, minNorm, minUV);
			glm::dvec3 f(1, 1, 1);

			secondaryRay(ray, current, minNorm, minUV, sx, sy, refDir, f, roughness, contrib, offset);

			double tmin = 0;
			double tmax = glm::length(minHit - ray.origin);

			if (_scene->atmosphereBounds(ray, tmin, tmax))
			{
				glm::dvec3 hit, col;

				//std::cout << "atmosphere bounds hit: " << tmin << ", " << tmax << "\n";

				if (raymarch(ray, hit, col, tmin, tmax))
				{
					minHit = hit;
					refDir = randomUnitVec(sx, sy);
					f = 1.0*col;
					color = col;
					contrib = col;
					roughness = 1;
					//std::cout << "atmosphere hit\n";
				}
			}

			for (Light* light : _scene->lights)
			{
				bool shadow = false;
				glm::dvec3 lightDir = light->getPoint() - (minHit+SHADOW_BIAS*minNorm);
				double maxt = glm::length(lightDir);

				//double cos_alpha = (light->rad / sqrt(vecLengthSquared(lightDir) + std::pow(light->rad, 2)));

				double hfrac = 1 / (M_PI*vecLengthSquared(light->pos - minHit)); //fraction of the hemisphere


				Ray shadow_ray(minHit + SHADOW_BIAS*minNorm, lightDir);

				shadow = !visible(shadow_ray, maxt);

				if (!shadow)
				{
					double d = glm::dot(minNorm, glm::normalize(light->pos - minHit));

					if (d < 0)
						d = 0;

					double l = pow(d, (1.0 / roughness));

					i = light->col*l*hfrac;
				}
			}

			glm::dvec3 caustic =  depth == 0 ? samplePhotons(minHit, minNorm, 42) : glm::dvec3(0, 0, 0);

			//return color*caustic;

			// continuation probability
			double q = compMax(contrib);

			if (depth <= MIN_DEPTH || drand() < q)
			{
				f *= depth <= MIN_DEPTH ? 1.0 : (1.0 / q);
				//diffuse*direct_light + diffuse*brdf*radiance + emmissive
				return color*i + f*radiance(Ray(minHit + offset*minNorm, refDir), ++depth, halton_sampler, halton_enum, sample, contrib) + current->material.emissive->get(minUV) + color*caustic;
			}
			else
				return glm::dvec3(0, 0, 0);///*1.0*depth / MAX_DEPTH * glm::dvec3(1, 1, 1);//*/ //current->material.diffuse->get(minUV)*i + current->material.emissive->get(minUV);
		}
		else
			return /*1.0*depth / MAX_DEPTH * glm::dvec3(1, 1, 1);//*/ glm::dvec3(0.04, 0.04, 0.04);
	}


	//checks if the given ray is visible, meaning that nothing overlaps it
	bool visible(Ray& ray, double mt)
	{
		bool hit = false;
		std::vector<Entity*> shadow_objects = _scene->intersect(ray, 0, mt-SHADOW_BIAS);

		std::vector<Entity*>::iterator shadow_it = shadow_objects.begin();

		while (!hit && shadow_it != shadow_objects.end())
		{
			Entity* t = *shadow_it;
			double t_shadow;

			if (t->intersect(ray, t_shadow))
			{
				hit = (t_shadow < mt)&&(t_shadow > 0);
			}
			++shadow_it;
		}

		if (hit)
			return false;

		double tmin = 0;
		double tmax = mt;

		if (_scene->atmosphereBounds(ray, tmin, tmax))
		{
			glm::dvec3 hit, col;

			if (raymarch(ray, hit, col, tmin, tmax)) return false;
		}

		return true;
	}

	void secondaryRay(Ray& ray, Entity* current, glm::dvec3& norm, glm::dvec2& UV, double sx, double sy, glm::dvec3& refDir, glm::dvec3& f, double& roughness, glm::dvec3& contrib, double& offset)
	{
		bool backface = false;

		if (glm::dot(norm, ray.dir) > 0)
		{
			norm *= -1.0;
			backface = true;
		}

		glm::dvec3 color = current->material.diffuse->get(UV);
		roughness = current->material.roughness;

		int type = rayType(current, ray, norm, UV);

		if (type == 1)
		{
			if (backface)
			{
				refDir = refr(ray.dir, norm, current->material.IOR);
			}
			else
			{
				refDir = refr(ray.dir, norm, 1.0 / current->material.IOR);
			}
			offset *= -1;

			contrib = glm::dvec3(1, 1, 1);

			f = 1.0*color;
		}
		else if (type == 0)
		{
			refDir = glm::reflect(ray.dir, norm);
			contrib = glm::dvec3(1, 1, 1);

			f = 1.0*color; // glm::dot(refDir, minNorm);
		}
		else
		{
			refDir = hemisphereSample_cos(norm, sx, sy, 2);

			if (current->material.roughness < .9)
			{
				refDir = sample_phong(glm::reflect(ray.dir, norm), norm, (1.0 / (current->material.roughness)) + 1, sx, sy);

				if (glm::dot(refDir, norm) < 0)
					refDir = glm::reflect(refDir, norm);

			}

			f = 1.0*color;

			glm::dvec3 inf = color;// *pow(dot, 1 / current->material.roughness);

			contrib *= inf;
			contrib = glm::mix(contrib, inf, 0.5);
		}
	}

	//traces a ray against the scene geometry, returns true on intersection
	bool trace(Ray& ray, glm::dvec3& minHit, glm::dvec3& minNorm, glm::dvec2& minUV, Entity*& obj)
	{
		glm::dvec3 hit, norm;
		glm::dvec2 uv;

		bool intersected = false;

		std::vector<const Octree::Node*> nodes = _scene->intersectSorted(ray, 0, INFINITY);

		std::vector<const Octree::Node*>::iterator nd = nodes.begin();

		//avgTests += objects.size();

		if (drand() < .5 && nodes.size() > 500)
		{
			std::cout << nodes.size() << ", ray: " << ray.origin.x << ", " << ray.origin.y << ", " << ray.origin.z << "\n";
		}

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
					if ((!intersected || vecLengthSquared(hit - ray.origin) < vecLengthSquared(minHit - ray.origin)))
					{
						current = tmp;
						minHit = hit;
						minNorm = norm;
						minUV = uv;
						intersected = true;

						if ((vecLengthSquared(hit - ray.origin) < std::pow(*curNode->maxt, 2)))
							term = true;
					}
				}
				++it;
			}
			++nd;
		}

		if(intersected)
			obj = current;

		return intersected;
	}

	//returns the type of the secondary ray, 0 for reflection, 1 for refraction, 2 for diffuse/glossy
	int rayType(Entity* entity, Ray& ray, glm::dvec3& norm, glm::dvec2& minUV)
	{
		int type = 2;

		double IOR = entity->material.IOR;

		double opacity = entity->material.diffuse->getAlpha(minUV) * entity->material.opacity;

		double r0 = std::pow((1 - IOR) / (1 + IOR), 2);

		//Schlicks approximation of the fresnel term
		double fs = r0 + (1 - r0)*std::pow(1 - glm::dot(glm::reflect(ray.dir, norm), norm), 5);

		if (entity->material.roughness < .001)
		{
			type = 0;
		}

		if (drand() > opacity)
		{
			if (drand() < fs)
				type = 0;
			else
				type = 1;
		}

		return type;
	}

	//uses raymarching to determine the intersection point of a ray with the atmosphere
	bool raymarch(Ray& r, glm::dvec3& hit, glm::dvec3& col, double mint, double maxt)
	{
		double t = mint + SHADOW_BIAS;
		double scatter;

		glm::dvec3 current = r.origin + mint*r.dir;

		while (t < maxt)
		{
			if (drand() < _scene->atmosphereDensity(current, col, scatter))
			{
				hit = current;
				return true;
			}

			current += RAYMARCH_STEPSIZE * r.dir;
			t += RAYMARCH_STEPSIZE;
		}

		return false;
	}

	//computes a radiance estimate from the photons surrounding the hit point
	glm::dvec3 samplePhotons(glm::dvec3 pos, glm::dvec3 dir, int count)
	{
		double dist = .0;
		glm::dvec3 res(0, 0, 0);
		double scale = 0;

		std::vector<Photon*> photons = _photon_map->getInRange(pos, scale, dist);

        /*if(photons.size() > 9*MAX_PHOTONS_PER_LEAF && drand() < .01)
        {
            std::cout << "photons: " << photons.size() << "\n";
        }*/

        count = std::min(count, (int)photons.size());

		std::partial_sort(photons.begin(), photons.begin() + count, photons.end(), [pos](const Photon* lhs, const Photon* rhs) {return vecLengthSquared(lhs->origin - pos) < vecLengthSquared(rhs->origin - pos); });

		double maxDist = 0;

		//std::cout << photons.size() << "\n";

        double lastDist = 0;

        if(photons.size() > 0)
            lastDist = vecLengthSquared(photons[0]->origin - pos);

		for (int i = 0; i < std::max(count - 1, 0); i++)
		{
			Photon* p = photons[i];

            /*if((vecLengthSquared(p->origin - pos) - lastDist) < 0)
            {
                std::cout << "non-ascending photon distance: " << lastDist << ", " << vecLengthSquared(p->origin - pos) << "\n";
            }

            lastDist = vecLengthSquared(p->origin - pos);*/

			res += .1*p->col*glm::dot(p->dir, dir);
		}

		if(photons.size() > 1)
        {
            maxDist = vecLengthSquared(photons[count - 2]->origin - pos);
			res /= (M_PI*maxDist);
        }

		return res;
	}

	//traces caustic photons from every light source
	void tracePhotons(int maxDepth, int count, Halton_sampler& halton_sampler, Halton_enum& halton_enum)
	{
		#pragma omp parallel
		{
			std::vector<Photon*> tmp;

			#pragma omp for
			for (int i = 0; i < count; i++)
			{
				for (Light* l : _scene->lights)
				{
					int tries = 0;
					bool stored = false;

					srand(i);

					while (!stored && tries < 500)
					{
						float sx = halton_sampler.sample(0, i * 500 + tries);
						float sy = halton_sampler.sample(1, i * 500 + tries);

						//sx = fmod(halton_enum.scale_x(sx), 1.0);
						//sy = fmod(halton_enum.scale_y(sy), 1.0);

						//std::cout << i * 500 + tries << ": " << sx << ", " << sy << "\n";

						glm::dvec3 pos = l->getPointInRange(1, sx, sy);
						glm::dvec3 dir = hemisphereSample_cos(glm::normalize(pos - l->pos), fmod(drand() + 5 * i, 1), fmod(drand() + 13 * i, 1), 2);

						Ray r(pos + SHADOW_BIAS*glm::normalize(pos - l->pos), dir);

						glm::dvec3 hit, norm;
						glm::dvec3 col = (1.0/count)*l->col;
						glm::dvec2 UV;
						Entity* current;

						int depth = 0;
						bool term = false;
						bool isCaustic = false;

						if (!trace(r, hit, norm, UV, current))
						{
							tries++;
							continue;
						}


						while (depth < maxDepth && !term)
						{
							double roughness = current->material.roughness;
							//std::cout << "tracing photon, depth: " << depth << "\n";
							if (roughness < 0.1)
							{
								//std::cout << "tracing caustics photon at depth: " << depth << "\n";
								if (!trace(r, hit, norm, UV, current))
								{
									term = true;
									continue;
								}

								roughness = current->material.roughness;
								glm::dvec3 refDir, f, contrib;
								double offset = SHADOW_BIAS;

								sx = halton_sampler.sample(0, maxDepth*(i * 500 + tries) + depth);
								sy = halton_sampler.sample(1, maxDepth*(i * 500 + tries) + depth);

								//sx = fmod(halton_enum.scale_x(sx), 1.0);
								//sy = fmod(halton_enum.scale_y(sy), 1.0);

								secondaryRay(r, current, norm, UV, fmod(drand() + 5 * i, 1), fmod(drand() + 13 * i, 1), refDir, f, roughness, contrib, offset);

								double tmin = 0;
								double tmax = glm::length(hit - r.origin);

								if (_scene->atmosphereBounds(r, tmin, tmax))
								{
									glm::dvec3 ahit, color;

									//std::cout << "atmosphere bounds hit: " << tmin << ", " << tmax << "\n";

									if (raymarch(r, ahit, color, tmin, tmax))
									{
										hit = ahit;
										refDir = randomUnitVec(fmod(rand() + 13 * i, 1), fmod(drand() + 7 * i, 1));
										f = 1.0*color;
										roughness = 1;
										//std::cout << "atmosphere hit\n";
									}
								}

								col *= f;

								r.origin = hit + offset*norm;
								r.setDir(refDir);

								isCaustic = true;
							}

							if (depth > 0 && isCaustic && roughness >= 0.1)
							{
								//std::cout << "photon stored\n";

								tmp.push_back(new Photon(hit, r.dir, col));
								term = true;
								stored = true;

								//_scene->push_back(new sphere(hit, 0.01, Material(new texture(col), new texture(glm::dvec3(0, 0, 0)), 1, 1)));

							}

							depth++;
						}
						tries++;
					}
				}
			}

			#pragma omp critical
			for (Photon* p : tmp)
			{
				_photon_map->push_back(p);
			}
		}
	}


    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    std::shared_ptr<Image> getImage() const { return _image; }

  private:
    bool _running = false;
    Octree* _scene;
	PhotonMap* _photon_map;
    Camera _camera;
    glm::dvec3 _light;
    std::shared_ptr<Image> _image;};
