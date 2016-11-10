#pragma once

#include <algorithm>
#include <memory>
//#include <future>

#include <glm/glm.hpp>

#include "camera.h"
#include "entities.h"
#include "image.h"
#include "octree.h"
#include <unistd.h>

class RayTracer {
  public:
    RayTracer() = delete;
    RayTracer(const Camera& camera, glm::dvec3 light)
        : _camera(camera), _light(light), _image(std::make_shared<Image>(0,0)){};

    void setScene(const Octree* scene) { _scene = scene; }

    void run(int w, int h) {
        // TODO Implement this
        _image = std::make_shared<Image>(w, h);

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for
        for (int y = 0; y < h; ++y) {
          if(_running)
          {
            for (int x = 0; x < w; ++x) {
                // TODO Implement this
                #pragma omp critical
                {
                _image->setPixel(x, y, {1*((double)y/h), 1*((double)x/w), 0});
                }
                usleep(1);
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
