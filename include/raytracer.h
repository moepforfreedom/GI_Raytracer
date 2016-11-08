#pragma once

#include <algorithm>
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
        : _camera(camera), _light(light), _image(0, 0){};

    void setScene(const Octree* scene) { _scene = scene; }

    void run(int w, int h) {
        // TODO Implement this
        _image = Image(w, h);

        // The structure of the for loop should remain for incremental rendering.
        #pragma omp parallel for
        for (int y = 0; y < h; ++y) {
            if(_running)
            {
            for (int x = 0; x < w; ++x) {
                // TODO Implement this
                _image.setPixel(x, y, {1, 0, 0});
                usleep(5);
            }
          }
        }
    }

    bool running() const { return _running; }
    void stop() { _running = false; }
    void start() { _running = true; }

    const Image& getImage() const { return _image; }

  private:
    bool _running = false;
    const Octree* _scene;
    Camera _camera;
    glm::dvec3 _light;
    Image _image;
};
