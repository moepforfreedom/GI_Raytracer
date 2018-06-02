## About

A simple raytracer written in C++

# Features
* Global Illumination via path tracing
* Reflection, refraction
* Area lights and soft shadows
* Adaptive quasi monte carlo sampling based on halton sequences
* Reflective and refractive caustics using photon mapping
* Several types of primitives including spheres, cones, triangles, boxes, etc.
* Simple materials with support for phong shading, diffuse/emissive maps
* Image-based and procedural textures
* Atmospheric effects such as fog and clouds
* Suport for importing OBJ meshes
* Octree-accelerated ray intersections
* Multithreading using OpenMP

## Requirements
* ``cmake``
* ``qt5``

## Building

~~~Bash
  cd <project-root>
  mkdir build; cd build
  cmake ..
  make
~~~

## Usage

~~~Bash
  ./global-illu <scene file>
~~~