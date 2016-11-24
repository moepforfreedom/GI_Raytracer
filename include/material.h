#pragma once

#include <glm/glm.hpp>

struct texture
{

    texture(glm::dvec3 col) : color(col)
    {

    }

    virtual glm::dvec3 get(glm::dvec2 uv)
    {
        return color;
    }

    glm::dvec3 color;
};

/// Represents the material properties of an entity. For now it only contains color, but it should
/// probably be extended to allow more options.
struct Material
{
        Material(texture* dif, texture* em) : diffuse(dif), emissive(em) {}

    texture* diffuse;
	texture* emissive;
};
