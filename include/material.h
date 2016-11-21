#pragma once

#include <glm/glm.hpp>

/// Represents the material properties of an entity. For now it only contains color, but it should
/// probably be extended to allow more options.
struct Material 
{
    constexpr explicit Material(glm::dvec3 color, glm::dvec3 emissive) : color(std::move(color)), emissive(std::move(emissive)) {}

    glm::dvec3 color;
	glm::dvec3 emissive;
};