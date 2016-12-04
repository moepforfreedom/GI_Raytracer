#pragma once

#include <glm/glm.hpp>
#include <QImage>
#include <iostream>

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

//Procedural checkerboard texture
struct checkerboard: texture
{
	glm::dvec3 a, b;
    checkerboard(int t, glm::dvec3 col0, glm::dvec3 col1): texture(glm::dvec3(0, 0, 0)), tiles(t), a(col0), b(col1)
    {
    }

    virtual glm::dvec3 get(glm::dvec2 uv)
    {
        if(((int)(uv.x * tiles) % 2 == 0) ^ ((int)(uv.y * tiles) % 2 == 0))
            return a;

        return b;
    }

    int tiles;
};

//Texture from image file
struct imageTexture : texture
{
	const char* fname;
	QImage image;
	glm::dvec2 tile;

	imageTexture(const char* name, glm::dvec2 t) : texture(glm::dvec3(0, 0, 0)), fname(name), image(name), tile(t)
	{
		std::cout << "bits per pixel: " << image.pixelFormat().redSize() << "\n";
	}

	virtual glm::dvec3 get(glm::dvec2 uv)
	{
		auto p = image.pixelColor((int)(uv.x * image.width() * tile.x) % image.width(), (int)(uv.y * image.height() * tile.y) % image.height());

		return{ p.red() / 255.0, p.green() / 255.0, p.blue() / 255.0 };//*qRed(p) / 255., qGreen(p) / 255., qBlue(p) / 255. };
	}

	int tiles;
};

/// Represents the material properties of an entity. For now it only contains color, but it should
/// probably be extended to allow more options.
struct Material
{
        Material(texture* dif, texture* em, double r) : diffuse(dif), emissive(em), roughness(r) {}

    texture* diffuse;
	texture* emissive;
	double roughness;
};
