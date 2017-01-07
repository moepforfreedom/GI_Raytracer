#pragma once

#define GLM_FORCE_SSE2 // or GLM_FORCE_SSE42 if your processor supports it
#define GLM_FORCE_ALIGNED

#include <glm/glm.hpp>
#include <QImage>
#include <iostream>

struct texture
{
    texture(glm::dvec3 col) : color(col)
    {

    }

    virtual glm::dvec3 get(glm::dvec2& uv)
    {
        return color;
    }

    virtual double getAlpha(glm::dvec2& uv)
    {
        return 1;
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

    virtual glm::dvec3 get(glm::dvec2& uv)
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
		std::cout << "loading texture: " << fname << "\n";
		std::cout << "bits per pixel: " << image.pixelFormat().redSize() << "\n";
	}

	virtual glm::dvec3 get(glm::dvec2& uv)
	{
		auto p = image.pixelColor((int)(uv.x * image.width() * tile.x) % image.width(), (int)(uv.y * image.height() * tile.y) % image.height());

		return{ p.red() / 255.0, p.green() / 255.0, p.blue() / 255.0 };//*qRed(p) / 255., qGreen(p) / 255., qBlue(p) / 255. };
	}

    virtual double getAlpha(glm::dvec2& uv)
	{
        if(!image.hasAlphaChannel())
            return 1;

		auto p = image.pixelColor((int)(uv.x * image.width() * tile.x) % image.width(), (int)(uv.y * image.height() * tile.y) % image.height());

		return p.alpha() / 255.0;
	}

	int tiles;
};

/// Represents the material properties of an entity.
struct Material
{
        Material(texture* dif, texture* em, double r, double o) : diffuse(dif), emissive(em), roughness(r), opacity(o) 
		{
			std::cout << "roughness: " << roughness << "\n";
		}

    texture* diffuse;
	texture* emissive;
	double roughness;
	double opacity;
	double IOR = 1.5;
};
