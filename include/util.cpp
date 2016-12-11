#include <glm/glm.hpp>
#include <random>
#include <iostream>
#include <ctime>
#include "util.h"


glm::dvec3 randomVec()
{
	return glm::dvec3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
}

glm::dvec3 randomUnitVec()
{
	return glm::normalize(glm::dvec3(2 * drand() - 1, 2 * drand() - 1, 2 * drand() - 1));
}

// generates a point of the hammersley point set, based on http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
glm::dvec2 hammersley2d(unsigned int i, unsigned int N)
{
	return glm::dvec2(float(i) / float(N), radicalInverse_VdC(i));
}

glm::dvec3 hemisphereSample_uniform(float u, float v)
{
	float phi = v * 2.0f * M_PI;
	float cosTheta = 1.0f - u;
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

glm::dvec3 hemisphereSample_cos(float u, float v, double power)
{
	float phi = v * 2.0f * M_PI;
	float cosTheta = pow(1.0f - u, (1.0f / power));
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

double PowerCosHemispherePdfW(glm::dvec3  aNormal, glm::dvec3  aDirection, double  aPower)
{
	float cosTheta = std::max(0.0, glm::dot(aNormal, aDirection));
	return (aPower + 1.0) * std::pow(cosTheta, aPower) * ((1.0 / M_PI) * 0.5);
}

void subrand(std::vector<double>&out, int n)
{
	double primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

	double lastRand = drand();

	double a = fmod(sqrt(primes[rand() % 11]), 1.0);

	out.clear();

	for (int i = 0; i < n; i++)
	{
		lastRand = fmod(lastRand + a, 1.0);

		out.push_back(lastRand);
	}
}

//generates a sequence of subrandom unit vectors
void subrandUnitVec(std::vector<glm::dvec3>&out, int n)
{
	double primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

	glm::dvec3 lastRand = randomVec();

	glm::dvec3 tmp;

	glm::dvec3 a = glm::dvec3(fmod(sqrt(primes[rand() % 11]), 1.0), fmod(sqrt(primes[rand() % 11]), 1.0), fmod(sqrt(primes[rand() % 11]), 1.0));

	out.clear();

	for (int i = 0; i < n; i++)
	{
		lastRand.x = hammersley2d(i, n).x;

		lastRand.y = hammersley2d(i, n).y;

		tmp = glm::dvec3(2.0*lastRand.x - 1, 2.0*lastRand.y - 1, 2.0*lastRand.z - 1);

		//std::cout << glm::normalize(tmp).x << ", " << glm::normalize(tmp).y << ", " << glm::normalize(tmp).z << "\n";

		double theta = acos(2 * lastRand.y - 1);

		out.push_back(glm::dvec3(sin(theta) * cos(2 * lastRand.x*M_PI), sin(theta) * sin(2 * lastRand.x*M_PI), cos(theta)));
	}
}

glm::dvec2 importance_sample_ggx(double x, double y, double a)
{
	float phi = 2.0 * M_PI * x;
	float theta = acos(sqrt((1.0 - y)/((a*a - 1.0) * y + 1.0)));
	return glm::dvec2(phi, theta);
}
