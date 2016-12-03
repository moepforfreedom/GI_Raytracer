#include <glm/glm.hpp>
#include <random>
#include <iostream>
#include <ctime>

#ifndef M_PI
#define M_PI 3.1415926535897
#endif
#define MAX_ENTITIES_PER_LEAF 16
#define MIN_LEAF_SIZE .03125
#define MAX_SUBDIV_RATIO 0.85
#define EPSILON 0.0001
#define DUPLICATE_THRESHOLD 150
#define SHADOW_BIAS 0.001
#define AA_JITTER 1
#define MAX_DEPTH 1
#define NOISE_THRESH 0.003
#define MIN_SAMPLES 8
#define SAMPLES 128
#define FOCAL_BLUR 0


//squared vector length
inline double vecLengthSquared(glm::dvec3 vec)
{
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

inline double drand()
{
	return (double)rand() / RAND_MAX;
}

inline glm::dvec3 randomVec()
{
	return glm::dvec3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
}

inline glm::dvec3 randomUnitVec()
{
	return glm::normalize(glm::dvec3(2*drand()-1, 2*drand()-1, 2*drand()-1));
}

inline int clamp(int min, int max, int val)
{
	if (val < min)
		return min;

	if (val > max)
		return max;

	return val;
}


inline float radicalInverse_VdC(uint bits) 
{
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

// generates a point of the hammersley point set, based on http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
inline glm::dvec2 hammersley2d(uint i, uint N) 
{
	return glm::dvec2(float(i) / float(N), radicalInverse_VdC(i));
}

inline glm::dvec3 hemisphereSample_uniform(float u, float v) 
{
	float phi = v * 2.0 * M_PI;
	float cosTheta = 1.0 - u;
	float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

inline glm::dvec3 hemisphereSample_cos(float u, float v) 
{
	float phi = v * 2.0 * M_PI;
	float cosTheta = sqrt(1.0 - u);
	float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

inline double PowerCosHemispherePdfW(const glm::dvec3  aNormal, glm::dvec3  aDirection, double  aPower)
{
	float cosTheta = std::max(0.0, glm::dot(aNormal, aDirection));
	return (aPower + 1.0) * std::pow(cosTheta, aPower) * ((1.0/M_PI) * 0.5);
}

//generates a sequence of subrandom numbers
inline void subrand(std::vector<double>&out, int n)
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
inline void subrandUnitVec(std::vector<glm::dvec3>&out, int n)
{
	double primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

	glm::dvec3 lastRand = randomVec();

	glm::dvec3 tmp;

	glm::dvec3 a = glm::dvec3(fmod(sqrt(primes[rand() % 11]), 1.0), fmod(sqrt(primes[rand() % 11]), 1.0), fmod(sqrt(primes[rand() % 11]), 1.0));

	out.clear();

	for (int i = 0; i < n; i++)
	{
		//lastRand.x = fmod(lastRand.x + a.x, 1.0);

		//lastRand.y = fmod(lastRand.y + a.y, 1.0);

		//lastRand.z = fmod(lastRand.z + a.z, 1.0);

		lastRand.x = hammersley2d(i, n).x;
		
		lastRand.y = hammersley2d(i, n).y;

		tmp = glm::dvec3(2.0*lastRand.x - 1, 2.0*lastRand.y - 1, 2.0*lastRand.z - 1);

		//std::cout << glm::normalize(tmp).x << ", " << glm::normalize(tmp).y << ", " << glm::normalize(tmp).z << "\n";

		double theta = acos(2 * lastRand.y - 1);

		out.push_back(glm::dvec3(sin(theta) * cos(2 * lastRand.x*M_PI), sin(theta) * sin(2 * lastRand.x*M_PI), cos(theta)));
	}
}
