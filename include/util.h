#ifndef _UTILS
#define _UTILS

#include <glm/glm.hpp>
#include <random>
#include <iostream>
#include <ctime>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.1415926535897
#endif
#define MAX_ENTITIES_PER_LEAF 16
#define MIN_LEAF_SIZE .03125
#define MAX_SUBDIV_RATIO 0.85
#define EPSILON 0.0001
#define DUPLICATE_THRESHOLD 150
#define SHADOW_BIAS 0.002
#define AA_JITTER 1
#define MAX_DEPTH 4
#define NOISE_THRESH 0.005
#define MIN_SAMPLES 8
#define SAMPLES 256
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

inline int clamp(int min, int max, int val)
{
	if (val < min)
		return min;

	if (val > max)
		return max;

	return val;
}

inline float radicalInverse_VdC(unsigned int bits)
{
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

glm::dvec3 randomVec();

glm::dvec3 randomUnitVec();

// generates a point of the hammersley point set, based on http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
glm::dvec2 hammersley2d(unsigned int i, unsigned int N);

glm::dvec3 hemisphereSample_uniform(float u, float v);

glm::dvec3 hemisphereSample_cos(float u, float v, double power);

double PowerCosHemispherePdfW(const glm::dvec3  aNormal, glm::dvec3  aDirection, double  aPower);

glm::dvec2 importance_sample_ggx(double x, double y, double a);

//generates a sequence of subrandom numbers
void subrand(std::vector<double>&out, int n);

//generates a sequence of subrandom unit vectors
void subrandUnitVec(std::vector<glm::dvec3>&out, int n);

#endif
