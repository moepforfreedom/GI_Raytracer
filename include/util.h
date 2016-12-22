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
#define MAX_SUBDIV_RATIO 0.75
#define EPSILON 0.0001
#define DUPLICATE_THRESHOLD 150
#define SHADOW_BIAS 0.002
#define AA_JITTER 1
#define MAX_DEPTH 8
#define NOISE_THRESH 0.004
#define MIN_SAMPLES 4
#define SAMPLES 32
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

//faster pow, based on http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b)
{
  union
  {
    double d;
    int x[2];
  } u = { a };

  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
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

inline glm::dvec3 refr(glm::dvec3 inc, glm::dvec3 norm, double eta)
{
	double d = glm::dot(norm, inc);
	double k = 1.0 - eta * eta * (1.0 - d * d);
    if (k < 0.000001)
        return glm::reflect(inc, norm);
    else
        return eta * inc - (eta * d + sqrt(k)) * norm;
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
