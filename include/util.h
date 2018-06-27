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
#define MAX_PHOTONS_PER_LEAF 16
#define MIN_LEAF_SIZE .0015
#define MAX_SUBDIV_RATIO 0.75
#define EPSILON 0.00001
#define DUPLICATE_THRESHOLD 150
#define SHADOW_BIAS 0.0001
#define AA_JITTER 1.2
#define MIN_DEPTH 4
#define MAX_DEPTH 32
#define NOISE_THRESH 0.0015
#define MIN_SAMPLES 8
#define SAMPLES 32
#define PHOTONS 75000
#define PHOTON_DEPTH 5
#define RAYMARCH_STEPSIZE 0.04
#define FOCAL_BLUR 0
#define GAMMA 2.2


//squared vector length
inline double vecLengthSquared(const glm::dvec3& vec)
{
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

inline glm::dvec3 fastNorm(glm::dvec3& vec)
{
	double l = vecLengthSquared(vec);

	return l == 1 ? vec : vec / sqrt(l);
}

inline double compMax(const glm::dvec3& vec)
{
	return std::max(std::max(vec.x, vec.y), vec.z);
}

inline double drand()
{
	return (double)rand() / RAND_MAX;
}

inline double clamp(double val, double min, double max)
{
	if (val < min)
		return min;

	if (val > max)
		return max;

	return val;
}

//gamma correction
inline glm::dvec3 gamma(const glm::dvec3& col, double g)
{
	return glm::dvec3(std::pow(col.x, 1.0 / g), std::pow(col.y, 1.0 / g), std::pow(col.z, 1.0 / g));
}

//faster pow approximation for small exponents, based on http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
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

inline double fastPrecisePow(double a, double b)
{
	// calculate approximation with fraction of the exponent
	int e = (int)b;
	union {
		double d;
		int x[2];
	} u = { a };
	u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
	u.x[0] = 0;

	// exponentiation by squaring with the exponent's integer part
	// double r = u.d makes everything much slower, not sure why
	double r = 1.0;
	while (e) {
		if (e & 1) {
			r *= a;
		}
		a *= a;
		e >>= 1;
	}

	return r * u.d;
}

inline double fastPrecisePow(double a, int b)
{
	int e = b;

	// exponentiation by squaring with the exponent's integer part
	// double r = u.d makes everything much slower, not sure why
	double r = 1.0;
	while (e) {
		if (e & 1) {
			r *= a;
		}
		a *= a;
		e >>= 1;
	}

	return r;
}

inline double fastLength(glm::dvec3& vec)
{
	return fastPow(vecLengthSquared(vec), 0.5);
}

//computes the binary inverse of a float number
inline float radicalInverse_VdC(unsigned int bits)
{
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

//computes refraction and total internal reflection for an input direction
inline glm::dvec3 refr(const glm::dvec3& inc, const glm::dvec3& norm, double eta)
{
	double d = glm::dot(norm, inc);
	double k = 1.0 - eta * eta * (1.0 - d * d);
    if (k < EPSILON)
        return glm::reflect(inc, norm);
    else
        return eta * inc - (eta * d + sqrt(k)) * norm;
}

inline glm::dvec3 randomUnitVec(double x, double y)
{
	double theta = acos(2 * y - 1);

	return glm::dvec3(sin(theta) * cos(2 * x*M_PI), sin(theta) * sin(2 * x*M_PI), cos(theta));
}

inline void intersectSIMD(float* __restrict t0, float* __restrict t1, const float* boxes, const float* ray, const float* invDir, const float* invlz)
{
	//float t[48]; //layout: [min1.x, min1.y, min1.z, max1.x, ...]

	for (int i = 0; i < 24; i++)
	{
		float tmp0 = (boxes[i] - ray[i])*invDir[i]; //t0
		float tmp1 = (boxes[i + 24] - ray[i])*invDir[i]; //t1

		t0[i] = invlz[i] * tmp1 + (1.0f - invlz[i])*tmp0;
		t1[i] = invlz[i] * tmp0 + (1.0f - invlz[i])*tmp1;
	}

}


glm::dvec3 randomVec();

// generates a point of the hammersley point set, based on http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
glm::dvec2 hammersley2d(unsigned int i, unsigned int N);

glm::dvec3 hemisphereSample_uniform(float u, float v);

glm::dvec3 hemisphereSample_cos(float u, float v, double power);

glm::dvec3 hemisphereSample_cos(glm::dvec3 normal, float u, float v, double power);

glm::dvec3 sphereCapSample_cos(glm::dvec3 normal, float u, float v, double power, double frac);

double PowerCosHemispherePdfW(const glm::dvec3  aNormal, glm::dvec3  aDirection, double  aPower);

glm::dvec3 sample_phong(const glm::dvec3 &outdir, const glm::dvec3 &n, double power, double sx, double sy);

glm::dvec2 importance_sample_ggx(double x, double y, double a);

//generates a sequence of subrandom numbers
void subrand(std::vector<double>& out, int n);

//generates a sequence of subrandom unit vectors
void subrandUnitVec(std::vector<glm::dvec3>& out, int n);

bool triBoxOverlap(glm::dvec3 boxcenter, glm::dvec3 boxhalfsize, glm::dvec3 triverts[3]);

#endif
