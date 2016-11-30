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
#define MIN_SAMPLES 128
#define SAMPLES 128
#define FOCAL_BLUR 0


//squared vector length
inline double vecLengthSquared(glm::dvec3 vec)
{
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

inline glm::dvec3 randomVec()
{
	return glm::dvec3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
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
