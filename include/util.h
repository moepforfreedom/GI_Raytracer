#include <glm/glm.hpp>
#include <random>
#include <iostream>

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
#define SAMPLES 1
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


//generates a sequence of subrandom numbers
inline std::vector<double> subrand(int n)
{
	double lastRand = drand();
	double a = drand();

	std::vector<double> res;
	
	for (int i = 0; i < n; i++)
	{
		lastRand = fmod(lastRand + a, 1.0);

		res.push_back(lastRand);		
	}

	return res;
}
