#include <glm/glm.hpp>
#include <random>

#ifndef M_PI
#define M_PI 3.1415926535897
#endif
#define MAX_ENTITIES_PER_LEAF 16
#define MIN_LEAF_SIZE .03125
#define MAX_SUBDIV_RATIO 0.85
#define EPSILON 0.0001
#define DUPLICATE_THRESHOLD 150
#define SHADOW_BIAS 0.001

//squared vector length
inline double vecLengthSquared(glm::dvec3 vec)
{
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

inline glm::dvec3 randomUnitVec()
{
	return glm::dvec3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
}
