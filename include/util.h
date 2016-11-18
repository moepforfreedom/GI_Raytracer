#include <glm/glm.hpp>
#ifndef M_PI
#define M_PI 3.1415926535897
#endif
#define MAX_ENTITIES_PER_LEAF 32
#define MIN_LEAF_SIZE .03125
#define EPSILON 0.0001

//squared vector length
inline double vecLengthSquared(glm::dvec3 vec)
{
	return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}


/*BoundingBox MeshBbox(std::vector<triangle*> tris)
{
	glm::dvec3 min(INFINITY, INFINITY, INFINITY), max(-INFINITY, -INFINITY, -INFINITY); //TODO: change this somehow

	for (int i = 0; i < tris.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			glm::dvec3 pos = tris[i]->vertices[j]->pos;

			if (pos.x < min.x)
				min.x = pos.x;
			if (pos.x > max.x)
				max.x = pos.x;

			if (pos.y < min.y)
				min.y = pos.y;
			if (pos.y > max.y)
				max.y = pos.y;

			if (pos.z < min.z)
				min.z = pos.z;
			if (pos.z > max.z)
				max.z = pos.z;
		}
	}

	return BoundingBox(min, max);
}*/
