#include <glm/glm.hpp>
#include <random>
#include <iostream>
#include <ctime>
#include "util.h"


glm::dvec3 randomVec()
{
	return glm::dvec3((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
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
	float cosTheta = fastPrecisePow(1.0f - u, (1.0f / power));
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	return glm::dvec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

glm::dvec3 hemisphereSample_cos(glm::dvec3 normal, float u, float v, double power)
{

	double z = std::abs(normal.z);

	glm::dmat3x3 rot(z + (1.0 / (1 + z))*-normal.y*-normal.y, (1.0 / (1 + z))*(normal.x*-normal.y), -normal.x,
		(1.0 / (1 + z))*(normal.x*-normal.y), z + (1.0 / (1 + z))*-normal.x*-normal.x, -normal.y,
		normal.x, normal.y, z);


	float phi = v * 2.0f * M_PI;
	float cosTheta = fastPrecisePow(1.0f - u, (1.0f / power));
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	glm::dvec3 res(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

	res = rot*res;

	if (normal.z < 0)
	{
		res.z *= -1.0;
	}

	return res;
}

glm::dvec3 sphereCapSample_cos(glm::dvec3 normal, float u, float v, double power, double frac)
{

	double z = std::abs(normal.z);

	glm::dmat3x3 rot(z + (1.0 / (1 + z))*-normal.y*-normal.y, (1.0 / (1 + z))*(normal.x*-normal.y), -normal.x,
		(1.0 / (1 + z))*(normal.x*-normal.y), z + (1.0 / (1 + z))*-normal.x*-normal.x, -normal.y,
		normal.x, normal.y, z);


	float phi = v * 2.0f * M_PI;
	float cosTheta = frac*fastPrecisePow(1.0f - u, (1.0f / power))+(1-frac);
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	glm::dvec3 res(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

	res = rot*res;

	if (normal.z < 0)
	{
		res.z *= -1.0;
	}

	return res;
}

double PowerCosHemispherePdfW(glm::dvec3  aNormal, glm::dvec3  aDirection, double  aPower)
{
	float cosTheta = std::max(0.0, glm::dot(aNormal, aDirection));
	return (aPower + 1.0) * std::pow(cosTheta, aPower) * ((1.0 / M_PI) * 0.5);
}

glm::dvec3 sample_phong(const glm::dvec3 &outdir, const glm::dvec3 &n, double power, double sx, double sy)
{
	double z = std::abs(outdir.z);

	glm::dmat3x3 rot(z + (1.0 / (1 + z))*-outdir.y*-outdir.y, (1.0 / (1 + z))*(outdir.x*-outdir.y), -outdir.x,
		(1.0 / (1 + z))*(outdir.x*-outdir.y), z + (1.0 / (1 + z))*-outdir.x*-outdir.x, -outdir.y,
		outdir.x, outdir.y, z);

	glm::dvec3 out = rot*hemisphereSample_cos(sx, sy, power);

	if (outdir.z < 0)
	{
		out.z *= -1.0;
	}

	return out;
}

//generates a sequence of subrandom numbers using additive recurrence
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



/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-M�ller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

int planeBoxOverlap(glm::dvec3 normal, float d, glm::dvec3 maxbox)
{
	int q;
	glm::dvec3 vmin, vmax;
	for (q = 0; q <= 2; q++)
	{
		if (normal[q]>0.0f)
		{
			vmin[q] = -maxbox[q];
			vmax[q] = maxbox[q];
		}
		else
		{
			vmin[q] = maxbox[q];
			vmax[q] = -maxbox[q];
		}
	}
	if (glm::dot(normal, vmin) + d>0.0f) return 0;
	if (glm::dot(normal, vmax) + d >= 0.0f) return 1;

	return 0;
}


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)             \
    p0 = a*v0[1] - b*v0[2];                    \
    p2 = a*v2[1] - b*v2[2];                    \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)              \
    p0 = a*v0[1] - b*v0[2];                    \
    p1 = a*v1[1] - b*v1[2];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2];   \
    if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)             \
    p0 = -a*v0[0] + b*v0[2];                   \
    p2 = -a*v2[0] + b*v2[2];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)              \
    p0 = -a*v0[0] + b*v0[2];                   \
    p1 = -a*v1[0] + b*v1[2];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2];   \
    if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)             \
    p1 = a*v1[0] - b*v1[1];                    \
    p2 = a*v2[0] - b*v2[1];                    \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)              \
    p0 = a*v0[0] - b*v0[1];                \
    p1 = a*v1[0] - b*v1[1];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1];   \
    if(min>rad || max<-rad) return 0;

bool triBoxOverlap(glm::dvec3 boxcenter, glm::dvec3 boxhalfsize, glm::dvec3 triverts[3])
{

	/*    use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */
	glm::dvec3 v0, v1, v2;
	float min, max, d, p0, p1, p2, rad, fex, fey, fez;
	glm::dvec3 normal, e0, e1, e2;

	// move everything so that the boxcenter is in (0,0,0)
	v0 = triverts[0] - boxcenter;
	v1 = triverts[1] - boxcenter;
	v2 = triverts[2] - boxcenter;

	/* compute triangle edges */
	e0 = v1 - v0;      /* tri edge 0 */
	e1 = v2 - v1;      /* tri edge 1 */
	e2 = v0 - v2;      /* tri edge 2 */


	//  test the 9 tests first (this was faster)
	fex = std::abs(e0[0]);
	fey = std::abs(e0[1]);
	fez = std::abs(e0[2]);
	AXISTEST_X01(e0[2], e0[1], fez, fey);
	AXISTEST_Y02(e0[2], e0[0], fez, fex);
	AXISTEST_Z12(e0[1], e0[0], fey, fex);

	fex = std::abs(e1[0]);
	fey = std::abs(e1[1]);
	fez = std::abs(e1[2]);
	AXISTEST_X01(e1[2], e1[1], fez, fey);
	AXISTEST_Y02(e1[2], e1[0], fez, fex);
	AXISTEST_Z0(e1[1], e1[0], fey, fex);

	fex = std::abs(e2[0]);
	fey = std::abs(e2[1]);
	fez = std::abs(e2[2]);
	AXISTEST_X2(e2[2], e2[1], fez, fey);
	AXISTEST_Y1(e2[2], e2[0], fez, fex);
	AXISTEST_Z12(e2[1], e2[0], fey, fex);

	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */
	FINDMINMAX(v0[0], v1[0], v2[0], min, max);
	if (min>boxhalfsize[0] || max<-boxhalfsize[0]) return false;

	/* test in Y-direction */
	FINDMINMAX(v0[1], v1[1], v2[1], min, max);
	if (min>boxhalfsize[1] || max<-boxhalfsize[1]) return false;

	/* test in Z-direction */
	FINDMINMAX(v0[2], v1[2], v2[2], min, max);
	if (min>boxhalfsize[2] || max<-boxhalfsize[2]) return false;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */
	normal = glm::cross(e0, e1);
	d = -glm::dot(normal, v0);  /* plane eq: normal.x+d=0 */
	if (!planeBoxOverlap(normal, d, boxhalfsize)) return 0;

	return true;
}
