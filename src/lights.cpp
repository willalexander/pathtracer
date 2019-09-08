#include <iostream>

#include "lights.h"

using namespace std;


wa_light::wa_light(float pxVal, float pyVal, float pzVal) :
	P(pxVal, pyVal, pzVal)
{
}

float3 wa_light::getP() const
{
	return P;
}

float3 wa_light::getCentre() const
{
	return P;
}


wa_areaLight::wa_areaLight(float px, float py, float pz, float ux, float uy, float uz, float vx, float vy, float vz, float cr, float cg, float cb) :
	wa_light(px, py, pz),
	U(ux, uy, uz),
	V(vx, vy, vz),
	N(normalized(cross(U, V))),
	surfaceArea(sqrt(ux * ux + uy * uy + uz * uz) * sqrt(vx * vx + vy * vy + vz * vz)),
	power(cr, cg, cb),
	irradiance(0.5 * power / surfaceArea),
	radiance(irradiance / M_PI)
{
}

float3 wa_areaLight::getCentre() const
{
	return P + 0.5*(U + V);
}

int wa_areaLight::geom_numQuads() const
{
	return 1;
}

void wa_areaLight::geom_placeVertices(wa_vertex *vertBuf) const
{
	float3 PU = P + U;
	float3 PV = P + V;
	float3 PUV = P + U + V;

	vertBuf[0].x = X(P); vertBuf[0].y = Y(P); vertBuf[0].z = Z(P);
	vertBuf[1].x = X(PU); vertBuf[1].y = Y(PU); vertBuf[1].z = Z(PU);
	vertBuf[2].x = X(PUV); vertBuf[2].y = Y(PUV); vertBuf[2].z = Z(PUV);
	vertBuf[3].x = X(PV); vertBuf[3].y = Y(PV); vertBuf[3].z = Z(PV);
}

void wa_areaLight::geom_placeIndices(wa_quad *quadBuf) const
{
	quadBuf[0].v0 = 0; quadBuf[0].v1 = 1; quadBuf[0].v2 = 2; quadBuf[0].v3 = 3;
}


float3 wa_areaLight::getN() const
{
	return N;
}

wa_colour wa_areaLight::getPower() const
{
	return power;
}

wa_colour wa_areaLight::getIrradiance() const
{
	return irradiance;
}

float3 wa_areaLight::samplePos() const
{
	return P + randFloat() * U + randFloat() * V;
}

bool wa_areaLight::radianceSample(float3 SP, float3 SN, wa_colour *Li, float3 *LP, float *reciprocal_pdf, float traceBias)
{
	*LP = samplePos();

	//If the radiance is being received by a surface, and this light sample point is not 'in front' of the surface, then it contributes nothing:
	//if ((SN != 0.0) && (dot(SN, (*LP-SP)) <= 0.0)) return false;

	//Radiance of this area light is constant across its surface:
	*Li = radiance;		
	
	//pdf is measured relative to solid angle at the surface, not relative to this light's surface area:
	*reciprocal_pdf = surfaceArea * surfaceArea_to_solidAngle(SP, *LP, N);

	return true;
}

void wa_areaLight::print()
{
	cout << "Area Light." << endl;
	cout << "P: " << P << endl;
	cout << "U: " << U << endl;
	cout << "V: " << V << endl;
	cout << "N: " << N << endl;
	cout << "power: " << power << endl;
}

void wa_areaLight::posAndDirSample(float3 *Ps, float3 *Ds)
{
	wa_sampler *sampler = (wa_sampler *)(new wa_samplerHemisphere_projectedSurfaceArea({ 0.0, 1.0, 0.0 }));

	*Ps = samplePos();		
	*Ds = sampler->generateSample(NULL);

	//Randomly choose the top or bottom face of the area light:
	if (randBool()) *Ds *= -1.0;
}





wa_sphereLight::wa_sphereLight(float px, float py, float pz, float radVal, float cr, float cg, float cb) :
	wa_light(px, py, pz),
	rad(radVal),
	power(cr, cg, cb),
	irradiance(power / (4.0 * M_PI * pow(rad, 2))),
	radiance(irradiance / M_PI)
{
}

int wa_sphereLight::geom_numQuads() const
{
	return 200;
}

void wa_sphereLight::geom_placeVertices(wa_vertex *vertBuf) const
{
	wa_vertex v0, v1, v2, v3;
	int quadNum;

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			quadNum = i * 20 + j;
			generateSpherePoints(P, rad, 10, 20, i, j, &v0, &v1, &v2, &v3);
			vertBuf[quadNum * 4 + 0] = v0;
			vertBuf[quadNum * 4 + 1] = v1;
			vertBuf[quadNum * 4 + 2] = v2;
			vertBuf[quadNum * 4 + 3] = v3;
		}	
	}
}

void wa_sphereLight::geom_placeIndices(wa_quad *quadBuf) const
{
	for (int i = 0; i < 200; i++)
	{
		quadBuf[i].v0 = i * 4 + 0;
		quadBuf[i].v1 = i * 4 + 1;
		quadBuf[i].v2 = i * 4 + 2;
		quadBuf[i].v3 = i * 4 + 3;
	}
}

float wa_sphereLight::getRad() const
{
	return rad;
}

wa_colour wa_sphereLight::getPower() const
{
	return power;
}

wa_colour wa_sphereLight::getIrradiance() const
{
	return irradiance;
}


//'sample' will generate a sample on one hemishpere of the light. The caller must provide the vector direction of the centre of the hemisphere:
void wa_sphereLight::sample(float3 D, float3 *sampleP, float3 *sampleN) const
{
	*sampleP = P + rad * hemisphereSampleUniform(D, NULL);
	*sampleN = normalized(*sampleP - P);
}

//'samplePos' will generate a sample on the surface of the light. 
void wa_sphereLight::samplePos(float3 *sampleP, float3 *sampleN) const
{
	*sampleP = P + rad * sphereSampleUniform();
	*sampleN = normalized(*sampleP - P);
}

//'sample' will generate a sample in the cone provided. The caller must provide the vector direction of the centre of the cone, and the cone angle:
void wa_sphereLight::sample(float3 D, float coneAngle, float3 *sampleD) const
{
	*sampleD = hemisphereSampleUniform(D, coneAngle, NULL);
}


bool wa_sphereLight::radianceSample(float3 SP, float3 SN, wa_colour *Li, float3 *LP, float *reciprocal_pdf, float traceBias)
{
	float3 SL = P - SP;							// Direction vector from surface point to centre of the light
	float SLdist = mag(SL);						// Distance from surface point to centre of light
	float3 SLn = normalized(SL);				// Unit vector in direction from surface point to light point
	float coneAngle = asin(rad / SLdist);		// The half-angle of the cone formed by this sphere light and the surface point
	float3 D;									// Direction vector from the surface point to a sample point on the surface of the light
	float dist;									// Distance from the surface point to the sample point on the light surface

	//If the surface point is inside or on the light source, then no contribution:
	if (SLdist < (rad + traceBias)) return false;

	//We sample only the cone subtended by the light at the shading point:
	sample(SLn, coneAngle, &D);

	//If the radiance is being received by a surface, and this light sample point is not 'in front' of the surface, then it contributes nothing:
	if ((SN != 0.0) && (dot(SN, D) <= 0.0)) return false;

	//Compute the corresponding point on the light surface:
	if (solveQuadraticEquation(1.0, -2.0 * dot(SL, D), pow(mag(SL), 2) - pow(rad, 2), NULL, &dist) == false) dist = mag(SL) * cos(coneAngle);
	*LP = SP + dist * D;

	//Radiance is constant all over the surface:
	*Li = radiance;

	//The pdf is relative to the solid angle at the surface, not relative to this light's surface area:
	*reciprocal_pdf = 2.0 * M_PI * (1.0 - cos(coneAngle));

	return true;
}


void wa_sphereLight::print()
{
	cout << "Sphere Light." << endl;
	cout << "P: " << P << endl;
	cout << "rad: " << rad << endl;
	cout << "power: " << power << endl;
}


void wa_sphereLight::posAndDirSample(float3 *Ps, float3 *Ds)
{
	float3 Ns;									//Normal vector at randomly sampled point on surface
	wa_sampler *sampler;						//Sampler object for sampling in a hemisphere
		
	samplePos(Ps, &Ns);
	sampler = (wa_sampler *)(new wa_samplerHemisphere_projectedSurfaceArea(Ns));

	*Ds = sampler->generateSample(NULL);
}
