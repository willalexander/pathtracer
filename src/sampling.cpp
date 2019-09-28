#include <iostream>

#include "sampling.h"

using namespace std;



wa_samplerHemisphere_uniform::wa_samplerHemisphere_uniform(float3 D_in) :
	D(D_in)
{
}

float3 wa_samplerHemisphere_uniform::generateSample(float *reciprocal_pdf)
{
	float3 result;								//The direction vector chosen as a sample
	float3 tangentA, tangentB;					//The azimuth axes at the base of the hemisphere to be sampled
	float theta, phi;							//Incident and azimuth angles of the sample direction wrt the unit hemisphere

	tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	tangentB = normalized(cross(tangentA, D));

	theta = acos(1.0 - randFloat());
	phi = randFloat() * 2.0 * M_PI;

	result = cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
	if (reciprocal_pdf != NULL) *reciprocal_pdf = 2.0 * M_PI;

	return result;
}


wa_samplerHemisphere_projectedSurfaceArea::wa_samplerHemisphere_projectedSurfaceArea(float3 D_in) :
	D(D_in)
{
}

float3 wa_samplerHemisphere_projectedSurfaceArea::generateSample(float *reciprocal_pdf)
{
	float3 result;								//The direction vector chosen as a sample
	float3 tangentA, tangentB;					//The azimuth axes at the base of the hemisphere to be sampled
	float theta, phi;							//Incident and azimuth angles of the sample direction wrt the unit hemisphere

	tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	tangentB = normalized(cross(tangentA, D));

	theta = asin(cbrt(randFloat()));
	phi = randFloat() * 2.0 * M_PI;

	result = cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
	if (reciprocal_pdf != NULL) *reciprocal_pdf = 2.0 * M_PI / (3.0 * sin(theta) * cos(theta));

	return result;
}



wa_samplerBeam::wa_samplerBeam(const float3& A_in, const float3& B_in) :
	A(A_in),
	B(B_in),
	len(mag(B - A))
{

}


wa_samplerBeamUniform::wa_samplerBeamUniform(const float3& A_in, const float3& B_in) :
	wa_samplerBeam(A_in, B_in)
{

}

float wa_samplerBeamUniform::generateSample(float *reciprocal_pdf)
{
	float random = randFloat();

	*reciprocal_pdf = len;
	return (random * len);
}


wa_samplerBeamExponentialAttenuation::wa_samplerBeamExponentialAttenuation(const float3& A_in, const float3& B_in, const float& st_in) :
	wa_samplerBeam(A_in, B_in),
	st(st_in)
{
}


float wa_samplerBeamExponentialAttenuation::generateSample(float *reciprocal_pdf)
{
	float random = randFloat();

	float result = (-1.0 / st) * log(1.0 - ((1.0 - exp(-1.0 * st * len)) * random));

	*reciprocal_pdf = 1.0 / ((st / (1.0 - exp(-1.0 * st * len))) * exp(-1.0 * st * result));

	return result;
}



wa_samplerBeamExponentialDistribution::wa_samplerBeamExponentialDistribution(const float3& A_in, const float3& B_in, const float& a_in) :
	wa_samplerBeam(A_in, B_in),
	a(a_in),
	recip_a(1.0 / a),
	exp_aL_minus_1(exp(a*len) - 1.0)
{
}


float wa_samplerBeamExponentialDistribution::generateSample(float *reciprocal_pdf)
{
	float random = randFloat();

	float result = recip_a * log(random * exp_aL_minus_1 + 1.0);

	*reciprocal_pdf = (exp_aL_minus_1) / (a * exp(a*result));     //          1.0 / ((st / (1.0 - exp(-1.0 * st * len))) * exp(-1.0 * st * result));

	return result;
}

wa_samplerBeamPowerDistribution::wa_samplerBeamPowerDistribution(const float3& A_in, const float3& B_in, const float& a_in) :
	wa_samplerBeam(A_in, B_in),
	a(a_in),
	aplus1((float)(a + 1)),
	L_aplus1(pow(len, aplus1)),
	powroot(1.0 / (float)(a + 1))
{
}

float wa_samplerBeamPowerDistribution::generateSample(float *reciprocal_pdf)
{
	float random = randFloat();

	float result = pow(random * L_aplus1, powroot); 

	*reciprocal_pdf = L_aplus1 / (aplus1 *  pow(result, a));

	return result;
}



