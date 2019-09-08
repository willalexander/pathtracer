#include "objects.h"

using namespace std;

wa_shader::wa_shader()
{
	kd = 0.5;
	kt = 0.0;
	ks = 0.5;

	proceduralVolume = 1;
	homogenous = true;

	isOpaque = true;
	isClear = true;
	isLight = false;
}

void wa_shader::initialize()
{
	kd = 0.5;
	kt = 0.0;

	sa = 0.0;
	se = 0.0;
	ss = 0.0;
	st = 0.0;

	isOpaque = true;
	isClear = true;
	isLight = false;

	diffuse = false;
	specular = false;
	sss = false;
	dielectric = false;
	dielectricsss = false;
}

float wa_shader::fresnel(float3 Ng, float3 D)
{
	/*//If a pure diffuse surface, set fresnel to 0.0:
	if ((ks.r == 0.0) && (ks.g == 0.0) && (ks.b == 0.0)) return 0.0;

	//If a pure specular surface, set fresnel to 1.0:
	if ((kd.r == 0.0) && (kd.g == 0.0) && (kd.b == 0.0)) return 0.0;

	float facingRatio = dot(N, D);										// Facing ratio of surface
	float fresnelFactor = pow(1.0 - facingRatio, 10.0);

	return (1.0 - fresnelFactor) * 0.05 + fresnelFactor * 1.0;*/

	return 0.0;
}

float wa_shader::BRDF(float3 Wo, float3 Wi)
{
	return (1.0 / M_PI);
}

float wa_shader::phase(float angle)
{
	return (1.0 / (4.0 * M_PI));
}

float3 wa_shader::refract(float3 Ng, float3 D, float *fresnel)
{
	float NdotD = dot(Ng,D);
	float3 tangent = normalized(D -NdotD*Ng);
	float cosTheta1 = NdotD;
	float sinTheta1 = sqrt(1.0 - pow(NdotD, 2));
	float sinTheta2, cosTheta2;

	//If we're entering the surface:
	if (NdotD > 0.0)
	{
		sinTheta2 = sinTheta1/ior;
		cosTheta2 = sqrt(1.0 - pow(sinTheta2, 2));

		*fresnel = 0.5*(pow((ior*cosTheta1 - cosTheta2) / (ior*cosTheta1 + cosTheta2), 2) + pow((cosTheta1 - ior * cosTheta2) / (cosTheta1 + ior * cosTheta2), 2));
		
		//If the ray is exactly perpendicular to the surface, this can lead to a division by zero. Set the fresnel value directly:
		if (std::isnan(*fresnel)) *fresnel = 0.0;

		return (sinTheta2 * tangent + cosTheta2 * Ng);
	}

	//If we're exiting the surface:
	sinTheta2 = sinTheta1 * ior;

	//Total internal reflection:
	if (sinTheta2 > 1.0)
	{
		*fresnel = 1.0;
		return 0.0;
	}

	cosTheta2 = sqrt(1.0 - pow(sinTheta2, 2));
	cosTheta1 *= -1.0;

	*fresnel = 0.5*(pow((cosTheta1 - ior*cosTheta2) / (cosTheta1 + ior*cosTheta2), 2) + pow((ior* - cosTheta2) / (ior*cosTheta1 + cosTheta2), 2));
	
	//If the ray is exactly perpendicular to the surface, this can lead to a division by zero. Set the fresnel value directly:
	if (std::isnan(*fresnel)) *fresnel = 0.0;

	return (sinTheta2 * tangent - 1.0*cosTheta2 * Ng);
}

wa_colour wa_shader::st_f(float3 P)
{
	if(homogenous) return st;

	//Linear transmittance from almost 0 to almost 1 along the x axis:
	wa_colour result = -1.0 * log((5.0 / 6.0) * X(P) + 0.5);

	return result;
}

wa_colour wa_shader::minimumExtinctionCoefficient()
{
	//Linear transmittance from 0 to 1 along the x axis:
	return 0.087;
}

wa_colour wa_shader::maximumExtinctionCoefficient()
{
	//Linear transmittance from 0 to 1 along the x axis:
	return 2.48;
}






wa_object::wa_object(string fileName_in, string shaderTypeName, wa_shader shader_in) : 
		fileName(fileName_in)
{
	if(shaderTypeName == "diffuse") shaderType = SHADER_TYPE_DIFFUSE;
	if(shaderTypeName == "sss") shaderType = SHADER_TYPE_SSS;

	shader = shader_in;
}

int wa_object::getID()
{
	return ID;
}

void wa_object::setID(int ID_in)
{
	ID = ID_in;
}

string wa_object::getFileName()
{
	return fileName;
}

int wa_object::getShaderType()
{
	return shaderType;
}

wa_shader wa_object::getShader()
{
	return shader;
}

void wa_object::setShader(wa_shader shader_in)
{
	shader = shader_in;
}