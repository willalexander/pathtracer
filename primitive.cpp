#include <iostream>

#include "primitive.h"



wa_primitive::wa_primitive()
{
}

wa_primitive::wa_primitive(float3 Pval, float3 Nval)
{
	P = Pval;
	N = Nval;
}

wa_primitive::~wa_primitive()
{
}

wa_colour wa_primitive::Le(float3 Wo)
{
	if (lightOrSurface == 0)
	{
		return (irradiance / M_PI);
	}

	else return wa_colour(0.0);
}

wa_colour wa_primitive::BRDF(float3  Wi, float3 Wo)
{
	return wa_colour((1.0 / M_PI));
}

wa_object *wa_primitive::getObject()
{
	return objectPointer;
}

void wa_primitive::setObject(wa_object *objectPointer_in)
{
	objectPointer = objectPointer_in;
}

float3 wa_primitive::getP()
{
	return P;
}

float3 wa_primitive::getN()
{
	return N;
}

void wa_primitive::setP(float3 val)
{
	P = val;
}

void wa_primitive::setN(float3 val)
{
	N = val;
}

