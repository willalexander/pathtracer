#pragma once

#include "fundamentals.h"

#define WEIGHT_EQUALITY 1


class wa_sampler
{
	public:
		virtual float3 generateSample(float *) = 0;
};

class wa_samplerHemisphere_uniform : public wa_sampler
{
	public:
		wa_samplerHemisphere_uniform(float3);
		virtual float3 generateSample(float *);

	private:
		const float3 D;
};

class wa_samplerHemisphere_projectedSurfaceArea : public wa_sampler
{
public:
	wa_samplerHemisphere_projectedSurfaceArea(float3);
	virtual float3 generateSample(float *);

private:
	const float3 D;
};


class wa_samplerBeam
{
	public:
		wa_samplerBeam(const float3&, const float3&);
		virtual float generateSample(float *) = 0;

		const float3 A;							// The start point of the beam line
		const float3 B;							// The end point of the beam line
		const float len;						// The distance from A to B
};


class wa_samplerBeamUniform : public wa_samplerBeam
{
	public:
		wa_samplerBeamUniform(const float3&, const float3&);
		virtual float generateSample(float *);
};

class wa_samplerBeamExponentialAttenuation : public wa_samplerBeam
{
	public:
		wa_samplerBeamExponentialAttenuation(const float3&, const float3&, const float&);
		virtual float generateSample(float *);

	private:
		const float st;
};
