#pragma once

#include "fundamentals.h"
#include "sampling.h"

#define LIGHT_AREA_TYPE 0
#define LIGHT_SPHERE_TYPE 1

class wa_light
{
	public:

		wa_light(float, float, float);

		virtual float3 getP() const;
		virtual float3 getCentre() const;
		virtual int geom_numQuads() const = 0;
		virtual void geom_placeVertices(wa_vertex *) const = 0;
		virtual void geom_placeIndices(wa_quad *) const = 0;
		virtual wa_colour getPower() const = 0;
		virtual wa_colour getIrradiance() const = 0;

		virtual void posAndDirSample(float3 *, float3 *) = 0;

		virtual bool radianceSample(float3, float3, wa_colour *, float3 *, float *, float) = 0;
		virtual float pdf(float3, float3) = 0;
		virtual void print() = 0;
		virtual int type() = 0;

		const float3 P;
};


class wa_areaLight : wa_light
{
	public:
		wa_areaLight(float, float, float, float, float, float, float, float, float, float, float, float);

		virtual float3 getCentre() const;
		float3 getN() const;
		virtual wa_colour getPower() const;
		virtual wa_colour getIrradiance() const;

		virtual int geom_numQuads() const;
		virtual void geom_placeVertices(wa_vertex *) const;
		virtual void geom_placeIndices(wa_quad *) const;

		virtual bool radianceSample(float3, float3, wa_colour *, float3 *, float *, float);
		virtual float pdf(float3, float3);

		virtual void print();
		virtual int type() { return LIGHT_AREA_TYPE; }

		virtual void posAndDirSample(float3 *, float3 *);

	//private:
		const float3 U;
		const float3 V;
		const float3 N;
		const float surfaceArea;
		const wa_colour power;
		const wa_colour irradiance;
		const wa_colour radiance;

		float3 samplePos() const;
};

class wa_sphereLight : public wa_light
{
	public:
		wa_sphereLight(float, float, float, float, float, float, float);

		float getRad() const;
		virtual wa_colour getPower() const;
		virtual wa_colour getIrradiance() const;

		virtual int geom_numQuads() const;
		virtual void geom_placeVertices(wa_vertex *) const;
		virtual void geom_placeIndices(wa_quad *) const;

		void sample(float3, float3 *, float3 *) const;
		void sample(float3, float, float3 *) const;
		virtual bool radianceSample(float3, float3, wa_colour *, float3 *, float *, float);
		virtual float pdf(float3, float3);

		virtual void print();
		virtual int type() { return LIGHT_SPHERE_TYPE; }

		virtual void posAndDirSample(float3 *, float3 *);

		private:
			const float rad;
			const wa_colour power;
			const wa_colour irradiance;
			const wa_colour radiance;

			void samplePos(float3 *, float3 *) const;
};