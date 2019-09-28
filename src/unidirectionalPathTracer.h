
#pragma once


#include "content.h"
#include "sampling.h"
#include <boost/thread.hpp>

#define CAMERA_RAY 0
#define DIFFUSE_RAY 1
#define SPECULAR_RAY 2


struct tile
{
	int l, r, b, t;
	int w, h;
};

class graphDatum
{
public:
	graphDatum(float, float, float, float, float, float);
	float x, y, z, r, g, b;
};

class udp_renderer
{
	public:
		udp_renderer(wa_content *, float *, int, int, int, int, wa_shader, int, float, float);
		~udp_renderer();
		
		void render(float *);
		std::vector<wa_path>& getRecordedPaths();
		std::vector<wa_line>& getDebugLines();

		float importanceSampleRayMarchLength_attenuation(float, float, float *);
		float importanceSampleRayMarchLength_lightProximity(float3, float3, float3, float *);
		float importanceSampleRayMarchLength_uniform(float, float *);

		void setRenderResolution(int);
		void reset();
		void setIterationCount(int);

		void abort();

		int generateTiles(int, tile *);
		void setRegion(int, int, int, int);
		
	private:
		void divideIntoTiles();
		void renderTile(int);
		
		
		wa_colour radiance(float3, float3, int, bool, wa_path *, int, wa_shader);
		wa_colour transmission(wa_object *, wa_shader, float3, float3, float3, float3, int, wa_path *, bool, wa_shader);
		wa_colour reflection(wa_object *, wa_shader, float3, float3, float3, float3, int, wa_path *, int, bool, bool, float, wa_shader);
		wa_colour directReflection(wa_object *, wa_shader, float3, float3, float3, wa_path *);
		wa_colour indirectReflection(wa_object *, wa_shader, float3, float3, float3, int, wa_path *, wa_shader);
		wa_colour volumeAttenuation(wa_shader, float3, float3, wa_path *);
		wa_colour volumeScattering(int, wa_shader, float3, float3, int, wa_path *, bool, float);
		wa_colour directLightingSingleScatter(float3, float3, int, int, wa_shader);
		wa_colour directLightHit(float3, float3, int, wa_shader, float3 *);
		wa_colour visibility(float3, float3, wa_path *, bool);

		float3 reflectVector(float3, float3);
		
		wa_colour facingRatio(wa_ray Wi, int depth);
		void lowerResolution(int, int *, int *);
	
		
		float *frameBuffer;
		wa_content *theContent;

		wa_camera *camera;
		std::vector<wa_light *> lights;

		const int numPixelSamples;

		const int numHemisphereSamples;
		const float reciprocal_numHemisphereSamples;

		const int numSphereSamples;
		const float reciprocal_numSphereSamples;

		const int numLightSamples;
		const float reciprocal_numLightSamples;

		const float3 worldBoundSphere_P;
		const float worldBoundSphere_rad;

		int numThreads;
		tile *tiles;
		float progress;
		boost::mutex mutex;
		int latestReportedProgress;
		
		wa_shader globalShader;

		std::vector<wa_path> recordedPaths;
		std::vector<wa_line> debugLines;

		float *internalBuf;
		int itCnt;

		int disW, disH;
		int renW, renH;
		int rf;

		bool abortRender;

		tile rr;

		int numLights;

		std::vector<graphDatum> graphData;

		int maxTraceDepth;
		float rayMarchDistanceDistributionExponent;
		float rayMarchDistanceSampleRate;
};


class wa_samplerBeamLightRadiance : public wa_samplerBeam
{
public:
	wa_samplerBeamLightRadiance(const float3&, const float3&, std::vector<wa_light *>&);
	virtual float generateSample(float *);
	~wa_samplerBeamLightRadiance();


private:
	const float3 D;							//The unit direction vector from A to B

	const int n;							//The number of light sources

	std::vector<float> c;					//The 'c' value for each light (the distance along the beam of the closest beam point to the light)
	std::vector<float> h;					//The 'h' value for each light (the closest distance from the beam to the light)
	std::vector<float> theta1;				// The 'theta1' value for each light, (angle between lines joining the light point to A and the closest beam point)
	std::vector<float> theta2;				// The 'theta2' value for each light, (angle between lines joining the light point to B and the closest beam point)

	std::vector<float> weight;				// The weight of the light (the probablilty of choosing it for consideration)

	float importanceSampleRayMarchLength_lightProximity(float, float, float, float);
	float importanceSampleRayMarchLength_lightProximity_pdf(float, float, float, float, float);
};