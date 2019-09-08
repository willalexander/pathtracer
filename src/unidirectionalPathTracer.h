
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


class udp_renderer
{
	public:
		udp_renderer(wa_content *, float *, int, int, int, int);
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
};

