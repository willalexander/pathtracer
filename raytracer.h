#pragma once

#include <iostream>
#include <string>

#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>


using namespace std;


class wa_raytracer
{
public:

	void init();

	void loadLightGeo(float, float, float, float, float, float, float, float, float);
	void loadGeo(string);

	void updateLight(int, float, float, float, float, float, float, float, float, float);

	void commit();

	bool trace(float, float, float, float, float, float, float *, float *, float *, float *, float *, float *, int *, float) const;

	RTCScene getEmbreeScene();


	void release();


private:

	RTCDevice device;
	RTCScene scene;
};