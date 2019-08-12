
#pragma once

#include <iostream>
#include <cmath>
#include <string>

#include "raytracer.h"
#include "lights.h"
#include "primitive.h"

using namespace std;

#define RAYTRACER_EMBREE 1


class wa_content
{
	public:
		wa_content() {}
		wa_content(wa_camera, std::vector<wa_light *>, std::vector<wa_object *>);

		void init();

		bool trace(wa_ray Wi, wa_primitive *intr, float) const;
		float traceShadow(float3, float3, float) const;

		wa_camera *getCamera();
		std::vector<wa_light *> getLights() const;

		void getImageResolution(int *, int *) const;

		~wa_content();


	private:
		void loadLightGeoIntoEmbree(wa_light *);
		void loadGeoFileIntoEmbree(string);
		
		
		int raytracer_type;
		
		
		wa_camera camera;
		std::vector<wa_light *> lights;
		std::vector<wa_object *> objects;

		wa_raytracer raytracer;
};

class analyticSphereObject
{
	public:
		float3 P;
		float r;

		unsigned int geomID;
};