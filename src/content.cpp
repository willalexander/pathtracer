
#include "content.h"


wa_content::wa_content(wa_camera camera_in, std::vector<wa_light *> lights_in, std::vector<wa_object *> objects_in)
	: camera(camera_in),
	lights(lights_in),
	objects(objects_in)
{ 
	raytracer_type = RAYTRACER_EMBREE;
}

void wa_content::init()
{
	if (raytracer_type == RAYTRACER_EMBREE)
	{
		raytracer.init();

		/*Load lights into the raytracer:*/
		for (int l = 0; l < lights.size(); l++)
		{
			objects[l]->setID(l);
			objects[l]->shader.objID = l;
			loadLightGeoIntoEmbree(lights[l]);
		}

		/*Load geometry into the raytracer:*/
		for (int g = lights.size(); g < objects.size(); g++)
		{
			objects[g]->setID(g);
			objects[g]->shader.objID = g;
			loadGeoFileIntoEmbree(objects[g]->getFileName());
		}

		raytracer.commit();
	}

	scene = raytracer.getEmbreeScene();
	rtcInitIntersectContext(&context);
}

void wa_content::loadLightGeoIntoEmbree(wa_light *l)
{
	//An area light is represented by a single geo quad:
	if(l->type() == LIGHT_AREA_TYPE)
	{
		wa_areaLight *al = (wa_areaLight *)(l);
		raytracer.loadLightGeo(X(l->P), Y(l->P), Z(l->P), X(al->U), Y(al->U), Z(al->U), X(al->V), Y(al->V), Z(al->V));
	}
}


void wa_content::loadGeoFileIntoEmbree(string fileName)
{
	raytracer.loadGeo(fileName);
}

bool wa_content::trace(wa_ray Wi, wa_primitive *intr_out, float bias) const
{
		float3 hitP, hitN;
		int hitID;
		bool hit = raytracer.trace(X(Wi.P), Y(Wi.P), Z(Wi.P), X(Wi.D), Y(Wi.D), Z(Wi.D), &(X(hitP)), &(Y(hitP)), &(Z(hitP)), &(X(hitN)), &(Y(hitN)), &(Z(hitN)), &hitID, bias);


		if(hit)
		{
			normalize(hitN);

			wa_primitive newPrim(hitP, hitN);
			newPrim.setObject(objects[hitID]);

			/*Light:*/
			if (hitID < lights.size())
			{
				newPrim.setLightOrSurface(0);
				newPrim.setIrradiance(lights[hitID]->getIrradiance());
			}

			/*Geo:*/
			else
			{
				newPrim.setLightOrSurface(1);
			}

			*intr_out = newPrim;
			return true;
		}
		return false;
}

bool wa_content::traceFAST(wa_ray *Wi, wa_primitive *intr)
{
	RTCRayHit rh;
	rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	rh.ray = { X(Wi->P), Y(Wi->P), Z(Wi->P), traceBias, X(Wi->D), Y(Wi->D), Z(Wi->D), 0.0, 10000000.0, 0, 0, 0 };

	rtcIntersect1(scene, &context, &rh);
	X(intr->P) = rh.ray.org_x + rh.ray.tfar * rh.ray.dir_x;
	Y(intr->P) = rh.ray.org_y + rh.ray.tfar * rh.ray.dir_y;
	Z(intr->P) = rh.ray.org_z + rh.ray.tfar * rh.ray.dir_z;

	if (rh.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		X(intr->N) = rh.hit.Ng_x;
		Y(intr->N) = rh.hit.Ng_y;
		Z(intr->N) = rh.hit.Ng_z;

		normalize(intr->N);

		intr->objectPointer = objects[rh.hit.geomID];
		intr->objID = rh.hit.geomID;

		/*Light:*/
		if (rh.hit.geomID < lights.size())	intr->irradiance = lights[rh.hit.geomID]->getIrradiance();

		/*Geo:*/
		else intr->lightOrSurface = 1;

		return true;
	}

	return false;
}


wa_camera *wa_content::getCamera()
{
	return &camera;
}

std::vector<wa_light *> wa_content::getLights() const
{
	return lights;
}

void wa_content::getImageResolution(int *resI_out, int *resJ_out) const
{
	*resI_out = camera.getResI();
	*resJ_out = camera.getResJ();
}

void wa_content::setTraceBias(float val)
{
	traceBias = val;
}

wa_content::~wa_content()
{
}