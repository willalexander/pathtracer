#pragma once

#include "objects.h"


class wa_primitive
{
	public:
		wa_primitive();
		wa_primitive(float3 Pval, float3 Nval);
		~wa_primitive();

		int getLightOrSurface() { return lightOrSurface; }
		void setLightOrSurface(int val) { lightOrSurface = val; }
	
		float getSurfaceArea() { return surfaceArea; }
		void setSurfaceArea(float val) { surfaceArea = val; }

		wa_colour getPower() { return power; }
		void setPower(wa_colour val) { power = val; }

		wa_colour getIrradiance() { return irradiance;  }
		void setIrradiance(wa_colour val) { irradiance = val;  }
		
		wa_colour Le(float3);
		wa_colour BRDF(float3, float3);

		wa_object *getObject();
		void setObject(wa_object *);

		float3 getP();
		float3 getN();

		void setP(float3 val);
		void setN(float3 val);

	private:
		int lightOrSurface;
		float surfaceArea;
		wa_colour power;
		wa_colour irradiance;

		wa_object *objectPointer;

		float3 P;
		float3 N;
};


