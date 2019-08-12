#pragma once

#include <iostream>
#include <string>

#include "fundamentals.h"

using namespace std;

#define SHADER_TYPE_DIFFUSE 0 
#define SHADER_TYPE_SSS 1

struct wa_shader
{
	wa_shader();
	float fresnel(float3, float3);

	wa_colour kd;
	wa_colour kt;
	wa_colour ks;

	wa_colour sa;
	wa_colour se;	
	wa_colour ss;
	wa_colour st;

	int proceduralVolume = 1;
	bool homogenous = true;

	bool isOpaque = true;
	bool isClear = true;
	bool isLight = false;
	bool diffuse = true;
	bool specular = false;
	bool dielectric = false;

	void initialize();
	float BRDF(float3, float3);
	float phase(float);

	wa_colour st_f(float3);

	wa_colour minimumExtinctionCoefficient();
	wa_colour maximumExtinctionCoefficient();
};

class wa_object
{	
	public:
		wa_object(string, string, wa_shader);

		int getID();
		void setID(int);

		string getFileName();
		int getShaderType();

		wa_shader getShader();
		void setShader(wa_shader);

	private:
		int ID;
		string fileName;
		int shaderType;

		wa_shader shader;
};