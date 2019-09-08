#pragma once

#include <iostream>
#include <vector>

#include "platform.h"
#include "float3.h"
#include "randomnumbers.h"

using std::ostream;

#define M_PI 3.1415926535


#define XAX float3{1.0, 0.0, 0.0}
#define YAX float3{0.0, 1.0, 0.0}
#define ZAX float3{0.0, 0.0, 1.0}


struct wa_vertex { float x, y, z, a; };
struct wa_triangle { int v0, v1, v2; };
struct wa_quad { int v0, v1, v2, v3; };

struct wa_ray
{
	float3 P, D;
};


class wa_colour
{
	friend ostream& operator<<(ostream&, const wa_colour&);
	friend wa_colour operator*(const float&, const wa_colour&);
	friend wa_colour operator*(const wa_colour&, const float&);
	friend wa_colour operator/(const wa_colour&, const float&);
	friend bool operator==(const wa_colour&, const float&);
	friend bool operator<(const wa_colour&, const float&);
	friend wa_colour exp(const wa_colour&);

	public:
		wa_colour();
		wa_colour(float val);
		wa_colour(float valR, float valG, float valB);

		wa_colour operator*(const float&);
		wa_colour operator*(const wa_colour&);
		wa_colour operator/(const float&);
		void operator+=(const wa_colour&);
		bool operator!=(const float&);
		void operator*=(const wa_colour&);
		void operator/=(const float&);
		wa_colour operator/(const wa_colour&);
		wa_colour operator+(const wa_colour&);
		wa_colour operator-(const wa_colour&);

		void set(float, float, float);
		void setC(int, float);

		void placeInFrameBuffer(float *, int, int, int);

		float average();

		void clamp(float, float);

		float c(int);

		float r, g, b;
};


class wa_camera
{
	friend ostream& operator<<(ostream&, const wa_camera&);

	public:
		wa_camera();
		wa_camera(int, float3, float3, float3, float, int, int, float);

		int getResI() const;
		int getResJ() const;
		wa_ray generateRayForPixel(int, int) const;
		float3 getP();
		float3 getD();
		void setP(float3);
		void setAt(float3);
		void setResolution(int, int);
		void getSensorInfo(float3 *, float3 *, float3 *, float3 *);
		
		float3 P;

		void rotateY(float);
		float3 rotateVecY(float3, float);


	private:
		int type;
		
		float3 D;
		float3 at;
		float3 up_simple;
		float fov;
		int resI, resJ;

		float sensorWidth;

		float3 up;
		float3 side;
		float sensorDist;
		float3 sensorOrigin;
		float pixelWidth;
};




class wa_path
{
	public:
		wa_path() {};

		int getSize();
		void addVert(float3);
		float3 getVert(int);
		
	private:
		std::vector<float3> verts;
};

class wa_line
{
	public:
		wa_line(float3 Aval, float3 Bval) : A(Aval), B(Bval) {};

		float3 getA() const { return A; };
		float3 getB() const { return B; };

	private:

		const float3 A;
		const float3 B;
};

struct wa_sphere
{
	float3 P;
	float r;
};

float randFloat();
bool randBool();

void generateRandTable();
float randFloat_table();

float3 hemisphereSample(float3, float *);
float3 hemisphereSampleUniform(float3, float *);
float3 hemisphereSampleUniform(float3, float, float *);
float3 sphereSampleUniform();
float3 sphereSampleBiased(float3, float, float *);
float3 sphereSampleBiased2(float3, float, float *);
float3 sphereSampleBiased3(float3, float, float *);
bool solveQuadraticEquation(float, float, float, float *, float *);

void generateSpherePoints(float3, float, int, int, int, int, wa_vertex *, wa_vertex *, wa_vertex *, wa_vertex *);

float surfaceArea_to_solidAngle(float3, float3, float3);




