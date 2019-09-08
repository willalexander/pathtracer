#include <iostream>

#include "fundamentals.h"

using namespace std;


float randTable[100];
int randCount = 0;

void generateRandTable()
{
	for (int i = 0; i < 100; i++)
	{
		randTable[i] = randFloat();
	}
}

float randFloat_table()
{
	if(randCount == 100) randCount = 0;
	return randTable[randCount++];
}


ostream &operator<<(ostream &output, const wa_colour &val)
{
	output << "wa_colour " << val.r << ", " << val.g << ", " << val.b;
	return output;
}

wa_colour operator*(const float& A, const wa_colour& B)
{
	wa_colour out(A * B.r, A * B.g, A * B.b);
	return out;
}

wa_colour operator*(const wa_colour& A, const float& B)
{
	wa_colour out(A.r * B, A.g * B, A.b * B);
	return out;
}

wa_colour operator/(const wa_colour& A, const float& B)
{
	wa_colour out(A.r / B, A.g / B, A.b / B);
	return out;
}

bool operator==(const wa_colour& A, const float& B)
{
	if ((A.r == B)&&(A.g == B)&&(A.b == B)) return true;
	return false;
}

wa_colour exp(const wa_colour& A)
{
	wa_colour out(exp(A.r), exp(A.g), exp(A.b));
	return out;
}

wa_colour::wa_colour()
{
	r = g = b = 0;
}

wa_colour::wa_colour(float val)
{
	r = g = b = val;
}

wa_colour::wa_colour(float valR, float valG, float valB)
{
	r = valR;
	g = valG;
	b = valB;
}

void wa_colour::placeInFrameBuffer(float *frameBuffer, int bufferWidth, int i, int j)
{
	frameBuffer[(j * bufferWidth + i) * 3 + 0] = r;
	frameBuffer[(j * bufferWidth + i) * 3 + 1] = g;
	frameBuffer[(j * bufferWidth + i) * 3 + 2] = b;
}

void wa_colour::set(float valR, float valG, float valB)
{
	r = valR;
	g = valG;
	b = valB;
}

wa_colour wa_colour::operator*(const float& B)
{
	wa_colour out(r * B, g * B, b * B);
	return out;
}

wa_colour wa_colour::operator*(const wa_colour& B)
{
	wa_colour out(r * B.r, g * B.g, b * B.b);
	return out;
}

void wa_colour::operator+=(const wa_colour& B)
{
	r += B.r;
	g += B.g;
	b += B.b;
}

wa_colour wa_colour::operator/(const float& B)
{
	wa_colour out(r / B, g / B, b / B);
	return out;
}

bool wa_colour::operator!=(const float& B)
{
	if ((r == B) && (g == B) && (b == B)) return false;
	return true;
}

void wa_colour::operator*=(const wa_colour& B)
{
	r *= B.r;
	g *= B.g;
	b *= B.b;
}

void wa_colour::operator/=(const float& B)
{
	r /= B;
	g /= B;
	b /= B;
}

wa_colour wa_colour::operator/(const wa_colour& A)
{
	wa_colour out(r / A.r, g / A.g, b / A.b);
	return out;
}

wa_colour wa_colour::operator+(const wa_colour& A)
{
	wa_colour out(r + A.r, g + A.g, b + A.b);
	return out;
}

wa_colour wa_colour::operator-(const wa_colour& A)
{
	wa_colour out(r - A.r, g - A.g, b - A.b);
	return out;
}

bool operator<(const wa_colour& A, const float& B)
{
	if ((A.r < B) && (A.g < B) && (A.b < B)) return true;
	return false;
}

float wa_colour::average()
{
	return (r + g + b) / 3.0;
}

void wa_colour::clamp(float min, float max)
{
	for(int co = 0; co < 3; co++)
	{
		if(c(co) < min) setC(co, min);
		if(c(co) > max) setC(co, max);
	}
}

/* Returns a single compoenent of the colour, indexed by an integer 0, 1, or 2:*/
float wa_colour::c(int component)
{
	if (component == 0) return r;
	if (component == 1) return g;
	if (component == 2) return b;
}

/* Sets a particular component of the colour*/
void wa_colour::setC(int component, float val)
{
	if(component == 0) r = val;
	if(component == 1) g = val;
	if(component == 2) b = val;
}



wa_camera::wa_camera()
	: P({0.0, 0.0, 0.0}),
	D({0.0, 0.0, 1.0}),
	up_simple({0.0, 1.0, 0.0}),
	fov(1.0)
{
	resI = 1024;
	resJ = 1024;
	sensorWidth = 1.0;


	type = 0;

	side = cross(D, up_simple);
	normalize(side);
	up = cross(side, D);
	//sensorDist = 0.5 * sensorWidth * ((float)(resJ) / (float)(resI)) / tan(0.5 * fov * (M_PI / 180.0));

	if (type == 1) sensorWidth = 1.0;

	cout << "sensorWidth: " << sensorWidth << endl;
	cout << "fov: " << fov << endl;
	cout << "tan: " << tan(0.5 * fov * (M_PI / 180.0)) << endl;

	sensorDist = 0.5 * sensorWidth / tan(0.5 * fov * (M_PI / 180.0));

	setResolution(1024, 1024);
}


wa_camera::wa_camera(int typeVal, float3 PVal, float3 DVal, float3 upVal, float fovVal, int resIVal, int resJVal, float sensorWidthVal) 
	: P(PVal),
	D(DVal),
	up_simple(upVal),
	fov(fovVal)
{
	resI = resIVal;
	resJ = resJVal;
	sensorWidth = sensorWidthVal;


	type = typeVal;

	side = cross(D, up_simple);
	normalize(side);
	up = cross(side, D);
	//sensorDist = 0.5 * sensorWidth * ((float)(resJ) / (float)(resI)) / tan(0.5 * fov * (M_PI / 180.0));

	if (type == 1) sensorWidth = sensorWidthVal;

	cout << "sensorWidth: " << sensorWidth << endl;
	cout << "fov: " << fov << endl;
	cout << "tan: " << tan(0.5 * fov * (M_PI / 180.0)) << endl;

	sensorDist = 0.5 * sensorWidth / tan(0.5 * fov * (M_PI / 180.0));

	setResolution(resIVal, resJVal);
}

void wa_camera::setResolution(int resIVal, int resJVal)
{
	resI = resIVal;
	resJ = resJVal;

	sensorOrigin = P + sensorDist * D - 0.5 * sensorWidth * (side + (((float)(resJ) / (float)(resI)) * up));
	pixelWidth = sensorWidth / (float)(resI);
}



ostream& operator<<(ostream& output, const wa_camera& B)
{
	output << "Camera:" << endl;
	output << "Type: " << B.type << endl;
	output << "P: " << B.P << endl;
	output << "D: " << B.D << endl;
	output << "up: " << B.up << endl;
	output << "fov: " << B.fov << endl;
	output << "resI, resJ: " << B.resI << ", " << B.resJ << endl;

	output << "side: " << B.side << endl;
	output << "sensorWidth: " << B.sensorWidth << endl;
	output << "sensorDist: " << B.sensorDist << endl;
	output << "sensorOrigin: " << B.sensorOrigin << endl;
	output << "pixelWidth: " << B.pixelWidth << endl;

	return output;
}

int wa_camera::getResI() const
{
	return resI;
}

int wa_camera::getResJ() const
{
	return resJ;
}

float3 wa_camera::getP()
{
	return P;
}

float3 wa_camera::getD()
{
	return D;
}

wa_ray wa_camera::generateRayForPixel(int i, int j) const
{
	wa_ray result;

	//Perspective camera:
	if (type == 0)
	{
		float3 D = sensorOrigin + pixelWidth * ((i + randFloat()) * side + (j + randFloat()) * up) - P;
		normalize(D);

		result = { P, D };
	}

	//Orthographic camera:
	if (type == 1)
	{
		result.P = P + (pixelWidth * (i + randFloat()) - 0.5 * sensorWidth) * side + (pixelWidth * (j + randFloat()) - 0.5 * sensorWidth) * up;
		result.D = D;
	}
	
	return result;
}

void wa_camera::getSensorInfo(float3 *camOrigin, float3 *sensorOrigin_out, float3 *iVec, float3 *jVec)
{
	*camOrigin = P;
	*sensorOrigin_out = sensorOrigin;

	*iVec = pixelWidth * side;
	*jVec = pixelWidth * up;
}

void wa_camera::setP(float3 P_in)
{
	P = P_in;
}

void wa_camera::setAt(float3 at_in)
{
	at = at_in;
}

int wa_path::getSize()
{
	return verts.size();
}

void wa_path::addVert(float3 v)
{
	verts.push_back(v);
}

float3 wa_path::getVert(int i)
{
	return verts[i];
}


void wa_camera::rotateY(float angle)
{
	//Rotate about the 'at' position by angle 'angle':
	float3 Prel = P - at;
	float3 Prel_new = rotateVecY(Prel, angle);
	P = at + Prel_new;

	D = rotateVecY(D, angle);
	up = rotateVecY(up, angle);
	side = rotateVecY(side, angle);

	sensorOrigin = P + sensorDist * D - 0.5 * sensorWidth * (side + (((float)(resJ) / (float)(resI)) * up));
}

float3 wa_camera::rotateVecY(float3 A, float angle)
{
	float3 B;

	X(B) = cos(angle) * X(A) + sin(angle) * Z(A);
	Y(B) = Y(A);
	Z(B) = -1.0 * sin(angle) * X(A) + cos(angle) * Z(A);

	return B;
}


float3 hemisphereSample(float3 D, float *weight)
{
	float3 out;

	float3 tangentA = cross(D, XAX);
	if ((tangentA.a[0] == 0.0) && (tangentA.a[1] == 0.0) && (tangentA.a[2] == 0.0)) tangentA = cross(D, YAX);
	float3 tangentB = cross(tangentA, D);

	float theta = randFloat() * 0.5 * M_PI;
	float phi = randFloat() * 2.0 * M_PI;

	out = cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);

	if(weight != NULL) *weight = sin(theta) * M_PI / 2.0;

	return out;
}

float3 hemisphereSampleUniform(float3 D, float *reciprocal_pdf)
{
	float3 out;

	float3 tangentA = normalized(cross(D, XAX));
	if((Y(D) == 0.0)&&(Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	float3 tangentB = normalized(cross(tangentA, D));

	float cosTheta = 1.0 - randFloat();
	float theta = acos(cosTheta);
	float phi = randFloat() * 2.0 * M_PI;

	//out = cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
	out = cosTheta * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
	if (reciprocal_pdf != NULL) *reciprocal_pdf = 2.0 * M_PI;

	return out;
}


/*Sample a partial hemisphere, from the central axis up to 'maxAngle':*/
float3 hemisphereSampleUniform(float3 D, float maxAngle, float *weight)
{
	float3 out;

	float3 tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	float3 tangentB = normalized(cross(tangentA, D));

	float theta = acos(1.0 - randFloat() * (1.0 - cos(maxAngle)));
	float phi = randFloat() * 2.0 * M_PI;

	out = cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
	if (weight != NULL) *weight = 1;

	return out;
}

//
// sphereSampleBiased() -	Generates a random sample in the unit sphere, biased towards the direction D.
//							The higher 'a', the more biased the sample
//
float3 sphereSampleBiased(float3 D, float a, float *reciprocal_pdf)
{
	float3 out;
	float3 tangentA, tangentB;
	float phi, theta;

	//Generate tangent vectors based on D:
	tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	tangentB = normalized(cross(tangentA, D));

	//phi values are distributed uniformly:
	phi = randFloat() * 2.0 * M_PI;

	//theta values are distibuted with a probability proportional to cos(theta/2):
	theta = 2.0 * acos(pow(1.0 - randFloat(), 1.0/3.0));


	//compute the PDF:
	*reciprocal_pdf = 1.0 / ( (3.0/(8.0*M_PI)) * cos(0.5*theta) );

	//std::cout << *reciprocal_pdf << std::endl;

	//Compute the direction vector:
	return cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
}

float3 sphereSampleBiased2(float3 D, float a, float *reciprocal_pdf)
{
	float3 out;
	float3 tangentA, tangentB;
	float phi, theta;

	//Generate tangent vectors based on D:
	tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	tangentB = normalized(cross(tangentA, D));

	//phi values are distributed uniformly:
	phi = randFloat() * 2.0 * M_PI;

	//theta values are distibuted uniformly:
	theta = 0.5 * M_PI * randFloat();

	//compute the PDF:
	*reciprocal_pdf = M_PI*M_PI * sin(theta);

	//Compute the direction vector:
	return cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
}

float3 sphereSampleBiased3(float3 D, float a, float *reciprocal_pdf)
{
	float3 out;
	float3 tangentA, tangentB;
	float phi, theta;

	//Generate tangent vectors based on D:
	tangentA = normalized(cross(D, XAX));
	if ((Y(D) == 0.0) && (Z(D) == 0.0)) tangentA = normalized(cross(D, YAX));
	tangentB = normalized(cross(tangentA, D));

	//phi values are distributed uniformly:
	phi = randFloat() * 2.0 * M_PI;

	//theta values are distibuted according to the pattern (PI/2 - theta) ^ a:
	theta = 0.5*M_PI - pow(pow(0.5*M_PI, a+1)*randFloat(), 1.0 / (a + 1.0));

	//compute the PDF:
	*reciprocal_pdf = 1.0 / (((a+1)/(4.0*pow(0.5*M_PI,a+2))) * pow(0.5*M_PI - theta, a) / sin(theta));

	//Compute the direction vector:
	return cos(theta) * D + sin(theta) * (cos(phi) * tangentA + sin(phi) * tangentB);
}

float3 sphereSampleUniform()
{
	float3 out;

	float theta = acos(1.0 - (2.0 * randFloat()));
	float phi = randFloat() * 2.0 * M_PI;

	out = cos(theta) * YAX + sin(theta) * (cos(phi) * XAX + sin(phi) * ZAX);

	return out;
}


bool solveQuadraticEquation(float a, float b, float c, float *s1, float *s2)
{
	float b2_4ac = pow(b, 2) - (4 * a * c);
	if (b2_4ac < 0) return false;

	if(s1 != NULL) *s1 = ((-1 * b) + sqrt(b2_4ac)) / (2 * a);
	if(s2 != NULL) *s2 = ((-1 * b) - sqrt(b2_4ac)) / (2 * a);

	return true;
}

void generateSpherePoints(float3 P, float rad, int latSeg, int lonSeg, int i, int j, wa_vertex *v0, wa_vertex *v1, wa_vertex *v2, wa_vertex *v3)
{
	float latDist = M_PI / (float)(latSeg);
	float lonDist = 2.0 * M_PI / (float)(lonSeg);

	float theta;
	float phi;
	
	float3 res[4];
	float3 D;

	for (int a = 0; a < 2; a++)
	{
		for (int b = 0; b < 2; b++)
		{
			theta = (i + a) * latDist;
			phi = (j + b) * lonDist;

			D = cos(theta) * YAX + sin(theta) * (cos(phi) * XAX + sin(phi) * ZAX);
			res[a * 2 + b] = P + rad * D;
		}
	}

	v0->x = X(res[0]); v0->y = Y(res[0]); v0->z = Z(res[0]);
	v1->x = X(res[1]); v1->y = Y(res[1]); v1->z = Z(res[1]);
	v2->x = X(res[3]); v2->y = Y(res[3]); v2->z = Z(res[3]);
	v3->x = X(res[2]); v3->y = Y(res[2]); v3->z = Z(res[2]);
}


//Convert from surface area at B with normal N, to solid angle subtended at A:
float surfaceArea_to_solidAngle(float3 A, float3 B, float3 N)
{
	float3 AtoB = B - A;

	return fabs(dot(normalized(AtoB), N)) / pow(mag(AtoB), 2);
}