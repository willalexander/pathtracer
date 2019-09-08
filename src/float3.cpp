#include "float3.h"

float3 normalized(float3 A)
{
	float magnitude = mag(A);
	float3 out = { A.a[0] / magnitude, A.a[1] / magnitude, A.a[2] / magnitude };

	return out;
}

bool operator==(const float3& A, const float& B)
{
	if ((X(A) == B) && (Y(A) == B) && (Z(A) == B)) return true;
	return false;
}

bool operator!=(const float3& A, const float& B)
{
	if ((X(A) != B) || (Y(A) != B) || (Z(A) != B)) return true;
	return false;
}

float3::float3()
{
	a[0] = a[1] = a[2] = 0;
}

float3::float3(float x)
{
	a[0] = a[1] = a[2] = x;
}

float3::float3(float x, float y, float z)
{
	a[0] = x;
	a[1] = y;
	a[2] = z;
}

ostream& operator<<(ostream& output, const float3& val)
{
	output << "Vector " << val.a[0] << " " << val.a[1] << " " << val.a[2];
	return output;
}

void operator*=(float3 &val, float a)
{
	X(val) = a * X(val);
	Y(val) = a * Y(val);
	Z(val) = a * Z(val);
}
