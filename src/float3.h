#pragma once

#include <iostream>
using std::ostream;

#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <boost/qvm/vec_traits_defaults.hpp>
using namespace boost::qvm;


struct float3
{
	friend float3 normalized(float3);
	friend bool operator!=(const float3&, const float&);

	float3();
	float3(float);
	float3(float, float, float);


	float a[3];
};

#define XAX float3{1.0, 0.0, 0.0}
#define YAX float3{0.0, 1.0, 0.0}
#define ZAX float3{0.0, 0.0, 1.0}


template <>
struct vec_traits<float3> : vec_traits_defaults<float3, float, 3>
{
	template <int I> static inline scalar_type & write_element(float3 & v) { return v.a[I]; }

	static inline scalar_type & write_element_idx(int i, float3 & v) { return v.a[i]; } //optional
};



ostream& operator<<(ostream&, const float3&);
void operator*=(float3&, float);