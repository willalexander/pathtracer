#pragma once

#define WINDOWS

struct wa_lightSimple
{
	float px, py, pz;
	float ux, uy, uz;
	float vx, vy, vz;

	bool dirty;
};