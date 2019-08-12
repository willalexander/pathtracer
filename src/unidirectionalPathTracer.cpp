
#include <iostream>
#include "stdlib.h"

#include "unidirectionalPathTracer.h"

using namespace std;

#define DEBUG_I 275
#define DEBUG_J 400
#define DEBUG_WIDTH -1

#define TRACE_BIAS 0.001
#define TRACE_DEPTH 3
#define RAYMARCH_STEPS_PER_METRE 0.02
#define LOW_VISIBILITY_THRESHOLD 0.01

#define NEXT_EVENT_ESTIMATION true

#define VERBOSE



udp_renderer::udp_renderer(wa_content *content_in, float *frameBufferPtr, int numPixelSamplesVal, int numHemisphereSamplesVal, int numSphereSamplesVal, int numLightSamplesVal)
	: theContent(content_in),
	lights(content_in->getLights()),
	numPixelSamples(numPixelSamplesVal),
	numHemisphereSamples(numHemisphereSamplesVal),
	reciprocal_numHemisphereSamples(1.0 / numHemisphereSamplesVal),
	numSphereSamples(numSphereSamplesVal),
	reciprocal_numSphereSamples(1.0 / numSphereSamplesVal),
	numLightSamples(numLightSamplesVal),
	reciprocal_numLightSamples(1.0 / numLightSamplesVal),
	worldBoundSphere_P(275.0, 275.0, 275.0),
	worldBoundSphere_rad(1100.0)
{
	frameBuffer = frameBufferPtr;

	camera = theContent->getCamera();

	disW = camera->getResI();
	disH = camera->getResJ();

	//By default, the render resolution is the same as the display resolution:
	renW = disW;
	renH = disH;
	rf = 1;
	setRegion(0, renW, 0, renH);

	//Global shader describes the default volume light transport properties of the whole world. Set to 'clear' by default!
	/*globalShader.isClear = false;
	globalShader.ss.set(0.001, 0.001, 0.001);
	globalShader.st.set(0.001, 0.001, 0.001);*/

	numThreads = 4;

	//Create the internal buffer and set it to black:
	internalBuf = (float *)(malloc(disW * disH * 3 * sizeof(float)));
	reset();
}

void udp_renderer::setRenderResolution(int downresFactorOf2)
{
	lowerResolution(downresFactorOf2, &renW, &renH);
	rf = (int)(pow(2.0, downresFactorOf2));

	camera->setResolution(renW, renH);

	setRegion(0, renW, 0, renH);
}

void udp_renderer::render(float *buffer)
{
	frameBuffer = buffer;

	wa_ray sampleRay;
	wa_colour pixelValue;

	//Define the boundaries of each render tile, one for each thread:
	divideIntoTiles();

	abortRender = false;

	//Render on threads: (Boost on Windows, pthread on Linux)
	progress = 0.0;

	boost::thread_group threads;

	for (int t = 0; t < numThreads; t++)
	{
		const auto tileFunc = [&, t]() { renderTile(t); };
		threads.create_thread(tileFunc);
	}
	
	threads.join_all();
}

void udp_renderer::renderTile(int t)
{
	wa_ray sampleRay;
	wa_colour pixelValue;

	int latestReportedProgress = 0;

	for (int j = tiles[t].b; j < tiles[t].t; j++)
	{
		for (int i = tiles[t].l; i < tiles[t].r; i++)
		{
			bool debugPix = false;
			if ((i >= (DEBUG_I - DEBUG_WIDTH)) && (i <= (DEBUG_I + DEBUG_WIDTH)) && (j >= (DEBUG_J - DEBUG_WIDTH)) && (j <= (DEBUG_J + DEBUG_WIDTH))) debugPix = true;

			pixelValue = 0.0;

			//Shoot a ray for each pixel sample:
			for (int s = 0; s < numPixelSamples; s++)
			{
				sampleRay = camera->generateRayForPixel(i, j);

				wa_path samplePath;
				if (debugPix) samplePath.addVert(sampleRay.P);
				bool hit;
				
				wa_colour sampleValue = radiance(sampleRay.P, -1.0 * sampleRay.D, 0, true, (debugPix) ? &samplePath : NULL, CAMERA_RAY) * (1.0 / (float)(numPixelSamples));
				pixelValue += sampleValue;

				if (debugPix) recordedPaths.push_back(samplePath);
				if (debugPix) cout << "Pix val: " << pixelValue << endl;
				if (debugPix) pixelValue = wa_colour(1.0, 0.0, 0.0);
			}

			float itCntF = (float)(itCnt);

			internalBuf[(j * disW + i) * 3 + 0] = (itCntF / (itCntF + 1))*internalBuf[(j * disW + i) * 3 + 0] + (1.0 / (itCntF + 1))*pixelValue.r;
			internalBuf[(j * disW + i) * 3 + 1] = (itCntF / (itCntF + 1))*internalBuf[(j * disW + i) * 3 + 1] + (1.0 / (itCntF + 1))*pixelValue.g;
			internalBuf[(j * disW + i) * 3 + 2] = (itCntF / (itCntF + 1))*internalBuf[(j * disW + i) * 3 + 2] + (1.0 / (itCntF + 1))*pixelValue.b;

			float tmp[3] = { internalBuf[(j * disW + i) * 3 + 0], internalBuf[(j * disW + i) * 3 + 1], internalBuf[(j * disW + i) * 3 + 2] };
		
			pixelValue.set(tmp[0], tmp[1], tmp[2]);

			//Put the pixel value in the framebuffer:
			for (int i1 = i * rf; i1 < (i + 1) * rf; i1++)
			{
				for (int j1 = j * rf; j1 < (j + 1) * rf; j1++)
				{
					if ((i1 >= disW)||(j1 >= disH)) continue;
					pixelValue.placeInFrameBuffer(frameBuffer, disW, i1, j1);
				}
			}

			pixelValue.placeInFrameBuffer(frameBuffer, disW, i, j);
		}

		progress += (1.0 / ((float)(numThreads) * (float)(tiles[t].t - tiles[t].b)));

		//The 'master' thread reports on the progress of all theads:
		if (t == 0)
		{
			int intProgress = 5 * (int)(progress * 20.0);
			if (intProgress > latestReportedProgress)
			{
				cout << "Progress: " << intProgress << "%" << endl;
				latestReportedProgress = intProgress;
			}
		}
	}
}

void udp_renderer::divideIntoTiles()
{
	tiles = (tile *)(malloc(numThreads * sizeof(tile)));
	int dim = (int)(sqrt((float)(numThreads)));

	for (int j = 0; j < dim; j++)
	{
		for (int i = 0; i < dim; i++)
		{
			tiles[j * dim + i].l = rr.l + (int)(i * ((float)(rr.w) / (float)(dim)));
			tiles[j * dim + i].r = rr.l + (int)((i + 1) * ((float)(rr.w) / (float)(dim)));
			tiles[j * dim + i].b = rr.b + (int)(j * ((float)(rr.h) / (float)(dim)));
			tiles[j * dim + i].t = rr.b + (int)((j + 1) * ((float)(rr.h) / (float)(dim)));
		}
	}
}

/*
	radiance()						Returns the radiance at point P0 in direction D:
*/
wa_colour udp_renderer::radiance(float3 P0, float3 D, int traceDepth, bool includeEmission, wa_path *debugPath, int rayType)
{
	float3 _D = -1.0 * D;			// The direction vector in the opposite direction to D
	wa_ray r1;						// wa_ray object for ray tracing
	wa_primitive intr;				// The intersection object for the first surface hit when tracing from P in direction -D						
	float3 P1;						// Position vector of the first surface intersected when ray tracing from P0
	float3 N1g, N1f;				// The geometric and face-forward normal vectors of the first surface intersected
	wa_object *object;				// The first surface intersected
	wa_shader shader;				// The shader of the first surface intersected
	wa_shader medium;				// The shader for the medium that P0 & P1 are in
	bool hit;						// Whether or not we hit a surface as a result of trace()
	float lambda;					// Solution to quadratic equation for ray intersection with the outer world bound
	wa_colour Lo(0.0);				// Result


	//By default, the current medium is the global atmosphere:
	medium = globalShader;

	//Trace a ray in the opposite direction of D to find the first visible surface:
	r1.P = P0;
	r1.D = _D;

	hit = theContent->trace(r1, &intr, TRACE_BIAS);

	//If nothing was hit, then limit the ray to the outer bound of the world, consider any volumetric scattering, and exit:
	if(!hit)
	{
		//If the ray origin is outside of the bound sphere, then stop and return black:
		if(mag(P0 - worldBoundSphere_P) >= worldBoundSphere_rad) return 0.0;

		solveQuadraticEquation(1.0, 2.0 * dot(P0 - worldBoundSphere_P, _D), pow(mag(P0 - worldBoundSphere_P), 2) - pow(worldBoundSphere_rad, 2), &lambda, NULL);
		P1 = P0 + lambda * _D;

		return volumeScattering(medium, P0, P1, traceDepth, debugPath);
	}

	P1 = intr.getP();
	N1g = intr.getN();
	object = intr.getObject();
	shader = object->getShader();

	if(debugPath) debugPath->addVert(P1);

	//If this intersection event is entering the surface:
	if (dot(N1g, D) > 0.0)
	{
		N1f = N1g;
		medium = globalShader;
	}

	//If this intersection event is exiting the surface:
	else
	{
		N1f = -1.0 * N1g;
		medium = shader;
	}

	//Emission:
	if (includeEmission) Lo += intr.Le(D);

	//Transmission from the other side of the surface:
	if ((traceDepth + 1) <= TRACE_DEPTH) Lo += transmission(object, shader, P1, N1f, D, traceDepth + 1, debugPath);

	//Reflection:
	Lo += reflection(object, shader, P1, N1f, D, traceDepth, debugPath, rayType);

	//Attenuation along the beam from P1 to P0:
	Lo *= volumeAttenuation(medium, P0, P1, debugPath);

	//Volumetric scattering along the beam from P1 to P0:
	Lo += volumeScattering(medium, P0, P1, traceDepth, debugPath);

	return Lo;
}


wa_colour udp_renderer::transmission(wa_object *object, wa_shader shader, float3 P, float3 N, float3 D, int traceDepth, wa_path *debugPath)
{
	//If this is not a transmissive surface, then return nothing:
	if (shader.kt == 0.0) return 0.0;

	//Light passes straight transmissive surfaces for now:
	return radiance(P, D, traceDepth, false, debugPath, CAMERA_RAY);
}

wa_colour udp_renderer::reflection(wa_object *object, wa_shader shader, float3 P, float3 N, float3 D, int traceDepth, wa_path *debugPath, int rayType)
{
	wa_colour Lo_diffuse(0.0);									// Outgoing radiance due to diffuse reflection
	wa_colour Lo_specular(0.0);									// Outgoing radiance due to specular reflection
	float fresnel = shader.fresnel(N, D);						// The ratio between diffuse & specular reflection due to Fresnel reflection

	//Diffuse:
	if ((shader.diffuse)||(shader.dielectric))
	{
		Lo_diffuse = directReflection(object, shader, P, N, D, debugPath);
		
		if((traceDepth + 1) <= TRACE_DEPTH) Lo_diffuse += indirectReflection(object, shader, P, N, D, traceDepth + 1, debugPath);
	}

	//Specular:
	if(((shader.specular)||(shader.dielectric))&&(rayType != DIFFUSE_RAY))
	{
		if ((traceDepth + 1) <= TRACE_DEPTH) Lo_specular += shader.ks * radiance(P, -1.0 * reflectVector(D, N), traceDepth + 1, true, debugPath, SPECULAR_RAY);
	}

	if(shader.diffuse) return Lo_diffuse;
	if(shader.specular) return Lo_specular;
	if(shader.dielectric) return ((1.0 - fresnel) * Lo_diffuse + fresnel * Lo_specular);
}

wa_colour udp_renderer::directReflection(wa_object *object, wa_shader shader, float3 P, float3 N, float3 D, wa_path *debugPath)
{
	wa_colour Li;					// The incoming radiance from a sample of a light
	float3 LP;						// The position vector or a sample on a light surface
	float reciprocal_Lpdf;			// The reciprocal of the pdf of sampling a light source
	wa_colour Li_visibility;		// The visibility between P and a sample point on the surface of a light source
	float3 Wi;						// The incoming light direction (a unit direction vector from P into space)
	wa_colour Lo(0.0);				// The result


	//Diffuse refection. Compute the contribution from each light source:
	for(int l = 0; l < lights.size(); l++)
	{
		for(int i = 0; i < numLightSamples; i++)
		{
			if (lights[l]->radianceSample(P, N, &Li, &LP, &reciprocal_Lpdf, TRACE_BIAS))
			{
				Li_visibility = visibility(P, LP, debugPath);

				//This light sample contributes nothing if there is too much geometry between it and this surface point
				if (Li_visibility == 0.0) continue;

				Wi = normalized(P - LP);

				Lo += reciprocal_numLightSamples * shader.kd * shader.BRDF(Wi, D) * fabs(dot(N, Wi)) * Li_visibility * reciprocal_Lpdf * Li;
			}
		}
	}

	return Lo;
}

wa_colour udp_renderer::indirectReflection(wa_object *object, wa_shader shader, float3 P, float3 N, float3 D, int traceDepth, wa_path *debugPath)
{
	float3 Wi;						// The unit direction vector in the hemisphere above the surface at P
	float reciprocal_Hpdf;			// Reciprocal of the pdf for sampling the hemisphere
	wa_colour Li;					// The incoming radiance towards P from a particular direction
	wa_colour Lo(0.0);				// The result

	for (int i = 0; i < numHemisphereSamples; i++)
	{
		Wi = -1.0 * hemisphereSampleUniform(N, &reciprocal_Hpdf);

		//Compute the incoming radiance from this direction by tracing a ray 
		Li = radiance(P, Wi, traceDepth, false, debugPath, DIFFUSE_RAY);

		//Compute the outgoing radiance due to this incoming radiance:
		Lo += reciprocal_numHemisphereSamples * shader.kd * shader.BRDF(Wi, D) * fabs(dot(N, Wi)) * reciprocal_Hpdf * Li;
	}

	return Lo;
}

wa_colour udp_renderer::volumeAttenuation(wa_shader medium, float3 A, float3 B, wa_path *debugPath)
{
	wa_colour transmittance;				// Transmittance of this medium between A and B:
	float3 IB = B - A;						// The vector of the internal beam (from A -> B)
	float IB_l = mag(IB);					// The distance from A to B
	float3 S;								// The position vector of a scatter event along the internal beam 
	int RMS;								// The number of ray march steps that will be taken along the internal beam
	float reciprocal_RMS;					// The reciprocal of the number of ray march steps
	wa_samplerBeam *beamSampler;			// Object for generating samples along the beam
	float RMD;								// The distance of each scatter event from A
	float reciprocal_RMLpdf;				// The reciprocal of the pdf for generating scatter points along the internal beam


	//No attenuation if the medium is clear:
	if (medium.isClear) return 1.0;

	//Compute the optical thickness of this medium between A and B:
	transmittance = exp(-1.0 * medium.st * mag(B - A));

	//Beer's law:
	return transmittance;
}

wa_colour udp_renderer::volumeScattering(wa_shader medium, float3 A, float3 B, int traceDepth, wa_path *debugPath)
{
	float reciprocal_Lpdf;					// The reciprocal of the pdf of light sampling **measured relative to soild angle at the point receiving light**
	float3 IB = B - A;						// The vector of the internal beam (from A -> B)
	float IB_l = mag(IB);					// The distance from A to B
	float3 S;								// The position vector of a scatter event along the internal beam 
	int RMS;								// The number of ray march steps that will be taken along the internal beam
	float reciprocal_RMS;					// The reciprocal of the number of ray march steps
	float RMD;								// The distance of each scatter event from A
	float reciprocal_RMLpdf;				// The reciprocal of the pdf for generating scatter points along the internal beam

	wa_colour S_Lo;							// The outgoing radiance at a scatter point on the internal beam
	wa_colour S_Lo_0;						// The contribution to S_Lo from an individual light sample or ray trace sample
	wa_colour Li;							// Incoming radiance at a scatter point
	float3 LP, LN;							// Position & normal vectors of sample points on light sources
	wa_colour Li_visibility;				// The visibility between a scatter point and a light source sample
	float3 Wi;								// a unit direction vector in the sphere sound a sattering point
	wa_colour Lo(0.0);						// The outgoing radiance at A
	wa_samplerBeam *beamSampler;			// Object for generating samples along the beam
	

	//No scattering unless this is a participating medium:
	if (medium.isClear) return 0.0;


	//Split the distance into substeps roughly of length 'RAYMARCH_DIST'
	RMS = max((int)(IB_l * RAYMARCH_STEPS_PER_METRE), 1);
	reciprocal_RMS = 1.0 / (float)(RMS);


	// Create the sampling object:
	beamSampler = (wa_samplerBeam *)(new wa_samplerBeamUniform(A, B));

	//Ray march along the internal beam:
	for (int s = 0; s < RMS; s++)
	{
		RMD = beamSampler->generateSample(&reciprocal_RMLpdf);

		S = A + RMD * normalized(IB);
		S_Lo = 0.0;

 
		//Direct lighting from light sources:
		for (int l = 0; l < lights.size(); l++)
		{
			for (int i = 0; i < numLightSamples; i++)
			{
				//Get a radiance sample from the light:
				if (lights[l]->radianceSample(S, 0.0, &Li, &LP, &reciprocal_Lpdf, TRACE_BIAS))
				{
					//Compute the visibility between this scatter point and the light:
					Li_visibility = visibility(S, LP, debugPath);

					S_Lo_0 = Li_visibility * reciprocal_Lpdf * Li;
					S_Lo_0 *= medium.ss * medium.phase(acos(-1.0 * dot(normalized(LP - S), normalized(IB))));
					S_Lo_0 *= reciprocal_numLightSamples;

					S_Lo += S_Lo_0;
				}
			}
		}


		//Indirect lighting from all paths in the sphere surrounding this scatter point:
		if ((traceDepth + 1) <= TRACE_DEPTH)
		{
			wa_primitive intr_next;

			for (int i = 0; i < numSphereSamples; i++)
			{
				//Sample the sphere:
				Wi = -1.0 * sphereSampleUniform();

				//trace a ray in this direction:
				Li = radiance(S, Wi, traceDepth + 1, false, debugPath, CAMERA_RAY);

				S_Lo_0 = Li * 4.0 * M_PI * reciprocal_numSphereSamples;
				S_Lo_0 *= medium.ss * medium.phase(acos(-1.0 * dot(Wi, normalized(IB))));

				S_Lo += S_Lo_0;
			}
		}

		S_Lo *= exp(-1.0 * medium.st * RMD);
		S_Lo *= reciprocal_RMLpdf * reciprocal_RMS;

		Lo += S_Lo;
	}

	delete beamSampler;

	return Lo;
}




/*
	visibility() -		returns the visibility between two points in space. (1.0 if they are mutually visible, 0.0 if they are completely occluded)	
*/
wa_colour udp_renderer::visibility(float3 A, float3 B, wa_path *debugPath)
{
	float AB_l;					// Distance between A and B
	float3 ABn;					// The unit direction vector from A to B
	float3 P0, P1;				// Current & next position along the line from A to B
	float d = 0.0;				// Current distance along the line from A to B
	wa_ray r1;					// wa_ray object for ray tracing
	wa_primitive intr;			// wa_primitive object for ray tracing
	wa_object *object;			// wa_object pointer for each intersection event 
	wa_shader shader;			// wa_shader object for each intersection
	wa_colour V = 1.0;			// the visibility between A and B (starts at 1.0) 
	wa_shader mediumCurrent;	// The volume light transport properties of the current medium
	bool enteringSurface;		// Whether or not we are entering a surface at each intersection event
	bool hit;					// Whether or not we hit a surface with the trace() call

	//Initialise some values:
	AB_l = mag(B - A);
	ABn = normalized(B - A);
	r1.D = ABn;
	P0 = A;
	mediumCurrent = globalShader;

	//Trace from A towards B, intersecting every object in between. Stop if we have reached B or if visibility has effectively become 0.0
	while(1)
	{
		if (debugPath)
		{
			int a = 0;
		}

		//Trace a ray from the current point towards B:
		r1.P = P0;
		hit = theContent->trace(r1, &intr, TRACE_BIAS);

		if(hit) d += mag(intr.getP() - P0);
		else d = 2.0 * AB_l;

		if (debugPath)
		{
			int a = 0;
		}

		//If a surface is hit and we haven't yet reached B, proceed:
		if(d <= (AB_l + TRACE_BIAS))
		{
			P1 = intr.getP();
			object = intr.getObject();
			shader = object->getShader();
			if (dot(ABn, intr.getN()) < 0.0) enteringSurface = true;
			else enteringSurface = false;

			//If at this intersection we are exiting a surface, then the current medium is the interior of that surface:
			if (enteringSurface) mediumCurrent = globalShader;
			else mediumCurrent = intr.getObject()->getShader();
		}

		//If no surface is hit, or we have hit a surface beyond B, then we have gone beyond B. So set B as the final position:
		else
		{
			P1 = B;
			d = AB_l;
		}

		//Visibility lost to due to attenuation between P0 and P1:
		V *= volumeAttenuation(mediumCurrent, P0, P1, NULL);

		//If we have reached B, then there is no further visibility to consider:
		if (d >= (AB_l - TRACE_BIAS)) return V;

		//Visibility lost due to the surface at P1:
		if (shader.isOpaque == true) return 0.0;

		//If the visibility has fallen too low, it is effectively 0.0:
		if (V < LOW_VISIBILITY_THRESHOLD) return 0.0;

		//Prepare for the next iteration:
		P0 = P1;
		if (enteringSurface) mediumCurrent = intr.getObject()->getShader();
		else mediumCurrent = globalShader;
	}

	return V;
}

wa_colour udp_renderer::facingRatio(wa_ray Wi, int depth)
{
	wa_primitive intr;

	if (theContent->trace(Wi, &intr, TRACE_BIAS) == true) return (fabs(dot(Wi.D, intr.getN())));

	return 0.0;
}


float3 udp_renderer::reflectVector(float3 V, float3 N)
{
	return -1.0 * V + 2.0 * dot(V, N) * N;
}


/*not implemented*/
float numericalIntegral_line(float(udp_renderer::*generator)(float, float), float3(*integrand)(float))
{
	return 0;
}

/*Using the provided sample generator, numerically integrate the provided function over the unit hemisphere:*/
float numericalIntegral_hemisphere(float3 (*generator)(float3, float *), float (*integrand)(float), int N)
{
	float integral = 0;
	float3 D = {0.0, 1.0, 0.0};
	float3 sampleDir;
	float weight;

	for (int i = 0; i < N; i++)
	{
		sampleDir = generator(D, &weight);
		integral += (2.0 * M_PI / (float)(N)) * weight * integrand((acos(Y(sampleDir))));
	}

	return integral;
}

/*Using the provided sample generator, numerically integrate the provided function over part of the unit hemisphere:*/
float numericalIntegral_partialHemisphere(float maxAngle, float3(*generator)(float3, float, float *), float(*integrand)(float), int N)
{
	float integral = 0;
	float3 D = { 0.0, 1.0, 0.0 };
	float3 sampleDir;
	float weight;

	float partialHemisphereArea = 2.0 * M_PI * (1.0 - cos(maxAngle));

	for (int i = 0; i < N; i++)
	{
		sampleDir = generator(D, maxAngle, &weight);
		integral += (partialHemisphereArea / (float)(N)) * weight * integrand((acos(Y(sampleDir))));
	}

	return integral;
}

std::vector<wa_path>& udp_renderer::getRecordedPaths()
{
	return recordedPaths;
}

std::vector<wa_line>& udp_renderer::getDebugLines()
{
	return debugLines;
}

//Importance sample ray march beam length with a probability density proportional to exponential attenuation of light arriving at the start of the beam:
float udp_renderer::importanceSampleRayMarchLength_attenuation(float st, float d, float *reciprocal_pdf)
{
	float random = randFloat();

	float result = (-1.0 / st) * log(1.0 - ((1.0 - exp(-1.0 * st * d)) * random));

	*reciprocal_pdf = 1.0 / ((st / (1.0 - exp(-1.0 * st * d))) * exp(-1.0 * st * result));

	return result;
}

//Importance sample ray march beam length with a probability density proportional to the intensity of light reaching the beam from each light source:
float udp_renderer::importanceSampleRayMarchLength_lightProximity(float3 A, float3 B, float3 LP, float *reciprocal_pdf)
{
	float h;						// Closest distance between light and beam
	float lambda;					// Point on beam A->B which is closest to the light
	float theta1, theta2;			// Angles at LP, between A/B and the closest point on the beam
	float3 D;						// Unit direction vector of the beam A->B
	float x;						// The result;

	D = normalized(B - A);
	lambda = dot(D, (LP - A));
	h = mag(LP - A - (lambda * D));
	theta1 = atan(lambda / h);
	theta2 = atan((mag(B - A) - lambda) / h);

	float E = randFloat();
	x = lambda - h * tan(theta1 * (1.0 - E) - theta2 * E);
	
	*reciprocal_pdf = (pow(h, 2) + pow(lambda - x, 2)) * (theta1 + theta2) / h;

	return x;
}

//Importance sample ray march beam length with uniform probability density:
float udp_renderer::importanceSampleRayMarchLength_uniform(float d, float *reciprocal_pdf)
{
	float random = randFloat();

	float result = random * d;

	*reciprocal_pdf = 1.0 / d;
	
	return result;
}

void udp_renderer::reset()
{
	itCnt = 0;

	for (int i = 0; i < disW*disH; i++)
	{
		internalBuf[i * 3 + 0] = 0.0;
		internalBuf[i * 3 + 1] = 0.0;
		internalBuf[i * 3 + 2] = 0.0;
	}
}

void udp_renderer::setIterationCount(int itCnt_in)
{
	itCnt = itCnt_in;
}

void udp_renderer::lowerResolution(int exp, int *w_out, int *h_out)
{
	int w_new = (int)((float)(disW) / pow(2.0, exp));
	int h_new = (int)((float)(disH) / pow(2.0, exp));

	if ((disW % w_new) != 0) w_new++;
	if ((disH % h_new) != 0) h_new++;

	*w_out = w_new;
	*h_out = h_new;
}

void udp_renderer::abort()
{
	abortRender = true;
}

int udp_renderer::generateTiles(int subdivLevel, tile *renderTiles)
{
	int dim = (int)(std::pow(2, subdivLevel));

	int tw = (int)((float)(renW) / (float)(dim));
	int th = (int)((float)(renH) / (float)(dim));

	for(int i = 0; i < dim; i++)
	{
		for(int j = 0; j < dim; j++)
		{
			renderTiles[j * dim + i].l = i * tw;
			renderTiles[j * dim + i].r = (i + 1) * tw;
			renderTiles[j * dim + i].b = j * th;
			renderTiles[j * dim + i].t = (j + 1) * th;

			if(i == (dim - 1)) renderTiles[j * dim + i].r = renW;
			if(j == (dim - 1)) renderTiles[j * dim + i].t = renH;
		}
	}

	return std::pow(dim, 2);
}

void udp_renderer::setRegion(int l, int r, int b, int t)
{
	rr.l = l;
	rr.r = r;
	rr.b = b;
	rr.t = t;

	rr.w = r - l;
	rr.h = t - b;
}

udp_renderer::~udp_renderer()
{
	free(internalBuf);
}