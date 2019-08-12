#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>

#include "raytracer.h"


struct vert { float x, y, z, a; };
struct tri { int v0, v1, v2; };
struct quad { int v0, v1, v2, v3; };


void wa_raytracer::init()
{
	/*Initialize Embree:*/
	device = rtcNewDevice(NULL);
	scene = rtcNewScene(device);
}

void wa_raytracer::loadLightGeo(float px, float py, float pz, float ux, float uy, float uz, float vx, float vy, float vz)
{
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_QUAD);

	vert *vertices = (vert *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(vert), 4);
	quad *quads = (quad *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, sizeof(quad), 1);


	vertices[0].x = px; vertices[0].y = py; vertices[0].z = pz;
	vertices[1].x = px + ux; vertices[1].y = py + uy; vertices[1].z = pz + uz;
	vertices[2].x = px + ux + vx; vertices[2].y = py + uy + vy; vertices[2].z = pz + uz + vz;
	vertices[3].x = px + vx; vertices[3].y = py + vy; vertices[3].z = pz + vz;

	quads[0].v0 = 0; quads[0].v1 = 1; quads[0].v2 = 2; quads[0].v3 = 3;

	rtcCommitGeometry(geom);
	unsigned int geomID = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);
}

void wa_raytracer::loadGeo(string fileName)
{
	cout << "Loading " << fileName.c_str() << "..." << endl;
	//auto start_time = microsec_clock::local_time();

	/*Load the geometry into tinyObj format. **ASSUMING THESE ARE TRIANGLES** */
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, fileName.c_str());

	//auto td = time_period(start_time, microsec_clock::local_time()).length();
	//cout << "tinyobj loading took: " << td << endl;


	//cout << "Num Verts: " << attrib.vertices.size() << endl;
	//cout << "Num Faces: " << shapes[0].mesh.num_face_vertices.size() << endl;

	//start_time = microsec_clock::local_time();
	/*Pass the data to Embree:*/
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

	vert *vertices = (vert *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(vert), attrib.vertices.size() / 3.0);
	tri *triangles = (tri *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(tri), shapes[0].mesh.num_face_vertices.size());

	for (int v = 0; v < attrib.vertices.size() / 3; v++)
	{
		vertices[v].x = attrib.vertices[v * 3 + 0];
		vertices[v].y = attrib.vertices[v * 3 + 1];
		vertices[v].z = attrib.vertices[v * 3 + 2];
	}

	int count = 0;
	for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++)
	{
		triangles[count].v0 = shapes[0].mesh.indices[f * 3 + 0].vertex_index;
		triangles[count].v1 = shapes[0].mesh.indices[f * 3 + 1].vertex_index;
		triangles[count].v2 = shapes[0].mesh.indices[f * 3 + 2].vertex_index;

		count++;
	}

	rtcCommitGeometry(geom);
	unsigned int geomID = rtcAttachGeometry(scene, geom);
	rtcReleaseGeometry(geom);

	//td = time_period(start_time, microsec_clock::local_time()).length();
	//cout << "Embree loading took: " << td << "\n\n" << endl;
}

void wa_raytracer::updateLight(int, float, float, float, float, float, float, float, float, float)
{

}

void wa_raytracer::commit()
{
	rtcCommitScene(scene);
}

bool wa_raytracer::trace(float px, float py, float pz, float dx, float dy, float dz, float *hPx, float *hPy, float *hPz, float *hNx, float *hNy, float *hNz, int *hID, float bias) const
{
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);

	RTCRayHit rh;
	rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	rh.ray = { px, py, pz, bias, dx, dy, dz, 0.0, 10000000.0, 0, 0, 0 };

	rtcIntersect1(scene, &context, &rh);
	float rhd = rh.ray.tfar;
	*hPx = rh.ray.org_x + rhd * rh.ray.dir_x;
	*hPy = rh.ray.org_y + rhd * rh.ray.dir_y;
	*hPz = rh.ray.org_z + rhd * rh.ray.dir_z;

	if (rh.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		*hNx = rh.hit.Ng_x;
		*hNy = rh.hit.Ng_y;
		*hNz = rh.hit.Ng_z;

		*hID = rh.hit.geomID;

		return true;
	}

	return false;
}

/*
float traceShadow(float3 A, float3 B, float bias) const;
void sphereBoundsFunc(const struct RTCBoundsFunctionArguments *args);
void sphereIntersectFunc(const RTCIntersectFunctionNArguments *args);
void sphereOccludedFunc(const RTCOccludedFunctionNArguments *args);
void sphereFilterFunction(const RTCFilterFunctionNArguments *args);
*/


void wa_raytracer::release()
{

}


RTCScene wa_raytracer::getEmbreeScene()
{
	return scene;
}