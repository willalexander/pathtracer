
#define NOMINMAX

#include <cstdlib>
#include <iostream>
#include <string>

#include "unidirectionalPathTracer.h"
#include "randomnumbers.h"

#include <boost/program_options.hpp>
#include <boost/date_time.hpp>
namespace po = boost::program_options;

using namespace std;



class objTri
{
	public:
		objTri(double ux, double uy, double uz, double vx, double vy, double vz, double wx, double wy, double wz)
		{
			v[0] = ux;
			v[1] = uy;
			v[2] = uz;

			v[3] = vx;
			v[4] = vy;
			v[5] = vz;

			v[6] = wx;
			v[7] = wy;
			v[8] = wz;
		}

		double v[9];
};

void wa_parseCommandLine(int argc, char** argv, po::options_description& desc, po::variables_map& vm) {
	desc.add_options()
		("help,h", "Output this help message")
		("bmp", po::value<std::string>()->default_value("out.bmp"), "Output BMP image")
		("input,i", po::value<std::string>()->default_value("in.waff"), "Input file")
		("psamples", po::value<int>()->default_value(1), "Number of samples per pixel")
		("hsamples", po::value<int>()->default_value(1), "Number of samples in the hemisphere")
		("lsamples", po::value<int>()->default_value(1), "Number of samples for each light")
		("convertOBJ", "Perform an OBJ conversion of a file, rather than a render")
		("obj", po::value<std::string>(), "input OBJ file");

	po::positional_options_description p;
	p.add("input", -1);

	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);
}


void getUpAndFovFromFile(char const *, float3 *, float *);
void getAreaLightAreasFromFile(char const *, std::vector<float>&);
void createLightList(char const *, std::vector<wa_light *> *);


wa_camera loadCamFromFile_nff(FILE *);
std::vector<wa_light *> loadLightsFromFile_nff(FILE *, std::vector<wa_object *>);
std::vector<wa_object *> loadGeoFromFile_waff(FILE *, std::vector<wa_object *>);
wa_content *loadContentFromFile_waff(std::string);

void writeImage(std::string, int, int, float *);
void exportJrSceneToObj(char const *, const std::vector<wa_light *>&, std::vector<wa_path>&, std::vector<wa_line>&);



int main(int argc, char** argv) 
{
	/*double integral = 0.0;
	float rpdf;

	for (long i = 0; i < 10000000; i++)
	{
		sphereSampleBiased3(float3(0.0, 1.0, 0.0), 1.5, &rpdf);
		integral += (1.0/10000000.0) * (double)(rpdf);
	}

	std::cout << "Integral: " << integral << std::endl;*/

	const std::string       name = "Command line options:";
	po::options_description desc(name);
	po::variables_map       vm;
	wa_parseCommandLine(argc, argv, desc, vm);
	std::string inp_file;
	int psamples, hsamples, lsamples;
	std::string objFile, fastFile;

	if(vm.count("convertOBJ"))
	{
		objFile = vm["obj"].as<std::string>();
		fastFile = objFile.substr(0, objFile.size() - 3) + "fast";
		std::cout << "Convert OBJ mode: " << objFile << " -> " << fastFile << std::endl;
		
		wa_raytracer tmpRayTracer;
		float *vertices;
		int *triangles;
		int numVertices, numTriangles;
		tmpRayTracer.convertObjFileToEmbreeFormat(objFile, &vertices, &triangles, &numVertices, &numTriangles);

		char *data = (char *)(malloc(numVertices * 16 + numTriangles * 12 + 8));
		std::memcpy(data + 8, vertices, numVertices*16);
		std::memcpy(data + 8 + numVertices*16, triangles, numTriangles*12);
		((int *)(data))[0] = numVertices;
		((int *)(data))[1] = numTriangles;
			
		std::cout << numVertices << " " << numTriangles << " " << ((int *)(data))[0] << " " << ((int *)(data))[1] << std::endl;

		free(vertices);
		free(triangles);

		std::ofstream of(fastFile, std::ios_base::binary);
		of.write(data, numVertices * 16 + numTriangles * 12 + 8);

		free(data);

		return 0;
	}

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count("input")) {
		inp_file = vm["input"].as<std::string>();
	}
	else {
		std::cerr << "No input file specified!" << std::endl;
		std::cerr << desc << std::endl;
		return EXIT_FAILURE;
	}

	if (vm.count("psamples")) psamples = vm["psamples"].as<int>();
	if (vm.count("hsamples")) hsamples = vm["hsamples"].as<int>();
	if (vm.count("lsamples")) lsamples = vm["lsamples"].as<int>();


	//Set up scene content with Embree as the ray tracer:
	wa_content *content = loadContentFromFile_waff(inp_file);
	


	//Create a frame buffer:
	int resI, resJ;
	content->getImageResolution(&resI, &resJ);
	float *frameBuffer = (float *)(malloc(resI * resJ * 3 * sizeof(float)));

	//Create my unidirectional path tracer:
	udp_renderer theRenderer(content, frameBuffer, psamples, hsamples, 1, lsamples);

	//Start the timer:
	using namespace boost::posix_time;
	const auto start_time = microsec_clock::local_time();

	//Render:
	theRenderer.render(frameBuffer);

	//Stop the timer:
	const auto td = time_period(start_time, microsec_clock::local_time()).length();
	std::cout << "Render time: " << td << endl;

	//Apply some tone mapping:
	float minVal = 1000000;
	float maxVal = 0;

	float range = maxVal - minVal;

	minVal = 0;
	range = 1;

	wa_colour finalValue;

	//Tone mapping & LUT:
	for (int j = 0; j < resJ; j++)
	{
		for (int i = 0; i < resI; i++)
		{
			for (int c = 0; c < 3; c++)
			{
				frameBuffer[(j * resI + i) * 3 + c] = pow((frameBuffer[(j * resI + i) * 3 + c] - minVal) / range, 0.45);
				if (frameBuffer[(j * resI + i) * 3 + c] < 0.0) frameBuffer[(j * resI + i) * 3 + c] = 0.0;
				if (frameBuffer[(j * resI + i) * 3 + c] > 1.0) frameBuffer[(j * resI + i) * 3 + c] = 1.0;
			}
		}
	}

	writeImage(vm["bmp"].as<std::string>(), resI, resJ, frameBuffer);

	free(frameBuffer);

	//Export the scene and any recorded paths to .obj so that it can be viewed in blender for debugging
	std::vector<wa_path> recordedPaths = theRenderer.getRecordedPaths();
	std::vector<wa_line> debugLines = theRenderer.getDebugLines();
	std::vector<wa_light *> tmpLights;
	exportJrSceneToObj(inp_file.c_str(), tmpLights, recordedPaths, debugLines);

	return 0;
}


void getUpAndFovFromFile(char const *fname, float3 *upPtr, float *fovPtr)
{
	FILE* fp = fopen(fname, "r");

	int c;
	float x, y, z;
	float fov_angle;
	while ((c = getc(fp)) != EOF)
	{
		if (c == 'v')
		{
			fscanf(fp, " from %f %f %f", &x, &y, &z);
			fscanf(fp, " at %f %f %f", &x, &y, &z);
			fscanf(fp, " up %f %f %f", &x, &y, &z);
			fscanf(fp, " angle %f", &fov_angle);

			break;
		}
	}

	fclose(fp);

	*upPtr = {x, y, z};
	*fovPtr = fov_angle;
}

void getAreaLightAreasFromFile(char const *fname, std::vector<float>& areaLightAreas)
{
	FILE* fp = fopen(fname, "r");

	int c, c1;
	float x, y, z, ux, uy, uz, vx, vy, vz;
	float fov_angle;
	while ((c = getc(fp)) != EOF)
	{
		if (c == 'l')
		{
			c1 = getc(fp);
			if (c1 == 'a')
			{
				fscanf(fp, "%f %f %f %f %f %f %f %f %f", &x, &y, &z,
					&ux, &uy, &uz, &vx, &vy, &vz);

				areaLightAreas.push_back(sqrt(ux*ux + uy * uy + uz * uz) * sqrt(vx*vx + vy * vy + vz * vz));

				continue;
			}
		}
	}

	fclose(fp);
}




void createLightList(char const *fname, std::vector<wa_light *> *lights)
{
	FILE* fp = fopen(fname, "r");

	int c, c1;
	float x, y, z, ux, uy, uz, vx, vy, vz, cr, cg, cb;
	float R;
	float fov_angle;
	while ((c = getc(fp)) != EOF)
	{
		if (c == 'l')
		{
			c1 = getc(fp);

			//Area light
			if (c1 == 'a')
			{
				fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f", &x, &y, &z,
					&ux, &uy, &uz, &vx, &vy, &vz, &cr, &cg, &cb);

				wa_areaLight *test = new wa_areaLight(x, y, z, ux, uy, uz, vx, vy, vz, cr, cg, cb);
				wa_light *theLight = (wa_light *)(test);
				lights->push_back(theLight);

				continue;
			}

			//Sphere light
			if (c1 == 's')
			{
				fscanf(fp, "%f %f %f %f %f %f %f", &R, &x, &y, &z, &cr, &cg, &cb);

				wa_sphereLight *theLight = new wa_sphereLight(x, y, z, R, cr, cg, cb);
				wa_light *baseLight = (wa_light *)(theLight);
				lights->push_back(baseLight);

				continue;
			}

			//Point light, so small sphere light in my implementation:
			if (c1 == ' ')
			{
				cr = cg = cb = 10000;
				fscanf(fp, "%f %f %f %f %f %f", &x, &y, &z, &cr, &cg, &cb);

				wa_sphereLight *theLight = new wa_sphereLight(x, y, z, 1.0, cr, cg, cb);
				wa_light *baseLight = (wa_light *)(theLight);
				lights->push_back(baseLight);
			}
		}
	}

	fclose(fp);
}

/*
	parses the (already open) file, retrieves the camera data, returns it by reference:
*/
wa_camera loadCamFromFile_nff(FILE *fp)
{
	int c;
	int type;
	float3 p, a, u;
	float fov_angle;
	float hither;
	int resI, resJ;
	float sensorWidth;

	while ((c = getc(fp)) != EOF)
	{
		if (c == 'v')
		{
			fscanf(fp, "%d", &type);
			fscanf(fp, " from %f %f %f", &X(p), &Y(p), &Z(p));
			fscanf(fp, " at %f %f %f", &X(a), &Y(a), &Z(a));
			fscanf(fp, " up %f %f %f", &X(u), &Y(u), &Z(u));
			fscanf(fp, " angle %f", &fov_angle);
			fscanf(fp, " hither %f", &hither);
			fscanf(fp, " resolution %d %d", &resI, &resJ);
			fscanf(fp, " width %f", &sensorWidth);

			break;
		}
	}

	float3 d = normalized(a - p);

	return wa_camera(type, p, d, u, fov_angle, resI, resJ, sensorWidth);
}

/*
	parses the (already open) file, retrieves the lights, returns them by reference:
*/
std::vector<wa_light *> loadLightsFromFile_nff(FILE *fp, std::vector<wa_object *> *objects)
{
	std::vector<wa_light *> lights;

	int c, c1;
	float x, y, z, ux, uy, uz, vx, vy, vz, cr, cg, cb;
	float R;
	float fov_angle;
	while ((c = getc(fp)) != EOF)
	{
		if (c == '#')
		{
			while ((c != '\n') && (c != EOF)) c = getc(fp);
			if (c == EOF) break;

			continue;
		}

		if (c == 'l')
		{
			cout << "Found Light." << endl;

			c1 = getc(fp);
			wa_light *theLight;
			wa_shader lightShader;
			lightShader.isLight = true;
			lightShader.kd = 0.5;

			//Area light
			if (c1 == 'a')
			{
				fscanf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f", &x, &y, &z,
					&ux, &uy, &uz, &vx, &vy, &vz, &cr, &cg, &cb);

				theLight = (wa_light *)(new wa_areaLight(x, y, z, ux, uy, uz, vx, vy, vz, cr, cg, cb));
			}

			//Sphere light
			else if (c1 == 's')
			{
				fscanf(fp, "%f %f %f %f %f %f %f", &R, &x, &y, &z, &cr, &cg, &cb);

				theLight = (wa_light *)(new wa_sphereLight(x, y, z, R, cr, cg, cb));
			}

			//Point light, so small sphere light in my implementation:
			else if (c1 == ' ')
			{
				cr = cg = cb = 10000;
				fscanf(fp, "%f %f %f %f %f %f", &x, &y, &z, &cr, &cg, &cb);

				theLight = (wa_light *)(new wa_sphereLight(x, y, z, 1.0, cr, cg, cb));
			}

			

			else continue;

			lights.push_back(theLight);
			objects->push_back(new wa_object(string(""), string("diffuse"), lightShader));
		}

		if (c == 'g')
		{
			ungetc(c, fp);
			break;
		}
	}

	return lights;
}

/*
	parses the (already open) file, retrieves the geometry data, returns it by reference:
*/
void loadGeoFromFile_waff(FILE *fp, std::vector<wa_object *> *objects)
{
	int c;
	char shaderTypeName[1000];
	int type;
	char fileName[1000];
	float dr, dg, db;
	float tr, tg, tb, ar, ag, ab, er, eg, eb, sr, sg, sb;
	float ior;
	wa_shader newShader;


	while ((c = getc(fp)) != EOF)
	{
		if (c == '#')
		{
			while ((c != '\n') && (c != EOF)) c = getc(fp);
			if (c == EOF) break;

			continue;
		}

		if (c == 'g')
		{
			fscanf(fp, "%s", shaderTypeName);

			newShader.initialize();

			string shaderTypeNameString(shaderTypeName);
			if (shaderTypeNameString == "diffuse")
			{
				newShader.diffuse = true;

				fscanf(fp, "%f %f %f %s", &dr, &dg, &db, fileName);

				newShader.kd.set(dr, dg, db);
				newShader.ks.set(0.0, 0.0, 0.0);
				tr = tg = tb = 0.0;
			}

			if (shaderTypeNameString == "specular")
			{
				fscanf(fp, "%f %f %f %s", &sr, &sg, &sb, fileName);

				newShader.kd.set(0.0, 0.0, 0.0);
				newShader.ks.set(sr, sg, sb);
				tr = tg = tb = 0.0;
				newShader.diffuse = false;
				newShader.specular = true;
			}

			if (shaderTypeNameString == "sss")
			{
				newShader.kd = 0.0;
				newShader.isOpaque = false;
				newShader.isClear = false;
				newShader.diffuse = newShader.specular = newShader.dielectric = false;

				fscanf(fp, "%d", &type);
				if (type == 0)
				{
					newShader.proceduralVolume = 0;
					newShader.homogenous = false;
					tr = tg = tb = 1.0;
				}

				else
				{
					fscanf(fp, "%f %f %f  %f %f %f  %f %f %f  %f %f %f", &tr, &tg, &tb, &ar, &ag, &ab, &er, &eg, &eb, &sr, &sg, &sb);
				}

				fscanf(fp, "%s", fileName);
			}

			if (shaderTypeNameString == "dielectric")
			{
				newShader.kd = 0.0;
				newShader.ks = 1.0;
				newShader.isOpaque = false;
				newShader.isClear = true;
				newShader.diffuse = newShader.specular = newShader.dielectricsss = false;
				newShader.dielectric = true;
				newShader.kt = tr = tg = tb = 1.0;

				fscanf(fp, "%d", &type);
				if (type == 0)
				{
					newShader.proceduralVolume = 0;
					newShader.homogenous = false;
					tr = tg = tb = 1.0;
				}

				else
				{
					fscanf(fp, "%f", &ior);
				}

				fscanf(fp, "%s", fileName);

				std::cout << "Dielectric: " << ior << std::endl;
			}

			if (shaderTypeNameString == "dielectricsss")
			{
				newShader.kd = 0.0;
				newShader.ks = 1.0;
				newShader.isOpaque = false;
				newShader.isClear = false;
				newShader.diffuse = newShader.specular = newShader.dielectric = false;
				newShader.dielectricsss = true;
				newShader.kt = 1.0;

				fscanf(fp, "%d", &type);
				if (type == 0)
				{
					newShader.proceduralVolume = 0;
					newShader.homogenous = false;
					tr = tg = tb = 1.0;
				}

				else
				{
					fscanf(fp, "%f  %f %f %f  %f %f %f  %f %f %f  %f %f %f", &ior, &tr, &tg, &tb, &ar, &ag, &ab, &er, &eg, &eb, &sr, &sg, &sb);
				}

				fscanf(fp, "%s", fileName);

				std::cout << "Scattering: " << sr << " " << sg << " " << sb << std::endl;
			}

			newShader.ior = ior;
			newShader.kt.set(tr, tg, tb);
			newShader.sa.set(ar, ag, ab);
			newShader.se.set(er, eg, eb);
			newShader.ss.set(sr, sg, sb);
			newShader.st.set(ar + sr, ag + sg, ab + sb);
			objects->push_back(new wa_object(string(fileName), string(shaderTypeName), newShader));

			cout << "Found geo" << endl;
		}
	}
}

/*
	Loads up camera, lights and geometry from a WAFF file:
*/
wa_content *loadContentFromFile_waff(std::string fileName)
{
	FILE* fp = fopen(fileName.c_str(), "r");

	std::vector<wa_object *> objects;

	wa_camera camera = loadCamFromFile_nff(fp);
	std::vector<wa_light *> lights = loadLightsFromFile_nff(fp, &objects);

	cout << camera << endl;
	for (int i = 0; i < lights.size(); i++)
	{
		lights[i]->print();
	}

	loadGeoFromFile_waff(fp, &objects);
	
	fclose(fp);

	wa_content *content = new wa_content(camera, lights, objects);
	content->init();

	return content;
}

//Write the final image out as a .bmp:
void writeImage(std::string filePath, int w, int h, float *pixels)
{
	char *data = (char *)(malloc(138 + w*h*4 * sizeof(char)));
	for (int i = 0; i < (138 + w * h * 4); i++) data[i] = 0;

	*((uint16_t *)(data)) = 0x4D42;								// BMP format
	*((uint32_t *)(data + 2)) = 138 + w*h*4;					// Total image size in bytes
	*((uint32_t *)(data + 10)) = 138;							// Offset to pixel data

	*((uint32_t *)(data + 14)) = 124;							// Header size (not including the parts above this line)
	*((int32_t *)(data + 18)) = w;
	*((int32_t *)(data + 22)) = h;
	*((uint16_t *)(data + 26)) = 1;								// planes
	*((uint16_t *)(data + 28)) = 32;							// bits per pixel
	*((uint32_t *)(data + 30)) = 3;								// uncompressed

	*((uint32_t *)(data + 54)) = 0x00ff0000;					// RGBA channels bit masks
	*((uint32_t *)(data + 58)) = 0x0000ff00;
	*((uint32_t *)(data + 62)) = 0x000000ff;
	*((uint32_t *)(data + 66)) = 0xff000000;
	*((uint32_t *)(data + 70)) = 0x73524742;					// sRGB


	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			int pixel[4] = { pixels[(j*w + i)*3 + 0] * 255.0, pixels[(j*w + i) * 3 + 1] * 255.0 , pixels[(j*w + i) * 3 + 2] * 255.0 , 255};
			
			data[138 + (j * w + i) * 4 + 0] = (char)(pixel[2]);
			data[138 + (j * w + i) * 4 + 1] = (char)(pixel[1]);
			data[138 + (j * w + i) * 4 + 2] = (char)(pixel[0]);
			data[138 + (j * w + i) * 4 + 3] = 255;
		}
	}

	std::ofstream of(filePath, std::ios_base::binary);
	of.write(data, 138 + w * h * 4);

	free(data);
}

void exportJrSceneToObj(char const *fname, const std::vector<wa_light *>& lights, std::vector<wa_path>& recordedPaths, std::vector<wa_line>& debugLines)
{
	//Now write the geometric info into an .obj file:
	FILE* fpo = fopen("debug_content.obj", "w");
	int vc = 1;

	//Write recorded path verts:
	for (int p = 0; p < recordedPaths.size(); p++)
	{
		for (int v = 0; v < recordedPaths[p].getSize(); v++)
		{
			float3 pathVert = recordedPaths[p].getVert(v);
			fprintf(fpo, "v %f %f %f\n", X(pathVert), Y(pathVert), Z(pathVert));
		}
	}

	for (int p = 0; p < recordedPaths.size(); p++)
	{
		fprintf(fpo, "l ");

		for (int v = 0; v < recordedPaths[p].getSize(); v++)
		{
			fprintf(fpo, "%i ", vc);
			vc++;
		}

		fprintf(fpo, "\n");
	}


	fclose(fpo);
}
