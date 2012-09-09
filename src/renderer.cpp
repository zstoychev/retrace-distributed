#include <SDL/SDL.h>
#include <SDL/SDL_net.h>
#include <time.h>
#include <list>
#include <utility>
#include "sdl.h"
#include "matrix.h"
#include "lights.h"
#include "camera.h"
#include "geometry.h"
#include "shading.h"
#include "environment.h"
#include "random_generator.h"
#include "cxxptl_sdl.h"
#include <string>
#include <cstdio>
#include "tcpsocket.h"
#include "cxxptl_sdl.h"
#include "renderer.h"
#include "util.h"

/// traces a ray in the scene and returns the visible light that comes from that direction
Color raytrace(const Ray& ray, Scene& scene)
{
	if (ray.depth > scene.settings.maxTraceDepth) return Color(0, 0, 0);
	if (ray.flags & RF_DEBUG) {
		printf("  Raytrace[start = ");
		ray.start.print();
		printf(", dir = ");
		ray.dir.print();
		printf("]\n");
	}
	IntersectionInfo closestInfo;
	closestInfo.distance = INF;
	Node* closestNode = NULL;
	for (int i = 0; i < (int) scene.nodes.size(); i++) {
		IntersectionInfo info;
		if (scene.nodes[i]->intersect(ray, info)) {
			if (info.distance < closestInfo.distance) {
				closestInfo = info;
				closestNode = scene.nodes[i];
			}
		}
	}
	// check if the closest intersection point is actually a light:
	bool hitLight = false;
	Color hitLightColor;
	for (int i = 0; i < (int) scene.lights.size(); i++) {
		double d = scene.lights[i]->intersect(ray);
		if (d < closestInfo.distance) {
			closestInfo.distance = d;
			hitLight = true;
			hitLightColor = scene.lights[i]->col;
		}
	}
	if (hitLight) return hitLightColor;
	if (!closestNode) {
		if (scene.environment != NULL)
			return scene.environment->getEnvironment(ray.dir);
		return Color(0, 0, 0);
	} else {
		if (ray.flags & RF_DEBUG) {
			printf("    Closest node is a %s at distance %.3lf\n", closestNode->geometry->getName(), closestInfo.distance);
			printf("      ip   = "); closestInfo.ip.println();
			printf("      norm = "); closestInfo.norm.println();
		}

		Color result = closestNode->shader->shade(ray, closestInfo);
		if (closestNode->active) {
			result *= 1.6f;
		}

		return result;
	}
}

bool lightIsVisible(Vector p, Vector l, Scene& scene);

/// same as raytrace(), but return the result of a single path-tracing path, using global illumination
/// @param ray - the ray to be traced
/// @param pathMultiplier - the current product of BRDFs along the path; should start with (1,1,1)
/// @param rnd - a random generator (for speed and convenience only)
Color pathtrace(const Ray& ray, const Color& pathMultiplier, Random& rnd, Scene& scene)
{
	Color result(0, 0, 0);
	if (ray.depth > scene.settings.maxTraceDepth) return result;
	IntersectionInfo closestInfo;
	closestInfo.distance = INF;
	Node* closestNode = NULL;
	for (int i = 0; i < (int) scene.nodes.size(); i++) {
		IntersectionInfo info;
		if (scene.nodes[i]->intersect(ray, info)) {
			if (info.distance < closestInfo.distance) {
				closestInfo = info;
				closestNode = scene.nodes[i];
			}
		}
	}
	// check if the closest intersection point is actually a light:
	bool hitLight = false;
	Color hitLightColor;
	for (int i = 0; i < (int) scene.lights.size(); i++) {
		double d = scene.lights[i]->intersect(ray);
		if (d < closestInfo.distance) {
			closestInfo.distance = d;
			hitLight = true;
			hitLightColor = scene.lights[i]->col;
		}
	}
	if (hitLight) {
		/*
		 * if the ray actually hit a light, check if we need to pass this light back along the path.
		 * If the last surface along the path was a diffuse one (Lambert/Phong), we need to discard the
		 * light contribution, since for diffuse material we do explicit light sampling too, thus the
		 * light would be over-represented and the image a bit too bright. We may discard light checks
		 * for secondary rays altogether, but we would lose caustics and light reflections that way.
		 */
		if (!(ray.flags & RF_DIFFUSE)) result = pathMultiplier * hitLightColor;
		return result;
	}
	if (!closestNode) {
		if (scene.environment != NULL)
			return scene.environment->getEnvironment(ray.dir) * pathMultiplier;
		return Color(0, 0, 0);
	}
	// We continue building the path in two ways:
	// 1) (a.k.a. "direct illumination"): connect the current path end to a random light.
	//    This approximates the direct lighting towards the intersection point.
	Color resultDirect(0.0f, 0.0f, 0.0f);
	if (scene.lights.size() > 0) {
		const Vector& x = closestInfo.ip; // intersection point
		const Vector& xn = closestInfo.norm;
		Vector norm = faceforward(xn, ray.dir); // intersection point normal
		
		// choose a random light:
		int lightIdx = rnd.randint(0, scene.lights.size() - 1);
		Light* light = scene.lights[lightIdx];
		
		// choose a random sample of that light:
		int sampleIdx = rnd.randint(0, light->getNumSamples() - 1);
		
		// sample the light and see if it came out nonzero:
		Vector lightSample;
		Color lightCol;
		light->getNthSample(sampleIdx, x, lightSample, lightCol);
		
		if (lightCol.intensity() > 0) {
			// w_out - the outgoing ray in the BRDF evaluation
			Vector w_out = lightSample - x;
			float w_out_lengthSqr = (float) w_out.lengthSqr();
			w_out.normalize();
			// evaluate the BRDF:
			Color brdf = closestNode->shader->eval(closestInfo, ray, w_out);
			// check light visibility at the intersection point:
			if (brdf.intensity() > 0 && lightIsVisible(x, lightSample, scene)) {
				// calculate the light contribution in a manner, consistent with classic path tracing:
				float solidAngle = light->solidAngle(x); // solid angle of the light, as seen from x.
				if (solidAngle > 1e-6) {
					// solidAngle / 2PI is the probability to hit that light when shooting a random ray
					lightCol = light->col * solidAngle / (2*(float)PI);
				} else {
					// this turns out basically the same as in the area case, when the area approaches zero
					lightCol = lightCol / (1 + w_out_lengthSqr);
				}
				// the probability to choose a particular light among all lights: 1/N
				float pdfLight = 1 / (float) scene.lights.size();
				// the probability to shoot a ray in a random direction: 1/2*pi
				float pdfHemisphere = 1 / (2*(float)PI);
				float pdf = pdfLight * pdfHemisphere; // combined probability
				// Kajia's rendering equation, evaluated at a single incoming/outgoing directions pair:
				               /* Li */   /*BRDFs@path*/ /*BRDF*/             /* (-w',n) */          /*ray probability*/
				resultDirect = lightCol * pathMultiplier * brdf * ((float) max(0.0, dot(w_out, norm)) /  pdf);
				
			}
		}
	}
	
	// 2) (a.k.a. "indirect illumination"): continue the path randomly, by asking the
	//    BRDF to choose a continuation direction
	Ray newRay;
	Color brdf; // brdf at the chosen direction
	float pdf;  // the probability to choose that specific newRay
	Color resultGi;
	// sample the BRDF:
	closestNode->shader->spawnRay(closestInfo, ray, newRay, brdf, pdf);
	if (pdf < 0) resultGi = Color(1, 0, 0); // bogus BRDF; mark in red
	else if (pdf == 0) resultGi = Color(0, 0, 0); // terminate the path, as required
	else resultGi = pathtrace(newRay, pathMultiplier * brdf / pdf, rnd, scene); // continue the path normally; accumulate the new term to the BRDF product
	result = resultDirect + resultGi;
	return result;
}


// trace a ray through pixel coords (x, y). In non-stereoscopic mode fetches a screen
// ray and raytrace()s it.
// In stereoscopic mode, a ray is traced through both cameras, and the results are
// mixed to create an anaglyph image.
Color renderSample(double x, double y, Scene& scene)
{
	if (scene.camera->stereoSeparation == 0) {
		return raytrace(scene.camera->getScreenRay(x, y), scene);
	} else {
		Color leftEye = raytrace(scene.camera->getScreenRay(x, y, -1), scene);
		Color rightEye= raytrace(scene.camera->getScreenRay(x, y, +1), scene);
		leftEye.desaturate(0.8f);
		rightEye.desaturate(0.8f);
		return leftEye * scene.camera->leftMask + 
		       rightEye* scene.camera->rightMask;
	}
}

// gets the color for a single pixel, without antialiasing
// (when DOF is enabled, no antialiasing is really needed)
Color renderPixelNoAA(int x, int y, Scene& scene)
{
	if (!scene.camera->dof && !scene.settings.gi) {
		return renderSample(x, y, scene);
	} else {
		Random& rnd = getRandomGen();
		int n = scene.settings.gi ? scene.settings.pathsPerPixel : scene.camera->samples;
		Color accum(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < n; i++) {
			double xoff = rnd.randdouble();
			double yoff = rnd.randdouble();
			Ray ray = scene.camera->getScreenRay(x + xoff, y + yoff);
			if (scene.settings.gi)
				accum += pathtrace(ray, Color(1, 1, 1), rnd, scene);
			else
				accum += raytrace(ray, scene);
		}
		return accum / float(n);
	}
}

// gets the color for a single pixel, with antialiasing. Assumes the pixel
// already holds some value.
// This simply adds four more AA samples and averages the result.
Color renderPixelAA(int x, int y, Scene& scene)
{
	const double offsets[5][2] = {
		{ 0, 0 }, {0.3, 0.3}, {0.6, 0}, {0, 0.6}, {0.6, 0.6}
	};
	Color accum(0, 0, 0);
	for (int samples = 1; samples < 5; samples++) {
		accum += renderSample(x + offsets[samples][0], y + offsets[samples][1], scene);
	}
	
	return accum;
}

/// checks if light (situated at point l) is visible at point p. This works
/// by tracing a ray along the two points and testing whether it is unobstructed.
bool lightIsVisible(Vector p, Vector l, Scene& scene)
{
	Vector LP = p - l;
	double len = LP.length();
	Ray ray;
	ray.start = l;
	ray.dir = LP;
	ray.dir.normalize(); // save the length of the LP
	ray.flags |= RF_SHADOW;
	for (int i = 0; i < (int) scene.nodes.size(); i++) if (!scene.nodes[i]->skip_shadow) {
		IntersectionInfo info;
		if (scene.nodes[i]->intersect(ray, info) 
		    && info.distance < len - 1e-6) return false; // a hit point was found and it's closer
		                                                 // to the light than length(LP); we're in shadow.
	}
	return true;
}

double getClosestGeom(Ray ray, Scene& scene)
{
	IntersectionInfo info;
	double closestDist = INF;
	for (int i = 0; i < (int) scene.nodes.size(); i++)
		if (scene.nodes[i]->intersect(ray, info))
			closestDist = min(closestDist, info.distance);
	return closestDist;
}

void renderNoAA(FragmentData& fd, Scene& scene) {
	const Rect& r = fd.fragment;
	for (int y = r.y0; y < r.y1; y++) {
		for (int x = r.x0; x < r.x1; x++) {
			fd.data[y - r.y0][x - r.x0] = renderPixelNoAA(x, y, scene);
		}
	}
}

void renderAA(FragmentData& fd, Scene& scene) {
	const Rect& r = fd.fragment;
	for (int y = r.y0; y < r.y1; y++) {
		for (int x = r.x0; x < r.x1; x++) {
			if (fd.needsAA.get(y - r.y0, x - r.x0)) {
				fd.data[y - r.y0][x - r.x0] = renderPixelAA(x, y, scene);
			}
		}
	}
}

extern const char* host;
extern int port;
extern const char* dataDir;
int numberOfThreads = get_processor_count();

std::list<IDToSceneRecord> idToScene;
typedef std::list<IDToSceneRecord>::iterator IDToSceneIt;
Mutex idToSceneMutex;

void freeIDToScene() {
	while (!idToScene.empty()) {
		delete idToScene.front().scene;
		idToScene.pop_front();
	}
}

IDToSceneIt Renderer::getSceneRecord(int id) {
	IDToSceneIt it;
	for (it = idToScene.begin(); it != idToScene.end(); it++) {
		if (it->id == id) {
			return it;
		}
	}

	return it;
}

IDToSceneIt Renderer::getZeroReferencedRecord() {
	IDToSceneIt it;
	for (it = idToScene.begin(); it != idToScene.end(); it++) {
		if (it->references == 0) {
			return it;
		}
	}

	return it;
}
	
bool Renderer::loadCurrentScene() {
	Lock lock(idToSceneMutex);
	IDToSceneIt it = getSceneRecord(currentSceneId);

	if (it != idToScene.end()) {
		currentScene = it->scene;
		it->references++;
		idToScene.push_back(*it);
		idToScene.erase(it);
		
		return true;
	} else {
		return false;
	}
}

void Renderer::handleRenderNoAA() {
	if (!isNoAAFDPosponed) {
		socket.read(&currentFD.fragment, sizeof(currentFD.fragment));
	} else {
		isNoAAFDPosponed = false;
	}

	if(currentScene == 0) {
		if (!loadCurrentScene()) {
			isNoAAFDPosponed = true;
			socket.writeLine("RQSCENE");
			return;
		}
	}

	renderNoAA(currentFD, *currentScene);
	socket.writeLine("RCVFRAGMENT");
	socket.write(&currentFD.data, sizeof(currentFD.data));

	socket.writeLine("RQFRAGMENT");
}

void Renderer::handleRenderAA() {
	if (!isAAFDPosponed) {
		socket.read(&currentFD.fragment, sizeof(currentFD.fragment));
		socket.read(&currentFD.needsAA, sizeof(currentFD.needsAA));
	} else {
		isAAFDPosponed = false;
	}

	if(currentScene == 0) {
		if (!loadCurrentScene()) {
			isAAFDPosponed = true;
			socket.writeLine("RQSCENE");
			return;
		}
	}

	renderAA(currentFD, *currentScene);

	socket.writeLine("RCVFRAGMENT");
	socket.write(&currentFD.data, sizeof(currentFD.data));

	socket.writeLine("RQFRAGMENT");
}

void Renderer::handleReceiveScene() {
	int length;
	socket.read(&length, sizeof(length));

	char* content = new char[length + 1];
	socket.read(content, length, 1000);
	content[length] = '\0';

	if (loadCurrentScene()) {
		return;
	}

	char filename[256];
	char prefix[256];
	sprintf(prefix, "%s/renderer%d_", dataDir, thread_index);
	findEmptyFN(filename, prefix, ".retrace");
	FILE* sceneOutput = fopen(filename, "w");
	if (sceneOutput == 0) {
		printf("Cannot access %s\n", filename);
		return;
	}
	if (fwrite(content, sizeof(char), length, sceneOutput) < (size_t) length) {
		printf("Error while writing received scene to %s\n", filename);
		fclose(sceneOutput);
		return;
	}
	fclose(sceneOutput);
	delete content;

	currentScene = new Scene();
	if (!currentScene->parseScene(filename)) {
		printf("Could not parse the scene!\n");
		delete currentScene;
		currentScene = 0;
		return;
	}
	currentScene->beginRender();
	currentScene->beginFrame();
		
	if (remove(filename) != 0) {
		printf("Cannot remove working file %s\n", filename);
	}

	{
		Lock lock(idToSceneMutex);
		IDToSceneIt it = getSceneRecord(currentSceneId);
		if (it != idToScene.end()) {
			delete currentScene;
			currentScene = it->scene;
			it->references++;
			idToScene.push_back(*it);
			idToScene.erase(it);
		} else {
			idToScene.push_back(IDToSceneRecord(currentSceneId, 1, currentScene));

			IDToSceneIt it;
			while (idToScene.size() > MAX_SIMULTANEOUS_SCENES
				&& (it = getZeroReferencedRecord()) != idToScene.end()) {
				delete it->scene;
				idToScene.erase(it);
			}

			printf("Scene #%d received.\n", currentSceneId);
		}
	}
}

void Renderer::handleChangeScene() {
	if (currentScene != 0) {
		Lock lock(idToSceneMutex);
		getSceneRecord(currentSceneId)->references--;
	}
	socket.read(&currentSceneId, sizeof(currentSceneId));
	currentScene = 0;

	printf("Renderer #%d's working scene changed to scene #%d\n", thread_index, currentSceneId);
}

Renderer::Renderer(int thread_index, Socket socket)
	: thread_index(thread_index), socket(socket), currentScene(0),
	currentSceneId(-1), isNoAAFDPosponed(false), isAAFDPosponed(false) {
}

void Renderer::run() {
	printf("Renderer #%d connected to %s:%d\n", thread_index, host, port);

	try {
		socket.writeLine("RQFRAGMENT");

		while (true) {
			char message[16];
			if (socket.readLine(message, sizeof(message))) {
				if (strcmp(message, "CHGSCENE") == 0) {
					handleChangeScene();
				} else if (strcmp(message, "RCVSCENE") == 0) {
					handleReceiveScene();
					if (isNoAAFDPosponed) {
						handleRenderNoAA();
					} else if (isAAFDPosponed) {
						handleRenderAA();
					}
				} else if (strcmp(message, "RENDERNOAA") == 0) {
					handleRenderNoAA();
				} else if (strcmp(message, "RENDERAA") == 0) {
					handleRenderAA();
				} else if (strcmp(message, "QUIT") == 0) {
					break;
				}
			}				
		}
	} catch(NetworkException&) {
		printf("Network exception in Renderer #%d's connection\n", thread_index);
	}
}

void RendererStarter::entry(int thread_index, int threads_count) {
	try {
		printf("Renderer #%d connecting to %s:%d...\n", thread_index, host, port);
		Socket socket(host, port);
		Renderer renderer(thread_index, socket);
		renderer.run();
	} catch (NetworkException) {
		printf("Renderer #%d cannot connect to %s:%d\n", thread_index, host, port);
	}
}

extern void printUsage();

int parseRendererArguments(int argc, char** argv) {
	if (argc < 3) {
		printUsage();
		return 1;
	}
	if (strcmp(argv[1], "-c") != 0) {
		printUsage();
		return 1;
	}
	host = argv[2];
	
	for (int i = 3; i < argc; i++) {
		if (strcmp(argv[i], "-port") == 0 && i < argc - 1) {
			port = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-dataDir") == 0 && i < argc - 1) {
			dataDir = argv[++i];
		} else if (strcmp(argv[i], "-numberOfThreads") == 0 && i < argc - 1) {
			numberOfThreads = atoi(argv[++i]);
		} else {
			printUsage();
			return 1;
		}
	}

	return 0;
}

int startRenderer(int argc, char** argv) {
	int result = parseRendererArguments(argc, argv);
	if (result != 0) {
		return result;
	}

	initializeSDL();

	ThreadPool pool;
	RendererStarter rendererStarter;
	pool.run(&rendererStarter, numberOfThreads);

	freeIDToScene();

	return 0;
}
