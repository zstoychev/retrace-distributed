#include <SDL/SDL.h>
#include "scene.h"
#include "random_generator.h"
#include "camera.h"
#include "shading.h"
#include <vector>
#include <cmath>

Node* activeNode = 0;
int activeNodeIndex = -1;
float activeNodeBrightnessMultiplier;

const double DEFAULT_CT = 5.0;
const double DEFAULT_CR = 30.0;
const double DEFAULT_CS = 0.2;
double Ct = DEFAULT_CT, Cr = DEFAULT_CR, Cs =DEFAULT_CS;
bool inverseTransformation = false;

const char * const UNIQUE_PREFIX = "$";

double lastRenderTime;

inline double getCt() {
	return Ct * (inverseTransformation ? -1.0 : 1.0);
}

inline double getCr() {
	return Cr * (inverseTransformation ? -1.0 : 1.0);
}

inline double getCs() {
	return inverseTransformation ? 1 / (1.0 + Cs) : (1.0 + Cs);
}

inline int getActiveNodeShaderIndex() {
	for (size_t i = 0; i < currentScene->shaders.size(); i++) {
		if (currentScene->shaders[i] == activeNode->shader) {
			return i;
		}
	}

	return -1;
}

template <typename T>
void findUniqueName(const std::vector<T*>& elements, const char* type, char* destination) {
	int number = 0;
	while (true) {
		sprintf(destination, "%s%s%d", UNIQUE_PREFIX, type, number);

		bool nameExists = false;

		for (size_t i = 0; i < elements.size(); i++) {
			if (strcmp(destination, elements[i]->name) == 0) {
				nameExists = true;
				break;
			}
		}

		if (nameExists) {
			number++;
		} else {
			break;
		}
	}
}

// handles the keyboard and mouse events in interactive mode
// the `dt' param is the frame render time, in seconds.
void handleKbdMouse(bool& mustExit, double dt)
{
	if (mustExit) return;
	static bool running = false;
	SDL_Event ev;
	while (SDL_PollEvent(&ev)) {
		Geometry* newGeometry = 0;
		switch (ev.type) {
			case SDL_QUIT:
				mustExit = true;
				return;
			case SDL_KEYDOWN:
			{
				switch (ev.key.keysym.sym) {
					case SDLK_ESCAPE:
						mustExit = true;
						return;
					case 'p':
					{
						// an utility function, that tells us where the camera is.
						printf("Camera position: (%.3lf, %.3lf, %.3lf)\n", currentScene->camera->pos.x, currentScene->camera->pos.y, currentScene->camera->pos.z);
						printf("   yaw: %.3lf\n", currentScene->camera->yaw);
						printf(" pitch: %.3lf\n", currentScene->camera->pitch);
						printf("  roll: %.3lf\n", currentScene->camera->roll);
						break;
					}
					case 'r':
					{
						// activate/deactivate running (camera movement is 4x faster when running)
						running = !running;
						break;
					}
					case '[':
						if (activeNodeIndex != -1) {
							activeNodeIndex = (currentScene->nodes.size() + activeNodeIndex - 1) % currentScene->nodes.size();
							activeNode->active = false;
							activeNode = currentScene->nodes[activeNodeIndex];
							activeNode->active = true;
						}
						break;
					case ']':
						if (activeNodeIndex != -1) {
							activeNodeIndex = (activeNodeIndex + 1) % currentScene->nodes.size();
							activeNode->active = false;
							activeNode = currentScene->nodes[activeNodeIndex];
							activeNode->active = true;
						}
						break;
					case 'q':
						activeNode->T.translate(Vector(getCt(), 0, 0));
						break;
					case 'a':
						activeNode->T.translate(Vector(0, getCt(), 0));
						break;
					case 'z':
						activeNode->T.translate(Vector(0, 0, getCt()));
						break;
					case 'w':
						activeNode->T.rotate(getCr(), 0, 0);
						break;
					case 's':
						activeNode->T.rotate(0, getCr(), 0);
						break;
					case 'x':
						activeNode->T.rotate(0, 0, getCr());
						break;
					case 'e':
						activeNode->T.scale(getCs(), 1, 1);
						break;
					case 'd':
						activeNode->T.scale(1, getCs(), 1);
						break;
					case 'c':
						activeNode->T.scale(1, 1, getCs());
						break;
					case '-':
						inverseTransformation = !inverseTransformation;
						break;
					case 't':
						Lambert* lambertShader;
						if ((lambertShader = dynamic_cast<Lambert*>(activeNode->shader)) != 0) {
							Random& rnd = getRandomGen();
							lambertShader->setColor(Color(rnd.randfloat(), rnd.randfloat(), rnd.randfloat()));
						}
						break;
					case SDLK_DELETE:
						if (activeNodeIndex != -1) {
							currentScene->nodes.erase(currentScene->nodes.begin() + activeNodeIndex);
							if (currentScene->nodes.empty()) {
								activeNodeIndex = -1;
							} else {
								activeNodeIndex %= currentScene->nodes.size();
								activeNode = currentScene->nodes[activeNodeIndex];
								activeNode->active = true;
							}
						}
						break;
					case 'v':
						activeNode->shader = currentScene->shaders[(currentScene->shaders.size() + getActiveNodeShaderIndex() - 1) % currentScene->shaders.size()];
						break;
					case 'b':
						activeNode->shader = currentScene->shaders[(getActiveNodeShaderIndex() + 1) % currentScene->shaders.size()];
						break;
					case 'f':
						newGeometry = new Sphere();
						break;
					case 'g':
						newGeometry = new Cube();
						break;
					case 'h':
						newGeometry = new Cylinder();
						break;
					case 'j':
						newGeometry = new Torus();
						break;
					case SDLK_F2:
						{
							extern const char* outputDir;
							std::string filename;
							filename += outputDir;
							filename += "/interactive.retrace";
							currentScene->saveToFile(filename.c_str());
							break;
						}
					default:
						break;
				}
				
				std::string nodeChangeKeys = "qazwsxedc";
				if (nodeChangeKeys.find(ev.key.keysym.sym) != std::string::npos) {
					activeNode->beginRender();
				} else if (ev.key.keysym.sym >= '1' && ev.key.keysym.sym <= '9') {
					int num = ev.key.keysym.sym - '0';
					double scale = (pow(2.0, num - 5.0));
					Ct = DEFAULT_CT * scale;
					Cr = DEFAULT_CR * scale;
					Cs = DEFAULT_CS * scale;
				} else if (newGeometry != 0) {
					findUniqueName(currentScene->geometries, "g", newGeometry->name);
					currentScene->geometries.push_back(newGeometry);

					Random& rnd = getRandomGen();
					Shader* newShader = new Lambert(Color(rnd.randfloat(), rnd.randfloat(), rnd.randfloat()));
					findUniqueName(currentScene->shaders, "s", newShader->name);
					currentScene->shaders.push_back(newShader);

					Node* newNode = new Node();
					newNode->geometry = newGeometry;
					newNode->shader = newShader;
					newNode->T.translate(currentScene->camera->pos + currentScene->camera->getCenterRay().dir * 10);
					findUniqueName(currentScene->nodes, "n", newNode->name);
					currentScene->nodes.push_back(newNode);

					newNode->beginRender();
					activeNodeIndex = currentScene->nodes.size() - 1;
					activeNode = newNode;
				}
				break;
			}
		}
	}
	const double KEYBOARD_SENSITIVITY = 1.0;
	const double MOUSE_SENSITIVITY = 0.01;
	Uint8 *keystate;
	int deltax, deltay;
	keystate = SDL_GetKeyState(NULL);
	double M = dt * (running ? 80:20);
	double R = dt * KEYBOARD_SENSITIVITY;
	// handle arrow keys (camera movement)
	if (keystate[SDLK_UP	]) currentScene->camera->move( M, 0);
	if (keystate[SDLK_DOWN	]) currentScene->camera->move(-M, 0);
	if (keystate[SDLK_LEFT	]) currentScene->camera->move(0, -M);
	if (keystate[SDLK_RIGHT	]) currentScene->camera->move(0, +M);
	// handle keypad keys (camera lookaround)
	if (keystate[SDLK_KP2	]) currentScene->camera->rotate(0, -R);
	if (keystate[SDLK_KP4	]) currentScene->camera->rotate(-R, 0);
	if (keystate[SDLK_KP6	]) currentScene->camera->rotate(+R, 0);
	if (keystate[SDLK_KP8	]) currentScene->camera->rotate(0, +R);
	
	// handle mouse movement (camera lookaround)
	SDL_GetRelativeMouseState(&deltax, &deltay);
	currentScene->camera->rotate(MOUSE_SENSITIVITY * deltax, -MOUSE_SENSITIVITY * deltay);
}
