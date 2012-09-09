#include "color.h"
#include "vector.h"
#include "scene.h"
#include "random_generator.h"
#include "tcpsocket.h"
#include <list>
#include "cxxptl_sdl.h"


extern const char* host;
extern int port;
extern const char* dataDir;
extern int numberOfThreads;

struct IDToSceneRecord {
	int id;
	int references;
	Scene* scene;

	IDToSceneRecord(int id, int references, Scene* scene)
		: id(id), references(references), scene(scene) {
	}
};

extern std::list<IDToSceneRecord> idToScene;
extern Mutex idToSceneMutex;

Color raytrace(const Ray& ray, Scene& scene);
bool lightIsVisible(Vector p, Vector l, Scene& scene);
Color pathtrace(const Ray& ray, const Color& pathMultiplier, Random& rnd, Scene& scene);
Color renderSample(double x, double y, Scene& scene);
Color renderPixelNoAA(int x, int y, Scene& scene);
Color renderPixelAA(int x, int y, Scene& scene);
bool lightIsVisible(Vector p, Vector l, Scene& scene);
double getClosestGeom(Ray ray, Scene& scene);
void renderNoAA(FragmentData& fd, Scene& scene);
void renderAA(FragmentData& fd, Scene& scene);

class Renderer {
	int thread_index;
	Socket socket;
	Scene* currentScene;
	int currentSceneId;
	FragmentData currentFD;
	bool isNoAAFDPosponed;
	bool isAAFDPosponed;
  
	std::list<IDToSceneRecord>::iterator Renderer::getSceneRecord(int id);
	std::list<IDToSceneRecord>::iterator Renderer::getZeroReferencedRecord();
	bool loadCurrentScene();
	void handleRenderNoAA();
	void handleRenderAA();
	void handleReceiveScene();
	void handleChangeScene();
  
public:
	Renderer(int thread_index, Socket socket);
	void run();
};

class RendererStarter : public Parallel {
public:
	virtual void entry(int thread_index, int threads_count);
};

int parseRendererArguments(int argc, char** argv);
int startRenderer(int argc, char** argv);
