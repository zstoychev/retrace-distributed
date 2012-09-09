#include "SDL/SDL.h"
#include <queue>
#include "cxxptl_sdl.h"
#include "scene.h"
#include "tcpsocket.h"

extern bool showUpdates;
extern int port;
extern Uint32 fragmentRenderingTimeout;
extern bool outputFiles;
extern const char* outputDir;
extern std::queue<std::string> scenesToRender;

extern volatile bool finishRendering;
extern volatile bool disablePreview;
extern volatile bool requestResolutionChange;
extern volatile bool inInteractive;
extern Event resolutionChangeEvent;
extern Event interactiveRequestChangeFrameEvent;
extern Event interactiveFrameChangedEvent;

extern std::queue<FragmentData*> finishedFragments;
extern Mutex finishedFragmentsMutex;
extern Condition finishedFragmentsCondition;

class RendererCommunicator;
extern std::queue<RendererCommunicator*> finishedCommunicators;
extern Mutex finishedCommunicatorsMutex;
extern Condition finishedCommunicatorsCondition;

class SceneManager;
extern SceneManager sceneManager;
extern Mutex sceneFragmentsMutex;
extern Condition sceneFragmentsCondition;

class UnexpectedBehaviourException : public std::exception {
};

template <typename T> void appendToQueue(std::queue<T>& destination, std::queue<T>& source);

class SceneManager {
	std::queue<FragmentData*>* sceneFragments[MAX_SIMULTANEOUS_SCENES];
	SceneProperties* renderingScenes[MAX_SIMULTANEOUS_SCENES];
	int renderingScenesCount;
	std::queue<std::string> queuedScenes;
	int totalScenes;
	bool disableAddingScenes;
  
	void addFragments(SceneProperties* properties, std::queue<FragmentData*>* sceneFragmentsQueue);
	void startRenderingNextScene();
	int getSceneIndex(SceneProperties* properties);
	void setCurrentScene();
  
public:
	SceneManager();
	~SceneManager();
	void addSceneFile(const char* fileName);
	void rerenderScene(SceneProperties* properties);
	void addAAFragments(SceneProperties* properties);
	void stopRenderingScene(SceneProperties* properties);
	SceneProperties* getNextScene();
	FragmentData* requestFragment();
	bool hasFragments();
	void returnFragment(FragmentData* fd);
};

class RendererCommunicator : public Runnable {
	Socket socket;
	FragmentData* fd;
	int renderingSceneID;
	SceneProperties* renderingScene;
	Uint32 startTime;
  
	void returnFragmentIfRequested();
	void handleFragmentRequest();
	void handleSceneRequest();
	void handleFragmentReceive();
	void handleQuit();
	void checkForTimeout();
  
public:
	RendererCommunicator(Socket socket);
	virtual void run();
};  

class FragmentsManager : public Runnable {
	bool processFragment(FragmentData* fd);
  
public:
	virtual void run();
};

class FinishedCommunicatorWatcher : public Runnable {
public:
	virtual void run();
};

class RenderServer : public Runnable {
public:
	virtual void run();
};

void handleMouse(SDL_MouseButtonEvent *mev);

int parseServerArguments(int argc, char** argv);
int startServer(int argc, char** argv);
