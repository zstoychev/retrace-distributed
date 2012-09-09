#include <vector>
#include <queue>
#include <exception>
#include <cstring>
#include "sdl.h"
#include "scene.h"
#include "camera.h"
#include "tcpsocket.h"
#include "cxxptl_sdl.h"
#include "constants.h"
#include "SDL/SDL.h"
#include "SDL/SDL_net.h"
#include "renderserver.h"

bool showUpdates = true;
extern int port;
Uint32 fragmentRenderingTimeout = 60000;
bool outputFiles = false;
const char* outputDir = ".";
std::queue<std::string> scenesToRender;

volatile bool finishRendering = false;
volatile bool disablePreview = false;
volatile bool inInteractive = false;
Event interactiveRequestChangeFrameEvent;
Event interactiveFrameChangedEvent;

std::queue<FragmentData*> finishedFragments;
Mutex finishedFragmentsMutex;
Condition finishedFragmentsCondition(finishedFragmentsMutex);

class RendererCommunicator;
std::queue<RendererCommunicator*> finishedCommunicators;
Mutex finishedCommunicatorsMutex;
Condition finishedCommunicatorsCondition(finishedCommunicatorsMutex);

SceneManager sceneManager;
Mutex sceneFragmentsMutex;
Condition sceneFragmentsCondition(sceneFragmentsMutex);

template <typename T> void appendToQueue(std::queue<T>& destination, std::queue<T>& source) {
	while (!source.empty()) {
		destination.push(source.front());
		source.pop();
	}
}

void setFinishRendering() {
	finishRendering = true;
	sceneFragmentsCondition.signalAll();
	finishedFragmentsCondition.signal();
	interactiveRequestChangeFrameEvent.signal();
}

void SceneManager::addFragments(SceneProperties* properties, std::queue<FragmentData*>* sceneFragmentsQueue) {
	properties->inAA = false;

	std::vector<Rect> buckets = getBucketsList(properties->scene.settings.frameWidth, properties->scene.settings.frameHeight);
	properties->finishedFragments = 0;
	properties->totalFragments = (int) buckets.size();
		
	for (size_t i = 0; i < buckets.size(); i++) {
		sceneFragmentsQueue->push(new FragmentData(properties, buckets[i]));
	}
}

void SceneManager::startRenderingNextScene() {
	if (!queuedScenes.empty() && !disableAddingScenes) {
		SceneProperties* properties = new SceneProperties;
		properties->id = totalScenes++;
		std::string filename = queuedScenes.front();
		queuedScenes.pop();

		if (!properties->scene.parseScene(filename.c_str(), &properties->content)) {
			printf("Could not parse the scene %s!\n", filename.c_str());
			return;
		}

		if (properties->scene.settings.interactive) {
			disableAddingScenes = true;
		}

		std::queue<FragmentData*>* sceneFragmentsQueue = new std::queue<FragmentData*>();
		addFragments(properties, sceneFragmentsQueue);

		{
			Lock lock(sceneFragmentsMutex);
			int index = renderingScenesCount;
			sceneFragments[index] = sceneFragmentsQueue;
			renderingScenes[index] = properties;
			renderingScenesCount++;

			if (index == 0) {
				setCurrentScene();
			}
		}

		sceneFragmentsCondition.signalAll();
	}
}

int SceneManager::getSceneIndex(SceneProperties* properties) {
	for (int index = 0; index < renderingScenesCount; index++) {
		if (renderingScenes[index] == properties) {
			return index;
		}
	}

	return -1;
}
	
void SceneManager::setCurrentScene() {
	::setCurrentScene(renderingScenes[0]);
}

SceneManager::SceneManager()
	: renderingScenesCount(0), totalScenes(0), disableAddingScenes(false) {
}

SceneManager::~SceneManager() {
	for (int i = 0; i < renderingScenesCount; i++) {
		while (!sceneFragments[i]->empty()) {
			delete sceneFragments[i]->front();
			sceneFragments[i]->pop();
		}
		delete sceneFragments[i];
		delete renderingScenes[i];
	}
}

void SceneManager::addSceneFile(const char* fileName) {
	queuedScenes.push(fileName);

	if (renderingScenesCount < MAX_SIMULTANEOUS_SCENES) {
		startRenderingNextScene();
	}
}

void SceneManager::rerenderScene(SceneProperties* properties) {
	std::queue<FragmentData*>* sceneFragmentsQueue = new std::queue<FragmentData*>();
	addFragments(properties, sceneFragmentsQueue);

	{
		Lock lock(sceneFragmentsMutex);
		int index = getSceneIndex(properties);
		delete sceneFragments[index];
		sceneFragments[index] = sceneFragmentsQueue;
		properties->id = totalScenes++;
	}

	sceneFragmentsCondition.signalAll();
}

void SceneManager::addAAFragments(SceneProperties* properties) {
	properties->inAA = true;

	std::vector<Rect> buckets = getBucketsList(properties->scene.settings.frameWidth, properties->scene.settings.frameHeight);
	properties->finishedFragments = 0;
	properties->totalFragments = (int) buckets.size();

	std::queue<FragmentData*> fragments;
	std::queue<FragmentData*> readyFragments;

	for (size_t i = 0; i < buckets.size(); i++) {
		FragmentData* fd = new FragmentData(properties, buckets[i]);
		const Rect& r = fd->fragment;

		bool foundPixelWithAANeed = false;

		for (int y = r.y0; y < r.y1; y++) {
			for (int x = r.x0; x < r.x1; x++) {
				fd->data[y - r.y0][x - r.x0] = fd->sceneProperties->vfb[y][x];

				Color neighs[4];
				neighs[0] = y > 0 ? properties->vfb[y-1][x] : properties->vfb[y][x];
				neighs[1] = y < frameHeight() - 1 ? properties->vfb[y+1][x] : properties->vfb[y][x];
				neighs[2] = x > 0 ? properties->vfb[y][x-1] : properties->vfb[y][x];
				neighs[3] = x < frameWidth() - 1 ? properties->vfb[y][x+1] : properties->vfb[y][x];
				Color average = (properties->vfb[y][x] + neighs[0] + neighs[1] + neighs[2] + neighs[3]) / 5;
				for (int i = 0; i < 4; i++) {
					if (colorDifference(neighs[i], average) > properties->scene.settings.aaThresh) {
						fd->needsAA.set(y - r.y0, x - r.x0, true);
						foundPixelWithAANeed = true;
					}
				}
			}
		}

		if (foundPixelWithAANeed) {
			fragments.push(fd);
		} else {
			readyFragments.push(fd);
		}
	}
		
	{
		Lock lock(sceneFragmentsMutex);
		int index = getSceneIndex(properties);
		appendToQueue(*sceneFragments[index], fragments);
	}
	sceneFragmentsCondition.signalAll();

	{
		Lock lock(finishedFragmentsMutex);
		appendToQueue(finishedFragments, readyFragments);
	}
}

void SceneManager::stopRenderingScene(SceneProperties* properties) {		
	std::queue<FragmentData*>* sceneFragmentsQueue;

	{
		Lock lock(sceneFragmentsMutex);
		int index = getSceneIndex(properties);
		sceneFragmentsQueue = sceneFragments[index];

		for (int i = index + 1; i < renderingScenesCount; i++) {
			sceneFragments[i - 1] = sceneFragments[i];
			renderingScenes[i - 1] = renderingScenes[i];
		}

		renderingScenesCount--;

		if (index == 0 && renderingScenesCount > 0) {
			setCurrentScene();
		}

		if (getNextScene() == 0) {
			disableAddingScenes = false;
		}
	}

	delete sceneFragmentsQueue;
	delete properties;

	startRenderingNextScene();
}

SceneProperties* SceneManager::getNextScene() {
	if (renderingScenesCount == 0) {
		return 0;
	} else {
		return renderingScenes[0];
	}
}

FragmentData* SceneManager::requestFragment() {
	FragmentData* result = 0;

	for (int i = 0; i < renderingScenesCount; i++) {
		if (!sceneFragments[i]->empty()) {
			result = sceneFragments[i]->front();
			sceneFragments[i]->pop();
			break;
		}
	}

	return result;
}

bool SceneManager::hasFragments() {
	for (int i = 0; i < renderingScenesCount; i++) {
		if (!sceneFragments[i]->empty()) {
			return true;
		}
	}

	return false;
}

void SceneManager::returnFragment(FragmentData* fd) {
	int index = getSceneIndex(fd->sceneProperties);
	sceneFragments[index]->push(fd);
}

void RendererCommunicator::returnFragmentIfRequested() {
	if (fd != 0) {
		{
			Lock lock(sceneFragmentsMutex);
			sceneManager.returnFragment(fd);
		}
		sceneFragmentsCondition.signal();

		fd = 0;
	}
}

void RendererCommunicator::handleFragmentRequest() {
	if (fd != 0) {
		throw UnexpectedBehaviourException();
	}

	{
		Lock lock(sceneFragmentsMutex);
		while (!sceneManager.hasFragments()) {
			if (finishRendering) {
				return;
			}
			sceneFragmentsCondition.wait();
		}

		fd = sceneManager.requestFragment();
	}

	if (fd->sceneProperties->id != renderingSceneID) {
		renderingSceneID = fd->sceneProperties->id;
		renderingScene = fd->sceneProperties;

		socket.writeLine("CHGSCENE");
		socket.write(&renderingScene->id, sizeof(renderingScene->id));
	}

	if (!fd->sceneProperties->inAA) {
		socket.writeLine("RENDERNOAA");
	} else {
		socket.writeLine("RENDERAA");
	}

	socket.write(&fd->fragment, sizeof(fd->fragment));
	if (fd->sceneProperties->inAA) {
		socket.write(&fd->needsAA, sizeof(fd->needsAA));

		if (!disablePreview && showUpdates && currentSP->id == fd->sceneProperties->id
			&& !fd->sceneProperties->scene.settings.interactive) {
			markRegion(fd->fragment);
		}
	}

	startTime = SDL_GetTicks();
}

void RendererCommunicator::handleSceneRequest() {
	if (renderingScene == 0) {
		throw UnexpectedBehaviourException();
	}

	socket.writeLine("RCVSCENE");
	int length = (int) renderingScene->content.size();
	socket.write(&length, sizeof(length));
	socket.write(renderingScene->content.c_str(), length);
}
	
void RendererCommunicator::handleFragmentReceive() {
	try {
		Color data[MAX_BUCKET_SIZE][MAX_BUCKET_SIZE];

		bool timeouted = !socket.read(data, sizeof(data), 1000);
			
		if (timeouted) {
			throw NetworkException();
		}

		if (fd == 0 ) {
			return;
		}

		if (fd->sceneProperties->inAA) {
			for (int i = 0; i < MAX_BUCKET_SIZE; i++) {
				for (int j = 0; j < MAX_BUCKET_SIZE; j++) {
					if (fd->needsAA.get(i, j)) {
						data[i][j] = (fd->data[i][j] + data[i][j]) / 5;
					} else {
						data[i][j] = fd->data[i][j];
					}
				}
			}
		}

		memcpy(fd->data,data, sizeof(data));

		{
			Lock lock(finishedFragmentsMutex);
			finishedFragments.push(fd);
		}
		finishedFragmentsCondition.signal();

		fd = 0;
	} catch(NetworkException&) {
		returnFragmentIfRequested();
		throw;
	}
}

void RendererCommunicator::handleQuit() {
	returnFragmentIfRequested();
	printf("QUIT request from %s\n", socket.getSocketAddressAsText().c_str());
}

void RendererCommunicator::checkForTimeout() {
	if ((SDL_GetTicks() - startTime) > fragmentRenderingTimeout) {
		returnFragmentIfRequested();
	}
}

RendererCommunicator::RendererCommunicator(Socket socket)
	: socket(socket) {
}

void RendererCommunicator::run() {
	try {
		fd = 0;
		renderingSceneID = -1;
		renderingScene = 0;

		while (!finishRendering) {
			char message[16];
			if (socket.readLine(message, sizeof(message), 100)) {
				if (strcmp(message, "RQSCENE") == 0) {
					handleSceneRequest();
				} else if (strcmp(message, "RQFRAGMENT") == 0) {
					handleFragmentRequest();
				} else if (strcmp(message, "RCVFRAGMENT") == 0) {
					handleFragmentReceive();
				} else if (strcmp(message, "QUIT") == 0) {
					handleQuit();
					break;
				} else {
					throw UnexpectedBehaviourException();
				}
			} else {
				checkForTimeout();
			}
		}			
	} catch(UnexpectedBehaviourException&) {
		printf("Unexpected behaviour from %s\n", socket.getSocketAddressAsText().c_str());
	} catch (NetworkException&) {
		printf("Network exception in the connection to %s\n", socket.getSocketAddressAsText().c_str());
	}

	printf("Closing the connection to %s\n", socket.getSocketAddressAsText().c_str());
	try {
		socket.writeLine("QUIT");
	} catch (NetworkException&) {
		printf("Network exception in the connection to %s\n", socket.getSocketAddressAsText().c_str());
	}
	socket.close();

	returnFragmentIfRequested();

	{
		Lock lock(finishedCommunicatorsMutex);
		finishedCommunicators.push(this);
	}
	finishedCommunicatorsCondition.signal();
}

bool FragmentsManager::processFragment(FragmentData* fd) {
	for (int y = fd->fragment.y0; y < fd->fragment.y1; y++) {
		for (int x = fd->fragment.x0; x < fd->fragment.x1; x++) {
			fd->sceneProperties->vfb[y][x] = fd->data[y - fd->fragment.y0][x - fd->fragment.x0];
		}
	}

	if (!disablePreview && showUpdates && !fd->sceneProperties->scene.settings.interactive
		&& fd->sceneProperties->id == currentSP->id) {
		if (!displayVFBRect(fd->fragment, fd->sceneProperties->vfb)) {
			setFinishRendering();
			return false;
		}
	}

	return true;
}

void FragmentsManager::run() {
	Uint32 start = SDL_GetTicks();

	while (!scenesToRender.empty()) {
		sceneManager.addSceneFile(scenesToRender.front().c_str());
		scenesToRender.pop();
	}

	bool fragmentsAdded = true;
	{
		Lock lock(sceneFragmentsMutex);
		if (!sceneManager.hasFragments()) {
			fragmentsAdded = false;
		}
	}
	if (!fragmentsAdded) {
		setFinishRendering();
		extern volatile bool wantToQuit;
		wantToQuit = true;
		return;
	}

	while (!finishRendering) {
		FragmentData* fd;

		{
			Lock lock(finishedFragmentsMutex);
			while (finishedFragments.empty()) {
				if (finishRendering) {
					return;
				}
				finishedFragmentsCondition.wait();
			}

			fd = finishedFragments.front();
			finishedFragments.pop();
		}

		if (!processFragment(fd)) {
			return;
		}

		SceneProperties* sp = fd->sceneProperties;
		delete fd;

		sp->finishedFragments++;

		if (sp->isReady()) {
			Scene& scene = sp->scene;
			if (!sp->inAA && scene.settings.wantAA
				&& !scene.camera->dof && !scene.settings.gi) {
				sceneManager.addAAFragments(sp);
			} else {
				if(!disablePreview && currentSP->id == sp->id && (!showUpdates || scene.settings.interactive)) {
					if (!displayVFB(sp->vfb)) {
						setFinishRendering();
						return;
					}
				}

				if (scene.settings.interactive && inInteractive) {
					interactiveRequestChangeFrameEvent.signal();
					interactiveFrameChangedEvent.wait();

					if (inInteractive) {
						sp->content = scene.getAsText();
						sceneManager.rerenderScene(sp);
					} else {
						scene.settings.interactive = false;
						if (!displayVFB(sp->vfb)) {
							setFinishRendering();
							return;
						}
					}
				}
				
				if (disablePreview || !scene.settings.interactive) {
					if (outputFiles) {
						char filenamePrefix[256];
						sprintf(filenamePrefix, "%s/scene%d_", outputDir, sp->id);
						takeScreenshot(sp, filenamePrefix);
					}

					sceneManager.stopRenderingScene(sp);

					SceneProperties* nextScene;
					{
						Lock lock(sceneFragmentsMutex);
						nextScene = sceneManager.getNextScene();
					}

					if (nextScene == 0) {
						extern double lastRenderTime;
						lastRenderTime = (SDL_GetTicks() - start) / 1000.0;
						printf("Render time: %0.2lf seconds\n", lastRenderTime);
						finishRendering = true;

						sceneFragmentsCondition.signalAll();
					} else {
						if (!disablePreview && showUpdates && !nextScene->scene.settings.interactive) {
							if (!displayVFB(nextScene->vfb)) {
								setFinishRendering();
								return;
							}
						}
					}
				}
			}
		}
	}
}

void FinishedCommunicatorWatcher::run() {
	while (!finishRendering) {
		Lock lock(finishedCommunicatorsMutex);
		while (finishedCommunicators.empty()) {
			if (finishRendering) {
				return;
			}
			finishedCommunicatorsCondition.wait();
		}
		delete finishedCommunicators.front();
		finishedCommunicators.pop();
	}
}

void RenderServer::run() {
	SimpleCachedThreadPool<42> pool;

	FragmentsManager fm;
	FinishedCommunicatorWatcher fcw;

	pool.run_async(&fm);
	pool.run_async(&fcw);

	try {
		ServerSocket server(port);
		printf("Server started.\n");
		
		while (!finishRendering) {
			Socket client;
			bool established = server.accept(client);
			if (established) {
				printf("Received connection from %s\n", client.getSocketAddressAsText().c_str());
				RendererCommunicator* rc = new RendererCommunicator(client);
				pool.run_async(rc);
			} else {
				SDL_Delay(100);
			}
		}
		server.close();
	} catch(NetworkException&) {
		printf("A network error has occured! Stopping server...");
		finishRendering = true;
		sceneFragmentsCondition.signalAll();
		finishedFragmentsCondition.signalAll();
	}

	finishedCommunicatorsCondition.signalAll();

	pool.wait();

	while (!finishedCommunicators.empty()) {
		delete finishedCommunicators.front();
		finishedCommunicators.pop();
	}
	while (!finishedFragments.empty()) {
		delete finishedFragments.front();
		finishedFragments.pop();
	}
}

void handleMouse(SDL_MouseButtonEvent *mev)
{
	
}

extern void printUsage();

int parseServerArguments(int argc, char** argv) {
	if (argc < 3) {
		printUsage();
		return 1;
	}
	if (strcmp(argv[1], "-s") == 0) {
		scenesToRender.push(argv[2]);
	} else if (strcmp(argv[1], "-ss") == 0) {
		FILE* scenesFile = fopen(argv[2], "r");
		if (!scenesFile) {
			printf("Cannot open file %s!", argv[2]);
			return 1;
		}

		char filename[256];
		while (fgets(filename, sizeof(filename), scenesFile) != 0) {
			int length = strlen(filename);
			if (filename[length - 1] == '\n') {
				filename[length - 1] = '\0';
			}
			scenesToRender.push(filename);
		}

		fclose(scenesFile);
	} else {
		printUsage();
		return 1;
	}
	
	for (int i = 3; i < argc; i++) {
		if (strcmp(argv[i], "-disablePreview") == 0) {
			disablePreview = true;
		} else if (strcmp(argv[i], "-disableUpdates") == 0) {
			showUpdates = false;
		} else if (strcmp(argv[i], "-port") == 0 && i < argc - 1) {
			port = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-fragmentRenderingTimeout") == 0 && i < argc - 1) {
			fragmentRenderingTimeout = atoi(argv[++i]);
		} else if (strcmp(argv[i], "-outputDir") == 0 && i < argc - 1) {
			outputFiles = true;
			outputDir = argv[++i];
		} else {
			printUsage();
			return 1;
		}
	}

	return 0;
}

int startServer(int argc, char** argv) {
	int result = parseServerArguments(argc, argv);
	if (result != 0) {
		return result;
	}

	if (!initializeSDL())  return 1;
	if (!disablePreview && !initGraphics()) {
		return 1;
	}

	disableAcceleratedStructures = true;

	SimpleCachedThreadPool<1> pool;
	RenderServer server;
	pool.run_async(&server);
	
	if (!disablePreview) {
		handleEvents();
		waitForUserExit();
		closeGraphics();
	}

	pool.wait();
	
	return 0;
}
