/***************************************************************************
 *   Copyright (C) 2009-2012 by Veselin Georgiev, Slavomir Kaslev et al    *
 *   admin@raytracing-bg.net                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <SDL/SDL.h>
#include <SDL/SDL_net.h>
#include <stdio.h>
#include "sdl.h"
#include "bitmap.h"
#include "scene.h"
#include "geometry.h"
#include "util.h"

SDL_Surface* screen = NULL;
SDL_Thread *render_thread;
SDL_mutex *render_lock;
extern volatile bool disablePreview;
extern volatile bool finishRendering;
extern volatile bool inInteractive;
volatile bool requestResolutionChange = false;
extern Event interactiveRequestChangeFrameEvent;
extern Event interactiveFrameChangedEvent;
Event resolutionChangeEvent;
volatile bool wantToQuit = false;

class MutexRAII {
	SDL_mutex* mutex;
public:
	MutexRAII(SDL_mutex* _mutex)
	{
		mutex = _mutex;
		SDL_mutexP(mutex);
	}
	~MutexRAII()
	{
		SDL_mutexV(mutex);
	}
};

bool initializeSDL() {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("Cannot initialize SDL: %s\n", SDL_GetError());
		return false;
	}
	
	if(SDLNet_Init() < 0) {
		printf("Cannot initialize SDL_NET: %s\n", SDL_GetError());
		return false;
	}

	return true;
}

bool changeResolution(int frameWidth, int frameHeight, bool fullscreen) {
	screen = SDL_SetVideoMode(frameWidth, frameHeight, 32, fullscreen ? SDL_FULLSCREEN : 0);
	if (!screen) {
		printf("Cannot set video mode %dx%d - %s\n", frameWidth, frameHeight, SDL_GetError());
		return false;
	}

	return true;
}

void changeResolutionIfNeeded() {
	{
		MutexRAII raii(render_lock);
		if (!disablePreview && frameWidth() != currentScene->settings.frameWidth || frameHeight() != currentScene->settings.frameHeight
			|| isInFullscreen() != currentScene->settings.interactive) {
			requestResolutionChange = true;
		}
	}

	if (requestResolutionChange) {
		resolutionChangeEvent.wait();
	}
}

/// try to create a frame window with the given dimensions
bool initGraphics()
{
	render_lock = SDL_CreateMutex();
	SDL_WM_SetCaption("retrace - rendering", NULL);
	return true;
}

/// closes SDL graphics
void closeGraphics(void)
{
	SDL_DestroyMutex(render_lock);
	SDL_Quit();
}

void setCurrentScene(SceneProperties* sp) {
	MutexRAII raii(render_lock);
	currentScene = &(sp->scene);
	currentSP = sp;

	if (currentScene->settings.interactive) {
		currentScene->beginRender();
		currentScene->beginFrame();

		extern Node* activeNode;
		extern int activeNodeIndex;

		activeNodeIndex = currentScene->nodes.empty() ? -1 : 0;
		activeNode = activeNodeIndex == -1 ? 0 : currentScene->nodes[activeNodeIndex];
		if (activeNode != 0) {
			activeNode->active = true;
		}
	}
}

/// displays a VFB (virtual frame buffer) to the real framebuffer, with the necessary color clipping
bool displayVFB(Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE])
{
	changeResolutionIfNeeded();

	{
		MutexRAII raii(render_lock);

		if (wantToQuit) return false;

		int rs = screen->format->Rshift;
		int gs = screen->format->Gshift;
		int bs = screen->format->Bshift;
		for (int y = 0; y < screen->h; y++) {
			Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
			for (int x = 0; x < screen->w; x++)
				row[x] = vfb[y][x].toRGB32(rs, gs, bs, currentScene->settings.outputGamma);
		}
		SDL_Flip(screen);
	}

	return true;
}

static void handleEvent(SDL_Event& ev)
{
	switch (ev.type) {
		case SDL_QUIT:
			wantToQuit = true;
			return;
		case SDL_KEYDOWN:
		{
			switch (ev.key.keysym.sym) {
				case SDLK_ESCAPE:
					wantToQuit = true;
					return;
				case SDLK_F12:
					takeScreenshot(currentSP + 0);
					break;
				default:
					break;
			}
			break;
		}
		case SDL_MOUSEBUTTONUP:
		{
			extern void handleMouse(SDL_MouseButtonEvent *mev);
			handleMouse(&ev.button);
			break;
		}
		default:
			break;
	}
}

/// waits the user to indicate he wants to close the application (by either clicking on the "X" of the window,
/// or by pressing ESC)
void waitForUserExit(void)
{
	char message[128];
	extern double lastRenderTime;
	sprintf(message, "retrace (rendertime: %.2lfs)", lastRenderTime);
	SDL_WM_SetCaption(message, NULL);
	SDL_Event ev;
	while (!wantToQuit) {
		while (!wantToQuit && SDL_WaitEvent(&ev)) {
			handleEvent(ev);
		}
	}
}

/// returns the frame width
int frameWidth(void)
{
	if (screen) return screen->w;
	return 0;
}

/// returns the frame height
int frameHeight(void)
{
	if (screen) return screen->h;
	return 0;
}

bool isInFullscreen() {
	if (screen) return (screen->flags & SDL_FULLSCREEN) != 0;
	return false;
}

bool takeScreenshot(SceneProperties* sp, const char* filenamePrefix)
{
	char fn[256];
	if (filenamePrefix == 0) {
		findEmptyFN(fn, "retrace_", ".bmp");
	} else {
		findEmptyFN(fn, filenamePrefix, ".bmp");
	}
	Bitmap bmp;
	bmp.generateEmptyImage(sp->scene.settings.frameWidth, sp->scene.settings.frameHeight);
	for (int y = 0; y < sp->scene.settings.frameHeight; y++)
		for (int x = 0; x < sp->scene.settings.frameWidth; x++)
			bmp.setPixel(x, y, sp->vfb[y][x]);
	bool res = bmp.saveBMP(fn, sp->scene.settings.outputGamma);
	if (res) printf("Rendered scene #%d saved to '%s'\n", sp->id, fn);
	else printf("Couldn't save rendered scene #%d to '%s'\n", sp->id, fn);
	return res;
}

void Rect::clip(int W, int H)
{
	x1 = min(x1, W);
	y1 = min(y1, H);
	w = max(0, x1 - x0);
	h = max(0, y1 - y0);
}

bool drawRect(Rect r, const Color& c)
{
	MutexRAII raii(render_lock);
	
	if (wantToQuit) return false;
	
	r.clip(frameWidth(), frameHeight());
	int rs = screen->format->Rshift;
	int gs = screen->format->Gshift;
	int bs = screen->format->Bshift;
	Uint32 clr = c.toRGB32(rs, gs, bs, currentScene->settings.outputGamma);
	for (int y = r.y0; y < r.y1; y++) {
		Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
		for (int x = r.x0; x < r.x1; x++)
			row[x] = clr;
	}
	SDL_UpdateRect(screen, r.x0, r.y0, r.w, r.h);
	
	return true;
}

bool displayVFBRect(Rect r, Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE])
{
	changeResolutionIfNeeded();

	{
		MutexRAII raii(render_lock);

		if (wantToQuit) return false;
	
		r.clip(frameWidth(), frameHeight());
		int rs = screen->format->Rshift;
		int gs = screen->format->Gshift;
		int bs = screen->format->Bshift;
		for (int y = r.y0; y < r.y1; y++) {
			Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
			for (int x = r.x0; x < r.x1; x++)
				row[x] = vfb[y][x].toRGB32(rs, gs, bs, currentScene->settings.outputGamma);
		}
		SDL_UpdateRect(screen, r.x0, r.y0, r.w, r.h);
	}
	
	return true;
}

bool markRegion(Rect r)
{
	MutexRAII raii(render_lock);

	if (wantToQuit) return false;
	
	r.clip(frameWidth(), frameHeight());
	const int L = 8;
	if (r.w < L || r.h < L) return true; // region is too small to be marked
	Uint32* row;
	const Uint32 BRACKET_COLOR = Color(0.0f, 1.0f, 0.0f).toRGB32();
	row = (Uint32*) ((Uint8*) screen->pixels + (r.y0) * screen->pitch);
	for (int x = r.x0; x < r.x0 + L; x++) row[x] = BRACKET_COLOR;
	for (int x = r.x1 - L - 1; x < r.x1; x++) row[x] = BRACKET_COLOR;
	row = (Uint32*) ((Uint8*) screen->pixels + (r.y1 - 1) * screen->pitch);
	for (int x = r.x0; x < r.x0 + L; x++) row[x] = BRACKET_COLOR;
	for (int x = r.x1 - L; x < r.x1; x++) row[x] = BRACKET_COLOR;
	for (int y = r.y0 + 1; y < r.y0 + L; y++) {
		row = (Uint32*) ((Uint8*) screen->pixels + (y) * screen->pitch);
		row[r.x0] = row[r.x1 - 1] = BRACKET_COLOR;
	}
	for (int y = r.y1 - L - 1; y < r.y1 - 1; y++) {
		row = (Uint32*) ((Uint8*) screen->pixels + (y) * screen->pitch);
		row[r.x0] = row[r.x1 - 1] = BRACKET_COLOR;
	}
	SDL_UpdateRect(screen, r.x0, r.y0, r.w, r.h);
	
	return true;
}

void handleEvents(void)
{
	Uint32 ticks;
	bool isInFirstFrame;

	while (!wantToQuit) {
		{
			MutexRAII raii(render_lock);
			if (finishRendering) break;
			if (requestResolutionChange) {
				if (!changeResolution(currentScene->settings.frameWidth, currentScene->settings.frameHeight, currentScene->settings.interactive)) {
					disablePreview = true;
				}
				inInteractive = currentScene->settings.interactive;

				if (inInteractive) {
					isInFirstFrame = true;
					SDL_ShowCursor(0);
				} else {
					SDL_ShowCursor(1);
				}

				resolutionChangeEvent.signal();

				requestResolutionChange = false;
			}
			if (!inInteractive) {
				SDL_Event ev;
				while (SDL_PollEvent(&ev)) {
					handleEvent(ev);
					if (wantToQuit) break;
				}
			}
		}

		if (!inInteractive) {
			SDL_Delay(100);
		} else {
			interactiveRequestChangeFrameEvent.wait();

			if (!isInFirstFrame) {
				MutexRAII raii(render_lock);
				currentScene->beginFrame();
				double renderTime = (SDL_GetTicks() - ticks) / 1000.0;
				extern void handleKbdMouse(bool& mustExit, double dt);
				bool mustExit = false;
				handleKbdMouse(mustExit, renderTime);
				if (mustExit) {
					inInteractive = false;
				}
			}
			isInFirstFrame = false;

			interactiveFrameChangedEvent.signal();

			ticks = SDL_GetTicks();
		}
	}

	while (!finishRendering) {
		SDL_Delay(100);
	}
}
