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
#ifndef __SDL_H__
#define __SDL_H__

#include "color.h"
#include "constants.h"
#include "scene.h"
#include <vector>

bool initializeSDL();
bool changeResolution(int frameWidth, int frameHeight, bool fullscreen);
void changeResolutionIfNeeded();
bool initGraphics();
void closeGraphics(void);
void setCurrentScene(SceneProperties* sp);
bool displayVFB(Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE]); //!< displays the VFB (Virtual framebuffer) to the real one.
void waitForUserExit(void); //!< Pause. Wait until the user closes the application
int frameWidth(void); //!< returns the frame width (pixels)
int frameHeight(void); //!< returns the frame height (pixels)
bool isInFullscreen();

/// Takes a screen shot; converts the VFB to RGB32 and writes it to a .BMP file.
bool takeScreenshot(SceneProperties* sp, const char* filenamePrefix = NULL);

// fills a rectangle on the screen with a solid color
// fails if the render thread is about to be killed
bool drawRect(Rect r, const Color& c);

// same as displayVFB, but only updates a specific region.
// fails if the thread has to be killed
bool displayVFBRect(Rect r, Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE]);

// marks a region (places four temporary green corners)
// fails if the thread is to be killed
bool markRegion(Rect r);

void handleEvents();

extern volatile bool rendering; // used in main/worker thread synchronization

#endif // __SDL_H__
