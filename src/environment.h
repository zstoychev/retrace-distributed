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
#ifndef __ENVIRONMENT_H__
#define __ENVIRONMENT_H__

#include "color.h"
#include "vector.h"
#include "scene.h"

class Environment: public SceneElement {
protected:
	double yaw; // yaw angle (rotate the whole environment around axis Y with this angle)
	double cY, sY; // cos(yaw), sin(yaw)
public:
	Environment() { yaw = 0; }
	virtual ~Environment() {}
	/// gets a color from the environment at the specified direction
	virtual Color getEnvironment(const Vector& dir) = 0;
	
	ElementType getElementType() const { return ELEM_ENVIRONMENT; }
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("yaw", &yaw);
	}

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Environment");
		se.addProperty("yaw", yaw);
	}
	
	void beginFrame() { cY = cos(toRadians(yaw)); sY = sin(toRadians(yaw));	}
};

class Bitmap;
class CubemapEnvironment: public Environment {
	Bitmap* maps[6];
	Color getSide(const Bitmap& bmp, double x, double y);
	bool loadMaps(const char* folder);
	char folder[256];
public:
 	/// loads a cubemap from 6 separate images, from the specified folder.
 	/// The images have to be named "posX.bmp", "negX.bmp", "posY.bmp", ...
 	/// (or they may be .pfm images, not .bmp).
 	/// The folder specification should include a trailing slash; 
 	/// e.g. "/images/cubemaps/cathedral/" is OK.
	CubemapEnvironment();
	~CubemapEnvironment();
	Color getEnvironment(const Vector& dir);
	
	void fillProperties(ParsedBlock& pb);

	virtual void fillSceneEntry(SceneEntry& se) const {
		Environment::fillSceneEntry(se);
		se.setType("CubemapEnvironment");
		se.addFilenameProperty("folder", folder);
	}
};

/// This class implements a spherical environment, specified with a bitmap texture
class SphericalEnvironment: public Environment {
	Bitmap* texture;
	char textureFileName[256];
public:
	SphericalEnvironment() { texture = NULL; }
	~SphericalEnvironment();
	virtual Color getEnvironment(const Vector& dir);
	// from SceneElement:
	void fillProperties(ParsedBlock& pb);

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("SphericalEnvironment");
		Environment::fillSceneEntry(se);
		se.addFilenameProperty("textureFileName", textureFileName);
	}
};

#endif // __ENVIRONMENT_H__
