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
#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "vector.h"
#include "scene.h"

/// This is a base class, describing a camera.
class Camera: public SceneElement {
	// these internal vectors describe three of the ends of the imaginary
	// ray shooting screen
	Vector upLeft, upRight, downLeft;

	// the frontal, right and up directions after yaw-pitch-roll are applied.
	Vector upDir, rightDir, frontDir;
	double radius; // radius of the disc, where the dof rays originate from
public:
	Camera();
	Vector pos; //!< position of the camera in 3D.
	double yaw; //!< Yaw angle in degrees (rot. around the Y axis, meaningful values: [0..360])
	double pitch; //!< Pitch angle in degrees (rot. around the X axis, meaningful values: [-90..90])
	double roll; //!< Roll angle in degrees (rot. around the Z axis, meaningful values: [-180..180])
	double fov; //!< The Field of view in degrees (meaningful values: [3..160])
	double aspect; //!< The aspect ratio of the camera frame. Should usually be frameWidth/frameHeight,
	               /// based on the output resolution, but could be different on systems with non-square pixels.

	bool dof;         //!< Depth-of-Field enabled?
	bool autofocus;   //!< Autofocus on or off; if off, focusDist has to be specified manually
	double fNumber;   //!< The f-number of the lens
	int samples;      //!< How much rays to shoot per pixel
	double focusDist; //!< The distance to the focal plane in world units.
	double stereoSeparation; //!< If nonzero, a stereoscopic image is to be rendered and the param specifies
	                         ///  the distance between the two virtual cameras
	Color leftMask, rightMask; ///!< Stereoscopic mode: masks for the left- and right-eye image. Defaults: red/cyan.
	
	void beginFrame(void); //!< must be called before each frame. Computes the corner variables, needed for getScreenRay()
	void setPosition(Vector pos, Vector lookAt, Vector upDir);
	
	/// generates a screen ray through a pixel (x, y - screen coordinates, not necessarily integer).
	/// The offset parameter means:
	/// -1: the ray should originate from the "left" camera (in stereoscopic mode)
	///  0: normal ray, originates from the specified camera position
	/// +1: the ray should originate from the "right" camera (in stereoscopic mode)
	Ray getScreenRay(double x, double y, int offset = 0);
	Ray getCenterRay() const;
	
	ElementType getElementType() const { return ELEM_CAMERA; }
	
	void move(double ff, double ss);
	void rotate(double xx, double yy);

	void fillProperties(ParsedBlock& pb);

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Camera");
		se.addProperty("dof", dof);
		se.addProperty("autofocus", autofocus);
		se.addProperty("stereoSeparation", stereoSeparation);
		se.addProperty("leftMask", leftMask);
		se.addProperty("rightMask", rightMask);
		if (dof) {
			se.addProperty("fNumber", fNumber);
		}
		se.addProperty("samples", samples);
		if (dof && focusDist < INF) {
			se.addProperty("focusDist", focusDist);
		}
		se.addProperty("aspect", aspect);
		se.addProperty("pos", pos);
		se.addProperty("yaw", yaw);
		se.addProperty("pitch", pitch);
		se.addProperty("roll", roll);
	}
};










#endif // __CAMERA_H__
