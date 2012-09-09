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
#ifndef __LIGHTS_H__
#define __LIGHTS_H__

#include "constants.h"
#include "scene.h"
#include "matrix.h"

class Light: public SceneElement {
public:
	Color col; //!< "base" color of the light
	
	Light();
	virtual ~Light();
	
	/// Gets the number of samples required for this light
	virtual int getNumSamples(void) = 0;
	
	/// Gets a sample light position
	/// @param sampleIdx - the index of the sample, in [0..getNumSamples)
	/// @param hitPos - input - the position, whose illumination we're considering.
	/// @param samplePos [out] - the output light sample position
	/// @param color [out] - the output light intensity (color). It might depend on the hitPos (e.g. if the hitPos is behind the lamp, the color would be 0)
	virtual void getNthSample(int sampleIdx, const Vector& hitPos, Vector& samplePos, Color& color) = 0;
	
	virtual double intersect(const Ray& ray) = 0; //!< Intersects a ray with the light. Returns the distance to intersection, or +INF if there is none.
	
	/// Returns the approximate solid angle of the light with respect to a particular point.
	/// I.e., if the light's contour gets projected to the unit sphere around this point, find out the
	/// area of the projection (the result is in [0..2*PI])
	virtual float solidAngle(const Vector& point) = 0;
	
	// From SceneElement:
	ElementType getElementType() const { return ELEM_LIGHT; }
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Light");
		se.addProperty("color", col);
	}
};

/// Represents the good old point light
class PointLight: public Light {
	Vector pos; //!< Point light's position
public:
	PointLight();
	int getNumSamples(void) { return 1; }
	void getNthSample(int sampleIdx, const Vector& hitPos, Vector& samplePos, Color& color)
	{
		color = col;
		samplePos = pos;
	}
	float solidAngle(const Vector& point) { return 0.0f; /* point light has no area */ }

	double intersect(const Ray& ray) { return +INF; /* you cannot intersect a point light */ }
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		Light::fillSceneEntry(se);
		se.setType("PointLight");
		se.addProperty("pos", pos);
	}
};

/**
 * @class RectLight
 * 
 * Represents a rectangle light. A "canonic" light is simply a 1x1 square, centered at (0, 0, 0),
 * looking downwards (points with y > 0 are not illuminated). A transform is provided
 * that places the light anywhere in the scene.
 */
class RectLight: public Light {
	int xSubd, ySubd; //!< sample density along the X and Y directions of the light (there will be xSubd*ySubd samples)
	Transform T, IT;  //!< Transform (to situate, orient and scale the canonical light somewhere in the scene)
	float area;       //!< Light's area (recalculated in beginFrame())
	Vector center;    //!< Light's center in world space (recalculated in beginFrame())
public:
	RectLight();
	int getNumSamples(void) { return xSubd * ySubd; }
	
	void getNthSample(int sampleIdx, const Vector& hitPos, Vector& samplePos, Color& color);
	
	double intersect(const Ray& ray);
	
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		Light::fillSceneEntry(se);
		se.setType("RectLight");
		se.addProperty("xSubd", xSubd);
		se.addProperty("ySubd", ySubd);
		se.addProperty(T);
	}
	void beginFrame(void);
	float solidAngle(const Vector& point);
};

#endif // __LIGHTS_H__
