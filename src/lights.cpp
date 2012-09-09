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


#include <string.h>
#include "lights.h"
#include "constants.h"
#include "matrix.h"
#include "random_generator.h"

Light::Light()
{
	col = Color(10000.0f, 10000.0f, 10000.0f);
}
Light::~Light()
{
}

void Light::fillProperties(ParsedBlock& pb)
{
	pb.getColorProp("color", &col);
}

/** @class PointLight */

PointLight::PointLight()
{
	pos = Vector(0.0, 0.0, 0.0);
}

void PointLight::fillProperties(ParsedBlock& pb)
{
	Light::fillProperties(pb);
	pb.getVectorProp("pos", &pos);
}

/** @class RectLight */

RectLight::RectLight()
{
	xSubd = ySubd = 2;
	T.resetTransform();
}

void RectLight::getNthSample(int sampleIdx, const Vector& hitPos, Vector& samplePos, Color& color)
{
	Random& rnd = getRandomGen();
	// stratified sampling: subdivide the unit square into xSubd * ySubd subrectangles
	// and take a random points from each one of them. I.e., getNthSample() gets a random sample
	// from one of those subrectangles. The subrectangle in question is determined by sampleIdx.
	float x = (sampleIdx % xSubd + rnd.randfloat()) / (float) xSubd;
	float y = (sampleIdx / xSubd + rnd.randfloat()) / (float) ySubd;
	
	// create the sample in local space - in the (-0.5, 0, -0.5) -- (0.5, 0, 0.5) square
	Vector sample(x - 0.5, 0, y - 0.5);
	samplePos = T.transformPoint(sample); // transform the sample into world space
	
	// transform the point that is to be shaded into local space first
	Vector hitPos_LS = IT.transformPoint(hitPos);
	if (hitPos_LS.y > 0) {
		color.makeZero(); // point is behind the light - unilluminated.
	} else {
		// otherwise, return light color, attenuated by the angle of incidence
		// (the cosine between the light's direction and the normed ray toward the hitpos)
		color = this->col * float((area * dot(Vector(0, -1, 0), hitPos_LS) / hitPos_LS.length()));
	}
}

double RectLight::intersect(const Ray& ray)
{
	Ray ray_LS = IT.transformRay(ray);
	// check if ray_LS (the incoming ray, transformed in local space) hits the oriented square 1x1, resting
	// at (0, 0, 0), pointing downwards:
	if (ray_LS.start.y >= 0) return +INF; // ray start is in the wrong subspace; no intersection is possible
	if (ray_LS.dir.y <= 0) return +INF; // ray direction points downwards; no intersection is possible
	double lengthToIntersection = -(ray_LS.start.y / ray_LS.dir.y); // intersect with XZ plane
	Vector p = ray_LS.start + ray_LS.dir * lengthToIntersection;
	if (fabs(p.x) < 0.5 && fabs(p.z) < 0.5) {
		// retransform the hit point to world space and return the true distance to intersection
		return (T.transformPoint(p) - ray.start).length(); // the hit point is inside the 1x1 square - return the intersection length
	} else
		return +INF;
}

void RectLight::fillProperties(ParsedBlock& pb)
{
	Light::fillProperties(pb);
	pb.getIntProp("xSubd", &xSubd, 1);
	pb.getIntProp("ySubd", &ySubd, 1);
	pb.getTransformProp(T, IT);
}

void RectLight::beginFrame(void)
{
	center = T.transformPoint(Vector(0, 0, 0));
	Vector a = T.transformPoint(Vector(-0.5, 0.0, -0.5));
	Vector b = T.transformPoint(Vector( 0.5, 0.0, -0.5));
	Vector c = T.transformPoint(Vector( 0.5, 0.0,  0.5));
	float width = (float) (b - a).length();
	float height = (float) (b - c).length();
	area = width * height; // obtain the area of the light, in world space
}

float RectLight::solidAngle(const Vector& point)
{
	Vector p = IT.transformPoint(point);
	if (p.y > 0) return 0; // if the point is behind the light - assume zero area, the point is not lit
	p.normalize();
	float cosA = (float) dot(Vector(0, -1, 0), p); // cosine between the direction to the point in question and the light's normal
	float dSqr = (float) (point - center).lengthSqr(); // quadratic attenuation
	return area * cosA / (1 + dSqr);
	/*     ^^^^   ^^^^   ^^^^^^^^^^
	 *the terms in the approximate formula above, explained:
	 * area - captures the light's size; bigger lights illuminate more brightly;
	 * cosA - captuers the fact, that a light that is rotated with respect to the illumination point and does not shine
	 *        directly towards it has a smaller spherical projection. In the limit case, where the two vectors are
	 *        orthogonal, the light's projection degenerates into a line, having zero area.
	 * (1 + dSqr) - captures the effect of attenuation due to distance. When the distance to the light is
	 *             increased twice, the lights's spherical projections becomes smaller, with each "side" decreasing in half.
	 *             Thus the projected area is four times smaller. The +1 term is to compensate for overexposure in points
	 *             that are too "near", and d becomes close to zero, thus the result as a whole becomes large, possibly
	 *             larger than 2*PI. Whereas, even that close, the maximum projection can only approach 2*PI, hence the
	 *             correction.
	 */
}
