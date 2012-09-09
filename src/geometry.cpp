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

#include <assert.h>
#include "geometry.h"
#include "util.h"
#include <stdio.h>
#include <string.h>


bool Plane::intersect(const Ray& ray, IntersectionInfo& info) const
{
	// intersect a ray with a XZ plane:
	// if the ray is pointing to the horizon, or "up", but the plane is below us,
	// of if the ray is pointing down, and the plane is above us, we have no intersection
	if ((ray.start.y > y && ray.dir.y > -1e-9) || (ray.start.y < y && ray.dir.y < 1e-9))
		return false;
	
	// calculate how much distance to the target height we need to cover from the start
	double toCover = this->y - ray.start.y;
	// calculate how much we should scale ray.dir so that we reach the target height (this may be negative) ...
	double scaling = toCover / ray.dir.y;
	if (scaling < 0) return false; // ... and if it is, then the intersection point is behind us; bail out
	
	// calculate the intersection
	Vector ip = ray.start + ray.dir * scaling;
	if (lengthX > 0 && fabs(ip.x) > lengthX * 0.5) return false;
	if (lengthZ > 0 && fabs(ip.z) > lengthZ * 0.5) return false;
	info.ip = ip;
	info.distance = scaling;
	if (bumps == 0) {
		info.norm = Vector(0, 1, 0);
	} else {
		const float freqX[3] = { 0.5f, 0.878123f, 1.35165f }, freqZ[3] = { 0.412123f, 1.0128738f, 1.187287f };
		const float anglesX[3] = { 0.15f, 0.67f, 2.2f }, anglesZ[3] = { 0.91f, 1.11123f, 2.67f };
		float nx = 0, ny = 0;
		float sx = (float) info.ip.x, sy = (float) info.ip.z;
		for (int i = 0; i < 3; i++) {
			nx += sin(i + freqX[i] * (sx * cos(anglesX[i]) + sy * sin(anglesX[i]))) * 0.1f * bumps;
			ny += sin(i + freqZ[i] * (sx * cos(anglesZ[i]) + sy * sin(anglesZ[i]))) * 0.1f * bumps;
		}
		float rtotal = 1.0f / (float) sqrt(1.0f + sqr(nx) + sqr(ny));
		info.norm = Vector(nx * rtotal, rtotal, ny * rtotal);
	}
	if (!(ray.flags & RF_SHADOW)) {
		info.u = info.ip.x;
		info.v = info.ip.z;
		info.g = this;
	}
	return true;
}

void Sphere::fillProperties(ParsedBlock& pb)
{
	pb.getDoubleProp("R", &R, 1e-9); // positive radius
	pb.getVectorProp("O", &O);
	char genType[100] = "-";
	uvgen = UVGEN_SPHERICAL;
	if (pb.getStringProp("uvgen", genType) && !strcmp(genType, "cube"))
		uvgen = UVGEN_CUBE;
}

bool Sphere::intersect(const Ray& ray, IntersectionInfo& info) const
{
	// compute the sphere intersection using a quadratic equation:
	Vector H = ray.start - O;
	double A = ray.dir.lengthSqr();
	double B = 2 * dot(H, ray.dir);
	double C = H.lengthSqr() - R*R;
	double Dscr = B*B - 4*A*C;
	if (Dscr < 0) return false; // no solutions to the quadratic equation - then we don't have an intersection.
	double x1, x2;
	x1 = (-B + sqrt(Dscr)) / (2*A);
	x2 = (-B - sqrt(Dscr)) / (2*A);
	double sol = x2; // get the closer of the two solutions...
	if (sol < 0) sol = x1; // ... but if it's behind us, opt for the other one
	if (sol < 0) return false; // ... still behind? Then the whole sphere is behind us - no intersection.
	
	info.distance = sol;
	info.ip = ray.start + ray.dir * sol;
	if (!(ray.flags & RF_SHADOW)) {
		info.norm = info.ip - O; // generate the normal by getting the direction from the center to the ip
		info.norm.normalize();
		info.g = this;
		if (uvgen == UVGEN_SPHERICAL) {
			info.u = (PI + atan2(info.ip.z - O.z, info.ip.x - O.x))/(2*PI);
			info.v = 1.0 - (PI/2 + asin((info.ip.y - O.y)/R)) / PI;
		} else {
			Vector v = info.norm;
			int maxDim = v.maxDimension();
			switch (maxDim) {
				case 0: v = Vector(signOf(v.x) * v.z, -v.y, v.x); break;
				case 1: v = Vector(v.z, signOf(v.y) * v.x, v.y); break;
				case 2: v = Vector(-signOf(v.z) * v.x, -v.y, v.z); break;
			}
			info.u = (v.x + 1) * 0.5;
			info.v = (v.y + 1) * 0.5;
		}
	}
	return true;
}

static inline void testIntersect(const Ray& ray, IntersectionInfo& info, const Vector& faceCenter, double c3, double start, double dir, const Vector& normal, double side)
{
	double toCover = c3 - start;
	if ((toCover > 0 && dir < 1e-9) || (toCover < 0 && dir > -1e-9)) return;
	double scaling = toCover / dir;
	if (scaling < 0 || info.distance < scaling) return;
	Vector ip = ray.start + ray.dir * scaling;
	/*
	double distanceFromCenter = fabs(faceCenter.x - ip.x);
	distanceFromCenter = max(distanceFromCenter, fabs(faceCenter.y - ip.y));
	distanceFromCenter = max(distanceFromCenter, fabs(faceCenter.z - ip.z));
	if (distanceFromCenter > side/2) return;
	*/
	// new - faster code:
	if (fabs(faceCenter.x - ip.x) > side*0.5) return;
	if (fabs(faceCenter.y - ip.y) > side*0.5) return;
	if (fabs(faceCenter.z - ip.z) > side*0.5) return;
	if (scaling < info.distance) {
		info.distance = scaling;
		info.ip = ip;
		if (!(ray.flags & RF_SHADOW)) {
			info.norm = normal;
			info.u = ip.x + ip.y;
			info.v = ip.z;
		}
	}
}

bool Cube::intersect(const Ray& ray, IntersectionInfo& info) const
{
	IntersectionInfo closest;
	closest.distance = INF;
	
	// check for intersection with the negative X and positive X sides
	if (ray.start.x < Vmax.x) testIntersect(ray, closest, Vector(Vmin.x, O.y, O.z), Vmin.x, ray.start.x, ray.dir.x, Vector(-1, 0, 0), side);
	if (ray.start.x > Vmin.x) testIntersect(ray, closest, Vector(Vmax.x, O.y, O.z), Vmax.x, ray.start.x, ray.dir.x, Vector(+1, 0, 0), side);
	
	// check for intersection with the negative Y and positive Y sides
	if (ray.start.y < Vmax.y) testIntersect(ray, closest, Vector(O.x, Vmin.y, O.z), Vmin.y, ray.start.y, ray.dir.y, Vector(0, -1, 0), side);
	if (ray.start.y > Vmin.y) testIntersect(ray, closest, Vector(O.x, Vmax.y, O.z), Vmax.y, ray.start.y, ray.dir.y, Vector(0, +1, 0), side);
	
	// check for intersection with the negative Z and positive Z sides
	if (ray.start.z < Vmax.z) testIntersect(ray, closest, Vector(O.x, O.y, Vmin.z), Vmin.z, ray.start.z, ray.dir.z, Vector(0, 0, -1), side);
	if (ray.start.z > Vmin.z) testIntersect(ray, closest, Vector(O.x, O.y, Vmax.z), Vmax.z, ray.start.z, ray.dir.z, Vector(0, 0, +1), side);
	
	if (closest.distance < INF) {
		if (ray.flags & RF_SHADOW) {
			info.distance = closest.distance;
			info.ip = closest.ip;
		} else {
			info = closest;
			info.g = this;
		}
		return true;
	} else return false;
}

bool CsgOp::intersect(const Ray& inray, IntersectionInfo& ret) const
{
	Ray ray = inray;
	IntersectionInfo li, ri;
	//
	// calculate the status of the starting point
	bool insideL = left->isInside(ray.start);
	bool insideR = right->isInside(ray.start);
	bool initState = boolOp(insideL, insideR);
	
	// compute the first two intersections
	bool intL = left->intersect(ray, li);
	bool intR = right->intersect(ray, ri);
	
	double totalDist = 0;
	while (intL || intR) {
		// check if the appropriate intersection is with the left geometry first
		bool leftFirst = (intL && !intR) || (intL && intR && li.distance < ri.distance);
		if (leftFirst) {
			totalDist += li.distance;
			ri.distance -= li.distance;
			insideL = !insideL;
		} else {
			totalDist += ri.distance;
			li.distance -= ri.distance;
			insideR = !insideR;
		}
		if (boolOp(insideL, insideR) != initState) {
			ret = leftFirst ? li : ri;
			ret.distance = totalDist;
			ret.g = this;
			return true;
		}
		if (leftFirst) {
			ray.start = li.ip + ray.dir * 1e-6;
			intL = left->intersect(ray, li);
			if (intL) li.distance += 1e-6;
		} else {
			ray.start = ri.ip + ray.dir * 1e-6;
			intR = right->intersect(ray, ri);
			if (intR) ri.distance += 1e-6;
		}
	}
	
	return false;
}

bool CsgDiff::intersect(const Ray& ray, IntersectionInfo& info) const
{
	if (!CsgOp::intersect(ray, info)) return false;
	/*
	 * Consider the following CsgDiff: a larger sphere with a smaller sphere somewhere on its side
	 * The result is the larger sphere, with some "eaten out" part. The question is:
	 * Where should the normals point, in the surface of the "eaten out" parts?
	 * These normals are generated by the smaller sphere, and point to the inside of the interior of
	 * the larger. They are obviously wrong.
	 *
	 * Solution: when we detect a situation like this, we flip the normals.
	 */
	 if ((ray.flags & RF_SHADOW) == 0 && right->isInside(info.ip - ray.dir * 1e-6) != right->isInside(info.ip + ray.dir * 1e-6))
		info.norm = -info.norm;
	 return true;
}

bool Node::intersect(const Ray& ray, IntersectionInfo& info) const
{
	// intersect a ray in its canonic space; if an intersection exists there,
	// transfer the info to world space
	if (geometry->intersect(invT.transformRay(ray), info)) {
		info.ip = T.transformPoint(info.ip);
		info.distance = (info.ip - ray.start).length();
		if (!(ray.flags & RF_SHADOW)) {
			info.norm *= invTransposedM; // apply the transpose of the inversed matrix m to the normal
			info.norm.normalize(); // if the transform contains scaling, normalization might be needed.
		}
		return true;
	}
	else return false;
	
}

bool Cylinder::intersect(const Ray& ray, IntersectionInfo& info) const
{
	const double EPS = 1e-6;
	double X0 = ray.start.x;
	double Z0 = ray.start.z;
	double C = sqr(X0) + sqr(Z0) - sqr(r);
	
	double minDist = INF;

	double X = ray.dir.x;
	double Z = ray.dir.z;
	double A = sqr(X) + sqr(Z);
	if (A > EPS) {
		double len2d = sqrt(A);
		X /= len2d;
		Z /= len2d;
		A = 1;
		double B2 = X * X0 + Z * Z0;
		double D4 = B2 * B2 - A * C;
		if (D4 >= 0) {
			double c = (-B2 - sqrt(D4));
			if (c < 0) c = (-B2 + sqrt(D4));
			if (c > 0) {
				double d = c / len2d;
				double y = d * ray.dir.y + ray.start.y;
				if (y >= -h && y <= +h && d < minDist) {
					minDist = d;
					info.ip = ray.start + d * ray.dir;
					info.distance = d;
					if (!(ray.flags & RF_SHADOW)) {
						info.norm = info.ip - (y * Vector(0, 1, 0));
						info.norm.normalize();
						info.u = info.ip.x;
						info.v = info.ip.z;
						info.g = this;
					}
				}
			}
		}
	}
	if (fabs(ray.dir.y) > 1e-6) {
		for (int top = -1; top <= 1; top += 2) {
			double ydiff = (h * top) - ray.start.y;
			if ((ydiff * ray.dir.y) <= 0.0) continue;
			double dist = ydiff / ray.dir.y;
			if (dist > 0 && dist < minDist) {
				Vector crossPoint = ray.start + dist * ray.dir;
				Vector temp = crossPoint;
				temp.y = 0;
				if (temp.lengthSqr() > r * r) continue;
				minDist = dist;
				info.ip = ray.start + dist * ray.dir;
				info.distance = dist;
				if (!(ray.flags & RF_SHADOW)) {
					info.norm = Vector(0, top, 0);
					info.u = info.ip.x;
					info.v = info.ip.z;
					info.g = this;
				}
			}
		}
	}
	return minDist < INF;
}

bool Torus::intersect(const Ray& ray, IntersectionInfo& info) const
{
	/*
	 * Ray-torus intersection overview:
	 *   Bound the torus with a cylinder of height 2*r and a radius of R + r.
	 *   Find the two intersections of the ray with that cylinder.
	 *   If they exist, partition the space between them to 100 small steps.
	 *   Trace the ray along these steps, and, at each step, calculate the torus equation.
	 *   The torus boundary is exactly where the equation is zero; if, at some point,
	 *   our current point is on the "negative" side, we've dug below the surface; use
	 *   binary search between the "previous" point and the current one, to find where the
	 *   equation becomes zero.
	 *
	 *   In fact, the method must work even when the ray starts from within the torus, in
	 *   which case we need to negate the torus equation.
	 */
	IntersectionInfo test;
	double close_dist, far_dist;
	bool negate;
	// first, check for the intersection point with the bounding cylinder:
	if (!boundingCyl.intersect(ray, test)) return false; // doesn't intersect the torus at all
	
	// check if the starting point is within the bounding cylinder, in which case no
	// second ray is to be traced
	if (boundingCyl.isInside(ray.start)) {
		close_dist = 0;
		far_dist = test.distance;
	} else {
		close_dist = test.distance;
		// find the other intersection point:
		Ray newRay = ray;
		newRay.start = test.ip + ray.dir * 1e-6;
		if (!boundingCyl.intersect(newRay, test)) return false; // shoudn't happen usually, but just in case
		far_dist = test.distance + close_dist + 1e-6;
	}
	
	#define TORUS(p) (sqr(R - sqrt(sqr(p.x) + sqr(p.z))) + sqr(p.y) - sqr(r))
	// check if our starting point is inside the torus already:
	Vector p = ray.start + ray.dir * close_dist;
	negate = (TORUS(p) < 0);
	
	// sample along the ray between distances close_dist and far_dist:
	double step = (far_dist - close_dist) * 0.01;
	Vector inc = ray.dir * step;
	for (double d = close_dist; d < far_dist + step; d += step) {
		double torus_equation = TORUS(p);
		if ((!negate && torus_equation < 0) || (negate && torus_equation > 0)) {
			// dug below the surface; now, use binary search to refine the intersection point
			double left  = d - step; // left is the distance to a point just before the intersection
			double right = d;        // right is the distance to a point just after the intersection
			double mid = (left + right) * 0.5;
			while (right - left > 1e-9) {
				mid = (left + right) * 0.5;
				p = ray.start + ray.dir * mid;
				torus_equation = TORUS(p);
				if ((!negate && torus_equation < 0) || (negate && torus_equation > 0)) {
					// midpoint is in torus; intersection is in (left..mid]
					right = mid;
				} else {
					// midpoint is outside the torus; intersection is in [mid..right)
					left = mid;
				}
			}
			info.ip = ray.start + ray.dir * mid;
			info.distance = mid;
			if (!(ray.flags & RF_SHADOW)) {
				Vector tubePoint = info.ip;
				tubePoint.y = 0;
				double length = sqrt(sqr(tubePoint.x) + sqr(tubePoint.z));
				tubePoint *= R / length;
				info.norm = info.ip - tubePoint;
				info.norm.normalize();
				info.u = atan2(info.ip.z, info.ip.x) * (R/r) / PI ;
				info.v = atan2(info.ip.y, length - R) / PI;
				info.g = this;
			}
			return true;
		}
		p += inc; // increment 0.01 * (far_dist - close_dist) along the ray.
	}
	return false;
}

bool Torus::isInside(const Vector& p) const
{
	return TORUS(p) < 0;
}
