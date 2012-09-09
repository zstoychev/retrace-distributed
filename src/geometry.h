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
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__


#include "vector.h"
#include "matrix.h"
#include "scene.h"

/// a structure, that holds all the info, which a Geometry::intersect() method
/// may need to save when an intersection is found.
class Geometry;
struct IntersectionInfo {
	Vector ip; //!< intersection point in the world-space
	double distance; //!< the distance to the intersection point along the ray
	Vector norm; //!< the normal of the geometry at the intersection point
	double u, v; //!< 2D UV coordinates for texturing, etc.
	const Geometry* g;
};


/// An abstract class that represents any intersectable primitive in the scene.
class Intersectable {
public:
	virtual ~Intersectable() {}
	/// Intersect a geometry with a ray. Returns true if an intersection is found,
	/// in which case the info structure is filled with details about the intersection.
	virtual bool intersect(const Ray& ray, IntersectionInfo& info) const = 0;
	/// Checks if the given point is "inside" the geometry, for whatever definition of
	/// inside is appropriate for the object. Returns a boolean value accordingly.
	virtual bool isInside(const Vector& p) const = 0;
};

/// An abstract class, that describes a geometry in the scene.
class Geometry: public SceneElement, public Intersectable {
public:
	virtual ~Geometry() {}
	
	ElementType getElementType() const { return ELEM_GEOMETRY; }
	
	virtual const char* getName() const = 0; //!< a virtual function, which returns the name of a geometry
};

class Shader;

/// A Node, which holds a geometry, linked to a shader.
class Node: public SceneElement, public Intersectable {
public:
	Node() { geometry = NULL; shader = NULL; T.resetTransform(); skip_shadow = false; active = false; }
	Geometry* geometry;
	Shader* shader;
	Transform T, invT;
	Matrix invTransposedM;
	bool skip_shadow;
	bool active;
	
	ElementType getElementType() const { return ELEM_NODE; }
	
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	void beginRender(void) { 
		invT = T.inverseTransform();
		invTransposedM = transposeMatrix(invT.m);
	}
	
	void fillProperties(ParsedBlock& pb)
	{
		skip_shadow = false;
		T.resetTransform();
		pb.getGeometryProp("geometry", &geometry);
		pb.getShaderProp("shader", &shader);
		pb.getBoolProp("skip_shadow", &skip_shadow);
		pb.getBoolProp("active", &active);
		pb.getTransformProp(T, invT);
	}

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Node");
		se.addProperty("geometry", geometry);
		se.addProperty("shader", (SceneElement*) shader);
		se.addProperty("skip_shadow", skip_shadow);
		se.addProperty("active", active);
		se.addProperty(T);
	}

	bool isInside(const Vector& p) const
	{
		return geometry->isInside(invT.transformPoint(p));
	}
};

/// A simple plane, parallel to the XZ plane (coinciding with XZ when y == 0)
class Plane: public Geometry {
	double y; /// the offset of the plane from the origin along the Y axis.
	float bumps;
	double lengthX, lengthZ; /// if nonzero, limit the plane to [lengthX x lengthZ]
public:
	Plane() { y = 0; lengthX = 0; lengthZ = 0;  bumps = 0;}
	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("y", &y);
		pb.getFloatProp("bumps", &bumps, 0, 100);
		pb.getDoubleProp("lengthX", &lengthX, 0, INF);
		pb.getDoubleProp("lengthZ", &lengthZ, 0, INF);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Plane");
		se.addProperty("y", y);
		se.addProperty("bumps", bumps);
		se.addProperty("lengthX", lengthX);
		se.addProperty("lengthZ", lengthZ);
	}
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	const char* getName() const { return "Plane"; }
	/// It has volume 0, so being "inside" the plane makes little sense.
	/// We could, of course, make it (y - p.y == 0), but we can't use it for anything meaningful in either case. 
	bool isInside(const Vector& p) const
	{
		return false;
	}
};

class Mesh;
/// A sphere
class Sphere: public Geometry {
	Vector O; /// center of the sphere
	double R; /// the sphere's radius
	enum UVGenType {
		UVGEN_SPHERICAL,
		UVGEN_CUBE,
	};
	UVGenType uvgen;
public:
	Sphere(Vector O = Vector(0, 0, 0), double R = 1): // canonic sphere
		O(O), R(R), uvgen(UVGEN_SPHERICAL) {}
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Sphere");
		se.addProperty("R", R);
		se.addProperty("O", O);
		se.addProperty("uvgen", uvgen == UVGEN_CUBE ? "cube" : "spherical");
	}
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	const char* getName() const { return "Sphere"; }
	bool isInside(const Vector& p) const
	{
		return (p - O).lengthSqr() < R*R;
	}
};

class Cube: public Geometry {
	Vector O; /// center of the cube
	double side; /// the cube's side
	Vector Vmin, Vmax; /// opposite corners of the cube
public:
	Cube() {
		O.makeZero(); side = 1;

		double halfSide = side/2;
		Vmin = O - Vector(halfSide, halfSide, halfSide);
		Vmax = O + Vector(halfSide, halfSide, halfSide);
	}
	void fillProperties(ParsedBlock& pb)
	{
		pb.getVectorProp("O", &O);
		pb.getDoubleProp("side", &side, 1e-9); // positive side
		double halfSide = side/2;
		Vmin = O - Vector(halfSide, halfSide, halfSide);
		Vmax = O + Vector(halfSide, halfSide, halfSide);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Cube");
		se.addProperty("O", O);
		se.addProperty("side", side);
	}
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	const char* getName() const { return "Cube"; }
	bool isInside(const Vector& p) const
	{
		return (p.x > Vmin.x && p.x < Vmax.x && p.y > Vmin.y && p.y < Vmax.y && p.z > Vmin.z && p.z < Vmax.z);
	}
};

// a cylinder, centered at (0, 0, 0), with height 2*h and radius r
class Cylinder: public Geometry {
	double h, r;
public:
	//
	Cylinder(double _h = 1.0, double _r = 1.0) { h = _h; r = _r; }
	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("h", &h, 1e-9); // must be positive
		pb.getDoubleProp("r", &r, 1e-9); // must be positive
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Cylinder");
		se.addProperty("h", h);
		se.addProperty("r", r);
	}
	bool isInside(const Vector& p) const
	{
		return (fabs(p.y) < h && (p.x * p.x + p.z * p.z) < r * r);
	}
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	const char* getName() const { return "Cylinder"; }
};

// a torus, centered at (0, 0, 0), radius R and inner (tube) radius r
class Torus: public Geometry {
	double R, r;
	Cylinder boundingCyl; // bounding cylinder for the torus
public:
	Torus() { R = 3; r = 1; boundingCyl = Cylinder(r + 1e-6, R + r + 1e-6); } 
	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("R", &R, 1e-9); // must be positive
		pb.getDoubleProp("r", &r, 1e-9); // must be positive
		boundingCyl = Cylinder(r + 1e-6, R + r + 1e-6);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Torus");
		se.addProperty("R", R);
		se.addProperty("r", r);
	}
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	const char* getName() const { return "Torus"; }
	bool isInside(const Vector& p) const;
};


class CsgOp: public Geometry {
protected:
	Intersectable* left, *right;
public:
	CsgOp() { left = right = NULL; }
	virtual bool boolOp(bool insideL, bool insideR) const = 0;
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.requiredProp("left");
		pb.requiredProp("right");
		pb.getIntersectableProp("left", &left);
		pb.getIntersectableProp("right", &right);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("CsgOp");
		se.addProperty("left", (SceneElement*) left);
		se.addProperty("right", (SceneElement*) right);
	}
	bool isInside(const Vector& p) const
	{
		return boolOp(left->isInside(p), right->isInside(p));
	}
};

class CsgUnion: public CsgOp {
public:
	bool boolOp(bool insideL, bool insideR) const { return insideL || insideR; }
	const char* getName() const { return "CsgUnion"; }
	virtual void fillSceneEntry(SceneEntry& se) const {
		CsgOp::fillSceneEntry(se);
		se.setType("CsgUnion");
	}
};

class CsgInter: public CsgOp {
public:
	bool boolOp(bool insideL, bool insideR) const { return insideL && insideR; }
	const char* getName() const { return "CsgInter"; }
	virtual void fillSceneEntry(SceneEntry& se) const {
		CsgOp::fillSceneEntry(se);
		se.setType("CsgInter");
	}
};

class CsgDiff: public CsgOp {
public:
	bool intersect(const Ray& ray, IntersectionInfo& info) const; // needed because we have to flip the normal in some cases.
	bool boolOp(bool insideL, bool insideR) const { return insideL && !insideR; }
	const char* getName() const { return "CsgDiff"; }
	virtual void fillSceneEntry(SceneEntry& se) const {
		CsgOp::fillSceneEntry(se);
		se.setType("CsgDiff");
	}
};

#endif // __GEOMETRY_H__
