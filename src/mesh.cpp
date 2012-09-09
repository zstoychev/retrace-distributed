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

#include <stdio.h>
#include <SDL/SDL.h>
#include "mesh.h"
#include "vector.h"

vector<string> tokenize(string s); // splits a string to tokens (the delimiter is assumed to be white space)
vector<string> split(string s, char separator); // same as tokenize, but using a fixed delimiter

static inline double getDouble(string s)
{
	if (s == "") return 0;
	double x;
	sscanf(s.c_str(), "%lf", &x);
	return x;
}

static inline int getInt(string s)
{
	if (s == "") return 0;
	int x;
	sscanf(s.c_str(), "%d", &x);
	return x;
}

// normalize a variable to be in the range [lim_min, lim_max)
//  if it's inside, just return it
//  otherwise, "wraparound" it so it fits the range
static inline double safeRange(double x, double lim_min, double lim_max)
{
	if (lim_min <= x && x < lim_max) return x;
	double f = (lim_max - lim_min);
	return lim_min + f * ((x - lim_min)/f - floor((x - lim_min)/f));
}

// parse a line like "1//2" "2//2" "3//2" to indices in the Triangle object
// this line comes from a subset of a "f" line in the .OBJ file.
Triangle::Triangle(string a, string b, string c)
{
	// a = "1",  "1/2", "1//2", "1/2/3"
	string things[3] = { a, b, c };
	for (int i = 0; i < 3; i++) {
		vector<string> parts = split(things[i], '/');
		t[i] = n[i] = 0;
		v[i] = getInt(parts[0]);
		if (parts.size() >= 2)
			t[i] = getInt(parts[1]);
		if (parts.size() >= 3)
			n[i] = getInt(parts[2]);
	}
}

// expresses C in terms of A and B, and returns the coefficients in p and q.
// i.e., find p, q such that C = p*A + q*B
static void solve2D(Vector A, Vector B, Vector C, double& p, double& q)
{
	double Dcr = A.x * B.y - A.y * B.x;
	if (fabs(Dcr) < 1e-9) {
		p = q = 0; // no sulution
		return;
	}
	// use Cramer's rule to determine p, q
	p = (C.x * B.y - C.y * B.x) / Dcr;
	q = (A.x * C.y - A.y * C.x) / Dcr;
}

bool Mesh::loadFromOBJ(const char* fn)
{
	Lock lock(cachedItemsMutex);
	
	if (filenameToObjFileProperties.find(fn) != filenameToObjFileProperties.end()) {
		ObjFileProperties*& properties = filenameToObjFileProperties[fn];
		vertices = properties->vertices;
		normals = properties->normals;
		uvs = properties->uvs;
		triangles = properties->triangles;
		bbox = properties->bbox;
		kdtreeroot = properties->kdtreeroot;
		hasNormals = properties->hasNormals;
	} else {
		FILE* f = fopen(fn, "rt");
		if (!f) return false;
		//
		// insert the dummies at index 0.
		// This is because all indices in the OBJ file are 1-based, and we want to have our
		// arrays 1-based as well for the following reasons:
		//  1) Consistency
		//  2) Simplification - normals/uvs may not exist, and having the Triangle's n[]/t[] arrays
		//     point to an existing (albeit zero) entry in the normals[]/uvs[] arrays can save us from
		//     some runtime checking
		vertices.push_back(Vector(0, 0, 0));
		normals.push_back(Vector(0, 0, 0));
		uvs.push_back(Vector(0, 0, 0));
		//
		char line[2048];
		double maxDist = 0;
		hasNormals = false;
		// read all lines of the input file
		while (fgets(line, sizeof(line), f)) {
			if (line[0] == '#') continue; // skip a comment line
			vector<string> tokens = tokenize(line);
			//
			if (tokens.size() < 3) continue; // skip lines having fewer than 3 tokens
			if (tokens[0] == "v") { // "v" line, e.g. "v 0.1111 -0.23123 5.150000"
				Vector v;
				v = Vector(getDouble(tokens[1]),
						   getDouble(tokens[2]),
						   getDouble(tokens[3]));
				maxDist = max(maxDist, v.length()); // update the maximum distance from the origin
				vertices.push_back(v);
			}
			if (tokens[0] == "vn") { // "vn" line, e.g. "vn 0 1 0"
				Vector vn;
				hasNormals = true;
				vn = Vector(getDouble(tokens[1]),
							getDouble(tokens[2]),
							getDouble(tokens[3]));
				normals.push_back(vn);
			}
			if (tokens[0] == "vt") { // "vt" line, e.g. "vt 0.5 0.125"
				Vector v;
				v = Vector(getDouble(tokens[1]),
						   getDouble(tokens[2]),
						   0);
				uvs.push_back(v);
			}
			if (tokens[0] == "f") { // "f" line, e.g. "f   1/1/1   2/2/1   3/3/1   4/4/1"
				int nTriangles = tokens.size() - 3;
				for (int i = 0; i < nTriangles; i++) {
					// enumerate all triangles in the input polygon.
					// The polygon has vertices 1, 2, ...N.
					// We create N - 2 triangles using vertices:
					//    (1, 2, 3); (1, 3, 4); (1, 4, 5); ... (1, i, i + 1) ... (1, N-1, N)
					Triangle t(tokens[1], tokens[2 + i], tokens[3 + i]);
					triangles.push_back(t);
				}
			}
		}
		//
		fclose(f);
		//
		// preprocess all triangles:
		for (int i = 0; i < (int) triangles.size(); i++) {
			Triangle& t = triangles[i];
			// calculate AB, AC and the geometric normal of the triangle
			t.AB = vertices[t.v[1]] - vertices[t.v[0]];
			t.AC = vertices[t.v[2]] - vertices[t.v[0]];
			t.gnormal = t.ABcrossAC = t.AB ^ t.AC;
			t.gnormal.normalize();
			// calculate dNdx, dNdy (for bump mapping)
			double px, py, qx, qy;
			Vector AB2D = uvs[t.t[1]] - uvs[t.t[0]]; // AB in texture space
			Vector AC2D = uvs[t.t[2]] - uvs[t.t[0]]; // AC in texture space
			solve2D(AB2D, AC2D, Vector(1, 0, 0), px, qx); // texture space's (1, 0) = px * AB2D + qx * AC2D
			solve2D(AB2D, AC2D, Vector(0, 1, 0), py, qy); // texture space's (0, 1) = py * AB2D + qy * AB2D
			t.dNdx = px * t.AB + qx * t.AC; // express the unit perturbances (dx and dy) in the 3D vertex space
			t.dNdy = py * t.AB + qy * t.AC; //  - we call them dNdx and dNdy
			t.dNdx.normalize(); // normalize the results, so the amount of bump doesn't depend on triangle size.
			t.dNdy.normalize();
		}
		//
		// calculate the bounding box of the mesh:
		bbox.makeEmpty();
		for (int i = 0; i < (int) triangles.size(); i++)
			for (int j = 0; j < 3; j++)
				bbox.add(vertices[triangles[i].v[j]]);
	}

	//
	// if the mesh is big enough and the relevant flag is set, build the K-d tree:
	if (!disableAcceleratedStructures && kdtreeroot == 0 && useKDTree && triangles.size() > 60) {
		buildTree();
	}

	if (!useKDTree) {
		kdtreeroot = 0;
	}

	if (filenameToObjFileProperties.find(fn) == filenameToObjFileProperties.end()) {
		ObjFileProperties* properties = new ObjFileProperties;
		properties->vertices = vertices;
		properties->normals = normals;
		properties->uvs = uvs;
		properties->triangles = triangles;
		properties->bbox = bbox;
		properties->kdtreeroot = kdtreeroot;
		properties->hasNormals = hasNormals;

		filenameToObjFileProperties[fn] = properties;
	} else if (useKDTree && filenameToObjFileProperties[fn]->kdtreeroot == 0) {
		filenameToObjFileProperties[fn]->kdtreeroot = kdtreeroot;
	}

	//
	// create the normals[] array - if needed:
	if (!hasNormals && autoSmooth) {
		hasNormals = true;
		normals.resize(vertices.size(), Vector(0, 0, 0)); // extend the normals[] array, and fill with zeros
		for (int i = 0; i < (int) triangles.size(); i++)
			for (int j = 0; j < 3; j++) {
				triangles[i].n[j] = triangles[i].v[j];
				normals[triangles[i].n[j]] += triangles[i].gnormal;
			}
		for (int i = 1; i < (int) normals.size(); i++)
			if (normals[i].lengthSqr() > 1e-9) normals[i].normalize();
	}

	return true;
}

bool intersectTriangleFast(const Ray& ray, const Vector& A, const Triangle& T, double& lambda2, double& lambda3, double& closestdist)
{
	const Vector& a = T.AB;
	const Vector& b = T.AC;
	Vector c = -ray.dir;
	Vector h = ray.start - A;
	/* 2. Solve the equation:
	 *
	 * A + lambda2 * AB + lambda3 * AC = ray.start + gamma * ray.dir
	 *
	 * which can be rearranged as:
	 * lambda2 * AB + lambda3 * AC - gamma * ray.dir = ray.start - A
	 *
	 * Which is a linear system of three rows and three unknowns, which we solve using Carmer's rule
	 */
	//
	// Find the determinant of the left part of the equation
	double Dcr = (T.ABcrossAC)*c;
	// check for zero; if it is zero, then the triangle and the ray are parallel
	if (fabs(Dcr) < 1e-9) return false;
	// find the reciprocal of the determinant. We would use this quantity later in order
	// to multiply by rDcr instead of divide by Dcr (division is much slower)
	double rDcr = 1.0 / Dcr;
	// calculate `gamma' by substituting the right part of the equation in the third column of the matrix,
	// getting the determinant, and dividing by Dcr)
	double gamma =   (T.ABcrossAC) * h * rDcr;
	// Is the intersection point behind us?  Is the intersection point worse than what we currently have?
	if (gamma <= 0 || gamma > closestdist) return false;
	lambda2 = (h ^ b) * c * rDcr;
	// Check if it is in range (barycentric coordinates)
	if (lambda2 < 0 || lambda2 > 1) return false;
	lambda3 = (a ^ h) * c * rDcr;
	
	// Calculate lambda3 and check if it is in range as well
	if (lambda3 < 0 || lambda3 > 1) return false;
	if (lambda2 + lambda3 > 1) return false;
	
	closestdist = gamma;
	return true;
}

void Mesh::fillInfo(const Ray& ray, IntersectionInfo& info, int triangleIndex, double lambda2, double lambda3, double gamma) const
{
	info.distance = gamma;
	info.ip = ray.start + ray.dir * gamma;
	if (!(ray.flags & RF_SHADOW)) {
		const Triangle& t = triangles[triangleIndex];
		if (normalSmoothing && hasNormals) {
			// if we require normal smoothing, we use the interpolation formula, based on the linear system we just solved.
			info.norm = normals[t.n[0]] 
			+ (normals[t.n[1]] - normals[t.n[0]]) * lambda2
			+ (normals[t.n[2]] - normals[t.n[0]]) * lambda3;
			// .. and we don't forget to normalize it later
			info.norm.normalize();
		} else {
			info.norm = t.gnormal;
		}
		// calculate the UV coordinates, using the same interpolation formula
		Vector uv = uvs[t.t[0]] 
			+ (uvs[t.t[1]] - uvs[t.t[0]]) * lambda2
			+ (uvs[t.t[2]] - uvs[t.t[0]]) * lambda3;
		info.u = uv.x;
		info.v = uv.y;
		if (bumpMap) {
			// use the UV coords to lookup into the bump map, fetch the dx and dy perturbances,
			// and perturb the normal vector
			float u = (float) safeRange(info.u, 0.0, 1.0); // limit UVs to [0..1)
			float v = (float) safeRange(info.v, 0.0, 1.0); // wraparound if necessary
			Color bumpAmount = 
				bumpMap->getFilteredPixel(bumpMap->getWidth() * u, bumpMap->getHeight() * v);
			float dx = bumpAmount.r;
			float dy = bumpAmount.g;
			Vector perturb = (t.dNdx * dx + t.dNdy * dy) * bumpStrength; // calculate bump amount...
			info.norm += perturb;                                        // ...and perturb the normal
			info.norm.normalize();
		}
	}
}

bool Mesh::intersect(const Ray& ray, IntersectionInfo& info) const
{
	// if a ray doesn't intersect the bounding volume, then it doesn't intersect the geometry either
	// (the bounding volume strictly encompasses the mesh)
	if (!bbox.testIntersect(ray)) return false;
	if (kdtreeroot) {
		return intersectKD(*kdtreeroot, bbox, ray, info);
	}
	double closestdist = INF;
	for (int i = 0; i < (int) triangles.size(); i++) {
		const Triangle& t = triangles[i];
		double lambda2, lambda3;
		if (!intersectTriangleFast(ray, vertices[t.v[0]], t, lambda2, lambda3, closestdist))
			continue;
		double gamma = closestdist;
		fillInfo(ray, info, i, lambda2, lambda3, gamma);
	}
	return closestdist < INF;
}

bool Mesh::isInside(const Vector& p) const
{
	/*
	 * How can we check about whether a point is inside a Mesh or not?
	 * For convex meshes, there is a good plane orientation testing method. However,
	 * for general meshes, the only algorithm is to shoot a random ray from that point
	 * towards somewhere and count the number of intersections. If odd, we're inside.
	 * If even, we're outside.
	 * This won't work for any "open" meshes.
	 *
	 * For consistency reasons, I'm actually using a fixed ray (1, 1, 1) instead of random
	 */
	 Ray ray;
	 ray.start = p;
	 ray.dir = Vector(0.693361274351, 0.693361274351, 0.693361274351); // sqrt(1/3)
	 int intersectionCount = 0;
	 IntersectionInfo info;
	 while (intersect(ray, info)) {
		intersectionCount++;
		ray.start = info.ip + ray.dir * 1e-6;
	 }
	 return (intersectionCount % 2) == 1;
}



void Mesh::buildTree(void)
{
	Uint32 start = SDL_GetTicks();
	vector<int> t;
	for (int i = 0; i < (int) triangles.size(); i++)
		t.push_back(i);
	kdtreeroot = new KDTree;
	build(*kdtreeroot, t, bbox, 0);
	Uint32 finish = SDL_GetTicks();
	printf("KDTree built: %d triangles in %.3lf seconds\n", (int) t.size(), (finish - start) / 1000.0);
}

inline bool intersectLeftBox(const BBox& box, const Vector& A, const Vector& B, const Vector& C, Axis axis) {
	return min(min(A[axis], B[axis]), C[axis]) <= box.vmax[axis] + 1e-6 && box.intersectTriangle(A, B, C);
}

inline bool intersectRightBox(const BBox& box, const Vector& A, const Vector& B, const Vector& C, Axis axis) {
	return max(max(A[axis], B[axis]), C[axis]) >= box.vmin[axis] - 1e-6 && box.intersectTriangle(A, B, C);
}

double Mesh::calculateFSplit(double sp, const vector<int>& triangleList, BBox box, Axis splitAxis) const {
	BBox boxLeft, boxRight;
	box.split((Axis) splitAxis, sp, boxLeft, boxRight);
	int numTrianglesLeft = 0;
	int numTrianglesRight = 0;

	for (int i = 0; i < (int) triangleList.size(); i++) {
		const Triangle& t = triangles[triangleList[i]];
		const Vector& A = vertices[t.v[0]];
		const Vector& B = vertices[t.v[1]];
		const Vector& C = vertices[t.v[2]];
		if (intersectLeftBox(boxLeft, A, B, C, (Axis) splitAxis)) {
			numTrianglesLeft++;
		}
		if (intersectRightBox(boxRight, A, B, C, (Axis) splitAxis)) {
			numTrianglesRight++;
		}
	}

	double surfaceAreaLeft = boxLeft.surfaceArea() / box.surfaceArea();
	double surfaceAreaRight = boxRight.surfaceArea() / box.surfaceArea();

	return 0.3 + numTrianglesLeft * surfaceAreaLeft + numTrianglesRight * surfaceAreaRight;
}

void Mesh::build(KDTree& node, const vector<int>& triangleList, BBox box, int depth)
{
	if (triangleList.size() < MAX_TRIANGLES_PER_LEAF ||
	    depth > MAX_TREE_DEPTH) {
		node.axis = AXIS_NONE;
		node.data.triangles = new vector<int>;
		(*node.data.triangles) = triangleList;
	} else {
		int splitAxis = depth % 3; // alternate splitting planes: X, Y, Z, X, Y, Z, ...
		double bbleft = box.vmin[splitAxis]; // find the left and right extents of the bbox along the chosen axis
		double bbright = box.vmax[splitAxis];
		
		double sp;
		
		double fSplit = INF;
		for (int i = 1; i < 20; i++) {
			double cur = bbleft + (bbright - bbleft) * i / 20;
			double curFSplit = calculateFSplit(cur, triangleList, box, (Axis) splitAxis);

			if (curFSplit < fSplit) {
				sp = cur;
				fSplit = curFSplit;
			}
		}
		
		double f = triangleList.size();
		if (f < fSplit) {
			node.axis = AXIS_NONE;
			node.data.triangles = new vector<int>;
			(*node.data.triangles) = triangleList;
			return;
		}

		BBox boxLeft, boxRight;
		box.split((Axis) splitAxis, sp, boxLeft, boxRight);
		vector<int> listLeft, listRight;
		for (int i = 0; i < (int) triangleList.size(); i++) {
			Triangle& t = triangles[triangleList[i]];
			const Vector& A = vertices[t.v[0]];
			const Vector& B = vertices[t.v[1]];
			const Vector& C = vertices[t.v[2]];
			if (intersectLeftBox(boxLeft, A, B, C, (Axis) splitAxis)) {
				listLeft.push_back(triangleList[i]);
			}
			if (intersectRightBox(boxRight, A, B, C, (Axis) splitAxis)) {
				listRight.push_back(triangleList[i]);
			}
		}
		node.axis = (Axis) splitAxis;
		node.sp = sp;
		node.data.children = new KDTree[2];
		build(node.data.children[0], listLeft, boxLeft, depth + 1);
		build(node.data.children[1], listRight, boxRight, depth + 1);
	}
}

bool Mesh::intersectKD(KDTree& node, const BBox& box, const FastRay& ray, IntersectionInfo& info) const
{
	if (node.axis == AXIS_NONE) {
		vector<int>& trios = *node.data.triangles;
		double foundlambda2 = 0, foundlambda3 = 0, dist = INF;
		int foundIdx = -1;
		for (int i = 0; i < (int) trios.size(); i++) {
			const Triangle& t = triangles[trios[i]];
			double l2, l3;
			if (intersectTriangleFast(ray, vertices[t.v[0]], t, l2, l3, dist)) {
				foundIdx = trios[i];
				foundlambda2 = l2;
				foundlambda3 = l3;
			}
		}
		if (foundIdx != -1) {
			Vector ip = ray.start + ray.dir * dist;
			if (box.inside(ip)) {
				fillInfo(ray, info, foundIdx, foundlambda2, foundlambda3, dist);
				return true;
			}
		}
		return false;
	} else {
		BBox boxChild[2];
		box.split(node.axis, node.sp, boxChild[0], boxChild[1]);
		int order[2] = { 0, 1 };
		if (ray.start[node.axis] > node.sp) swap(order[0], order[1]);
		// if the ray intersects the common wall between the two sub-boxes, then it invariably
		// intersects both boxes (we can skip the testIntersect() checks):
		// (see http://raytracing-bg.net/?q=node/68 )
		if (box.intersectWall(node.axis, node.sp, ray)) {
			if (intersectKD(node.data.children[order[0]], boxChild[order[0]], ray, info)) return true;
			return intersectKD(node.data.children[order[1]], boxChild[order[1]], ray, info);
		} else {
			// if the wall isn't hit, then we intersect exclusively one of the sub-boxes;
			// test one, if the test fails, then it's in the other:
			if (boxChild[order[0]].testIntersect(ray))
				return intersectKD(node.data.children[order[0]], boxChild[order[0]], ray, info);
			else
				return intersectKD(node.data.children[order[1]], boxChild[order[1]], ray, info);
		}
	}
	return false;
}

