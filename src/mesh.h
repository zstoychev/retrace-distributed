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
#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include <string>
#include "geometry.h"
#include "bitmap.h"
#include "bbox.h"
using namespace std;

// A node of the K-d tree. It is either a in-node (if axis is AXIS_X, AXIS_Y, AXIS_Z),
// in which case the 'sp' holds the split position, and data.children is an array
// of two children.
// If axis is AXIS_NONE, then it is a leaf node, and data.triangles holds a list of
// indices to triangles.
struct KDTree {
	Axis axis;
	double sp; // splitpos
	union {
		vector<int>* triangles;
		KDTree* children; // children[0], children[1]
	} data;
	
	KDTree()
	{
		axis = AXIS_NONE;
		sp = 0;
		data.triangles = NULL;
	}
	
	~KDTree()
	{
		if (axis == AXIS_NONE) {
			delete data.triangles;
		} else {
			delete[] data.children;
		}
	}
};

struct ObjFileProperties {
	vector<Vector> vertices;
	vector<Vector> normals;
	vector<Vector> uvs;
	vector<Triangle> triangles;
	BBox bbox;
	KDTree* kdtreeroot;
	bool hasNormals;
};

class Mesh: public Geometry {
	vector<Vector> vertices; //!< Vertices of the mesh (1-based, index 0 is dummy)
	vector<Vector> normals;  //!< Normals. 1-based, may contain only one element, the dummy at index 0
	vector<Vector> uvs;      //!< Texture coords (only x,y of the Vectors are used). 1-based, might be empty.
	vector<Triangle> triangles; //!< A list of all triangles that make up the mesh
	
	BBox bbox; //!< A bounding volume of the mesh
	KDTree* kdtreeroot; //!< A K-d tree for faster searching
	
	bool normalSmoothing; //!< Whether normal smoothing is to be enabled or not
	bool hasNormals;      //!< Whether the mesh has normals
	Bitmap* bumpMap;      //!< A pointer to a bump map texture
	float bumpStrength;   //!< The strength of the bumpmap. Applies scaling to the normal perturbation.
	bool useKDTree;       //!< Whether to build and use the K-d tree.
	bool autoSmooth;      //!< Automatically create smoothed normals if needed
	char fn[256];
	char bumpMapFile[256];

	// fils the IntersectionInfo for the specified triangle. lambda2, lambda3 and gamma are the paramters, found by
	// intersectTriangle() (gamma==closestdist)
	void fillInfo(const Ray& ray, IntersectionInfo& info, int triangleIndex, double lambda2, double lambda3, double gamma) const;

	void buildTree(void); // builds the K-d tree (called from loadFromOBJ())
	// a recursive internal function to build the K-d tree.
	void build(KDTree& node, const vector<int>& triangleList, BBox box, int depth);
	double calculateFSplit(double sp, const vector<int>& triangleList, BBox box, Axis splitAxis) const;
	// a recursive internal function to search the K-d tree for intersections
	// intersect a `node', known to be bounded by `box' with the givne `ray', and if an
	// intersection is found, update the `info' and return true
	bool intersectKD(KDTree& node, const BBox& box, const FastRay& ray, IntersectionInfo& info) const;

public:
	Mesh() {
		normalSmoothing = true;
		bumpMap = NULL;
		bumpStrength = 1;
		kdtreeroot = NULL;
		useKDTree = true;
		autoSmooth = false;
	}
	~Mesh() 
	{
		//if (kdtreeroot) delete kdtreeroot;
	}
	bool loadFromOBJ(const char* filename);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getBoolProp("useKDTree", &useKDTree);
		pb.getBoolProp("autoSmooth", &autoSmooth);
		char fullFN[256];
		if (pb.getFilenameProp("file", fn, fullFN))
			if (!loadFromOBJ(fullFN) || vertices.size() <= 1) {
				char msg[256];
				sprintf(msg, "Couldn't load OBJ file `%s' (bad format?)\n", fullFN);
				pb.signalError(msg);
			}
		pb.getBoolProp("normalSmoothing", &normalSmoothing);
		pb.getBitmapFileProp("bumpMap", &bumpMap, bumpMapFile);
		pb.getFloatProp("bumpStrength", &bumpStrength);
		if (bumpMap) {
			Lock lock(cachedItemsMutex);
			bumpMap->differentiate();
		}
	}

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Mesh");
		se.addProperty("useKDTree", useKDTree);
		se.addProperty("autoSmooth", autoSmooth);
		se.addFilenameProperty("file", fn);
		se.addProperty("normalSmoothing", normalSmoothing);
		if (bumpMap != 0) {
			se.addFilenameProperty("bumpMap", bumpMapFile);
			se.addProperty("bumpStrength", bumpStrength);
		}
	}
	
	bool intersect(const Ray& ray, IntersectionInfo& info) const;
	bool isInside(const Vector& p) const;
	
	const char* getName() const { return "Mesh"; }
	
};

/**
 * Intersects a ray with a triangle. Fills lambda2, lambda3 and closestdist if an intersection is found.
 * @param ray - the ray to intersect the triangle with
 * @param A - vertex A of the triangle
 * @param T - the triangle itself (needs AB, AC)
 * @param lambda2 [out] - the barycentric coordinate along AB is written here, if intersection is found
 * @param lambda3 [out] - the barycentric coordinate along AC is written here, if intersection is found
 * @param closestdist [in/out] - the distance to the intersection point is written here, if intersection
                                 is found, and the distance found is less than the current value of
                                 `closestdist'
 * @returns: true if an intersection is found, false otherwise
 */
bool intersectTriangleFast(const Ray& ray, const Vector& A, const Triangle& T, double& lambda2, double& lambda3, double& closestdist);


#endif // __MESH_H__
