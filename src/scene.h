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
#ifndef __SCENE_H__
#define __SCENE_H__

#include <vector>
#include <limits.h>
#include <cstring>
#include <map>
#include <string>
#include "color.h"
#include "vector.h"
#include "matrix.h"
#include "cxxptl_sdl.h"

enum ElementType {
	ELEM_GEOMETRY,
	ELEM_SHADER,
	ELEM_NODE,
	ELEM_TEXTURE,
	ELEM_ENVIRONMENT,
	ELEM_LIGHT,
	ELEM_CAMERA,
	ELEM_SETTINGS,
};

class SceneParser;
class Geometry;
class Intersectable;
class Shader;
class Node;
class Texture;
class Environment;
class Camera;
class Light;
class Bitmap;
struct Transform;

class ParsedBlock;
class SceneEntry;

struct Scene;

/// An abstract base class for each element of the scene (Camera, Geometries,...)
/// i.e, anything, that could be described using our scene definition language
class SceneElement {
public:
	Scene* scene;
	char name[64]; //!< A name of this element (a string like "sphere01", "myCamera", etc)
	SceneElement(); //!< A constructor. It sets the name to the empty string.
	virtual ~SceneElement() {} //!< a virtual destructor
	
	virtual ElementType getElementType() const = 0; //!< Gets the element type
	
	/**
	 * @brief set all the properties of a scene element from a parsed block
	 *
	 * This is a callback, called by the sceneparser, when it has finished
	 * parsing a block of properties, related to some scene element.
	 *
	 * Consider the following (part of) scene:
	 *
	 * Sphere mySphere01 {
	 *    center (12.5, 0.1, 0.3)
	 *    radius 5.0
	 * }
	 *
	 * A class Sphere should inherit from SceneElement and implement fillProperties() in the following manner
	 *
	 * class Sphere: public SceneElement {
	 *	Vector pos;
	 *	double radius;
	 * public:
	 * 	void fillProperties(ParsedBlock& pb) {
	 *		pb.getVectorProp("center", &pos);
	 *		pb.getDoubleProp("radius", &radius);
	 *	}
	 * };
	 *
	 * In fact, one would usually want to have a constructor that initializes the two
	 * properties to some default values.
	 * Also, for the sake of consistency, one should really rename that `pos' to `center' in
	 * the class's definition.
	 *
	 * (the implementation of SceneElement::fillProperties() does nothing)
	 */
	virtual void fillProperties(ParsedBlock& pb);

	virtual void fillSceneEntry(SceneEntry& se) const {}
	
	/**
	 * @brief a callback that gets called before the rendering commences
	 *
	 * If you need to setup some internal data structures before rendering has begun,
	 * you should place it here. This callback is executed after scene parsing, and before
	 * the rendering commences. You might actually use other SceneElement's, but the beginFrame()
	 * function is called for the different classes of SceneElement's in this specific order:
	 *
	 * 1) Lights
	 * 2) Geometries
	 * 3) Textures
	 * 4) Shaders
	 * 5) Nodes
	 * 6) Camera
	 * 7) GlobalSettings
	 *
	 * The order of calling beginFrame within the same group is undefined.
	 *
	 * All these callbacks are called by the Scene::beginRender() function.
	 */
	virtual void beginRender();
	
	/**
	 * @brief same as beginRender(), but gets called before each frame
	 *
	 * the difference between beginRender() and beginFrame() is that beginRender() is only
	 * called once, after parsing is done, whereas beginFrame is called before every frame
	 * (e.g., when rendering an animation).
	 */
	virtual void beginFrame();
	
	friend class SceneParser;
};

class ParsedBlock {
public:
	virtual ~ParsedBlock() {}
	// All these methods are intended to be called by SceneElement implementations
	// each method accepts two parameters: a name, and a value pointer.
	// Each method does one of these things:
	//  - the property with the given name is found and parsed successfully, in which
	//    case the value is filled in, and the method returns true.
	//  - The property is found and it wasn't parsed successfully, in which case
	//    a syntax error exception is raised.
	//  - The property is missing in the scene file, the value is untouched, and the method
	//    returns false (if this is an error, you can signal it with pb.requiredProp(name))
	//
	// Some properties also have min/max ranges. If they are specified and the parsed value
	// does not pass the range check, then a SyntaxError is raised.
	virtual bool getIntProp(const char* name, int* value, int minValue = INT_MIN, int maxValue = INT_MAX) = 0;
	virtual bool getBoolProp(const char* name, bool* value) = 0;
	virtual bool getFloatProp(const char* name, float* value, float minValue = -LARGE_FLOAT, float maxValue = LARGE_FLOAT) = 0;
	virtual bool getDoubleProp(const char* name, double* value, double minValue = -LARGE_DOUBLE, double maxValue = LARGE_DOUBLE) = 0;
	virtual bool getColorProp(const char* name, Color* value, float minCompValue = -LARGE_FLOAT, float maxCompValue = LARGE_FLOAT) = 0;
	virtual bool getVectorProp(const char* name, Vector* value) = 0;
	virtual bool getGeometryProp(const char* name, Geometry** value) = 0;
	virtual bool getIntersectableProp(const char* name, Intersectable** value) = 0;
	virtual bool getShaderProp(const char* name, Shader** value) = 0;
	virtual bool getTextureProp(const char* name, Texture** value) = 0;
	virtual bool getNodeProp(const char* name, Node** value) = 0;
	virtual bool getStringProp(const char* name, char* value) = 0; // the buffer should be 256 chars long
	
	// useful for scene assets like textures, mesh files, etc.
	// the value will hold the full filename to the file.
	// If the file/dir is not found, a FileNotFound exception is raised.
	virtual bool getFilenameProp(const char* name, char* value, char* fullPath) = 0;
	
	// Does the same logic as getFilenameProp(), but also loads the bitmap
	// file from the specified file name. The given bitmap is first deleted if not NULL.
	virtual bool getBitmapFileProp(const char* name, Bitmap** value, char* filename) = 0;
	
	// Gets a transform from the parsed block. Namely, it searches for all properties named
	// "scale", "rotate" and "translate" and applies them to T. IT is the inverse of T and
	// is computed after the last modification to T.
	virtual void getTransformProp(Transform& T, Transform& IT) = 0;
	
	virtual void requiredProp(const char* name) = 0; // signal an error (missing property of the given name)
	
	virtual void signalError(const char* msg) = 0; // signal an error with a specified message
	virtual void signalWarning(const char* msg) = 0; // signal a warning with a specified message
	
	// some functions for direct parsed block access:
	virtual int getBlockLines() = 0;
	virtual void getBlockLine(int idx, int& srcLine, char head[], char tail[]) = 0;
	virtual SceneParser& getParser() = 0;
};

class SceneParser {
public:
	virtual ~SceneParser() {}
	// All these methods are intended to be called by SceneElement implementations:
	virtual Shader* findShaderByName(const char* name) = 0;
	virtual Texture* findTextureByName(const char* name) = 0;
	virtual Geometry* findGeometryByName(const char* name) = 0;
	virtual Node* findNodeByName(const char* name) = 0;
	
	
	/**
	 * resolveFullPath() tries to find a file (or folder), by appending the given path to the directory, where
	 * the scene file resides. The idea is that all external files (textures, meshes, etc.) are
	 * stored in the same dir where the scene file (*.retrace) resides, and the paths to that external
	 * files do not mention any directories, just the file names.
	 *
	 * @param path (input-output) - Supply the given file name (as given in the scene file) here.
	 *                              If the function succeeds, this will return the full path to a file
	 *                              with that name, if it's found.
	 * @returns true on success, false on failure (file not found).
	 */
	virtual bool resolveFullPath(char* path) = 0;
};

struct SyntaxError {
	char msg[128];
	int line;
	SyntaxError();
	SyntaxError(int line, const char* msg1, const char* msg2);
	SyntaxError(int line, const char* format, ...);
};

struct FileNotFoundError {
	char filename[245];
	int line;
	FileNotFoundError();
	FileNotFoundError(int line, const char* filename);
};

/// Utility function: gets three doubles from a string in just the same way as the sceneparser will do it
/// (stripping any commas, parentheses, etc)
/// If there is error in parsing, a SyntaxError exception is raised.
void get3Doubles(int srcLine, char* expression, double& d1, double& d2, double& d3);

/// Splits a string (given in an expression) to a head token and a remaining. The head token is written into
/// `frontToken', and the remaining is copied back to expression
bool getFrontToken(char* expression, char* frontToken);

/// Splits a string (given in an expression) to a back token and a remaining. The back token is written into
/// `backToken', and the remaining is copied back to expression
bool getLastToken(char* expression, char* backToken);

void stripPunctuation(char* expression); //!< strips any whitespace or punctuation in front or in the back of a string (in-place)

struct SceneEntryProperty {
	char name[128];
	char value[256];

	SceneEntryProperty(const char* name, const char* value) {
		strncpy(this->name, name, sizeof(this->name) - 1);
		this->name[sizeof(this->name) - 1] = '\0';
		strncpy(this->value, value, sizeof(this->value) - 1);
		this->value[sizeof(this->value) - 1] = '\0';
	}
};

class SceneEntry {
	const char* type;
	const char* name;
	std::vector<SceneEntryProperty> properties;
public:
	SceneEntry() : type("Unknown"), name("") {}

	void setType(const char* type) {
		this->type = type;
	}

	const char* getType() const {
		return type;
	}

	void setName(const char* name) {
		this->name = name;
	}

	const char* getName() const {
		return name;
	}

	const std::vector<SceneEntryProperty>& getProperties() const {
		return properties;
	}

	void addProperty(const char* name, const char* value) {
		properties.push_back(SceneEntryProperty(name, value));
	}

	void addProperty(const char* name, double value) {
		char valueAsString[256];
		sprintf(valueAsString, "%f", value);
		properties.push_back(SceneEntryProperty(name, valueAsString));
	}

	void addProperty(const char* name, int value) {
		char valueAsString[256];
		sprintf(valueAsString, "%d", value);
		properties.push_back(SceneEntryProperty(name, valueAsString));
	}

	void addProperty(const char* name, const Vector& value) {
		char valueAsString[256];
		sprintf(valueAsString, "(%f, %f, %f)", value.x, value.y, value.z);
		properties.push_back(SceneEntryProperty(name, valueAsString));
	}

	void addProperty(const char* name, const Color& value) {
		char valueAsString[256];
		sprintf(valueAsString, "(%f, %f, %f)", value.r, value.g, value.b);
		properties.push_back(SceneEntryProperty(name, valueAsString));
	}

	void addProperty(const char* name, const SceneElement* value) {
		properties.push_back(SceneEntryProperty(name, value->name));
	}

	void addFilenameProperty(const char* name, const char* value) {
		char valueAsString[256];
		sprintf(valueAsString, "\"%s\"", value);
		properties.push_back(SceneEntryProperty(name, valueAsString));
	}

	void addProperty(const Transform& value) {
		addProperty("transformationMatrix0", Vector(value.m.m[0][0], value.m.m[0][1], value.m.m[0][2]));
		addProperty("transformationMatrix1", Vector(value.m.m[1][0], value.m.m[1][1], value.m.m[1][2]));
		addProperty("transformationMatrix2", Vector(value.m.m[2][0], value.m.m[2][1], value.m.m[2][2]));
		addProperty("translate", value.offset);
	}
};

/// This structure holds all global settings of the scene - frame size, antialiasing toggles, thresholds, etc...
struct GlobalSettings: public SceneElement {
	int frameWidth, frameHeight; //!< Frame sizes

	Color ambientLight;          //!< ambient color
	bool wantAA, wantPrepass;    //!< Antialiasing flag and prepass (a quick low-resolution rendering) flag
	double aaThresh;             //!< The antialiasing color difference threshold (see renderScene)
	
	int maxTraceDepth;           //!< Maximum recursion depth
	
	bool dbg;                    //!< A debugging flag (if on, various raytracing-related procedures will dump debug info to stdout).
	bool gi;                     //!< Global illumination (default: off)
	int pathsPerPixel;           //!< Paths per pixel when GI is on
	float outputGamma;           //!< Gamma correction: assumed output gamma. Defaults to 1 (no color change)
	bool interactive;            //!< Interactive mode on
	int renderThreads;           //!< How many threads to use in multi-threading.
	
	GlobalSettings();
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("GlobalSettings");
		se.addProperty("frameWidth", frameWidth);
		se.addProperty("frameHeight", frameHeight);
		se.addProperty("ambientLight", ambientLight);
		se.addProperty("maxTraceDepth", maxTraceDepth);
		se.addProperty("dbg", dbg);
		se.addProperty("wantPrepass", wantPrepass);
		se.addProperty("wantAA", wantAA);
		se.addProperty("aaThresh", aaThresh);
		se.addProperty("gi", gi);
		se.addProperty("pathsPerPixel", pathsPerPixel);
		se.addProperty("outputGamma", outputGamma);
		se.addProperty("interactive", interactive);
		se.addProperty("renderThreads", renderThreads);
	}
	ElementType getElementType() const { return ELEM_SETTINGS; }
};

struct Scene {
	std::vector<Geometry*> geometries;
	std::vector<Shader*> shaders;
	std::vector<Node*> nodes;
	std::vector<Node*> subNodes; // also Nodes, but without a shader attached; don't represent an scene object directly
	std::vector<Texture*> textures;
	std::vector<Light*> lights;
	Environment* environment;
	Camera* camera;
	GlobalSettings settings;
	
	Scene();
	~Scene();
	
	bool parseScene(const char* sceneFile, std::string* fileContent = 0); //!< Parses a scene file and loads the scene from it. Returns true on success.
	std::string getAsText() const;
	void saveToFile(const char* fileName) const;
	void beginRender(); //!< Notifies the scene so that a render is about to begin. It calls the beginRender() method of all scene elements
	void beginFrame(); //!< Notifies the scene so that a new frame is about to begin. It calls the beginFrame() method of all scene elements
};

struct SceneProperties {
	int id;
	Scene scene;
	Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE];
	std::string content;
	int finishedFragments;
	int totalFragments;
	bool inAA;

	SceneProperties() {
		memset(vfb, 0, sizeof(vfb));
	}

	bool isReady() {
		return finishedFragments >= totalFragments;
	}
};

extern Scene* volatile currentScene;
extern SceneProperties* volatile currentSP;

struct Rect {
	int x0, y0, x1, y1, w, h;
	Rect() {}
	Rect(int _x0, int _y0, int _x1, int _y1)
	{
		x0 = _x0, y0 = _y0, x1 = _x1, y1 = _y1;
		h = y1 - y0;
		w = x1 - x0;
	}
	void clip(int maxX, int maxY); // clips the rectangle against image size
};

class AANeed {
public:
	unsigned char needsAA[MAX_BUCKET_SIZE][(MAX_BUCKET_SIZE + 7) / 8];

	AANeed() {
		memset(needsAA, 0, sizeof(needsAA));
	}

	bool get(int i, int j) {
		return (int(needsAA[i][j / 8]) & (1 << (j % 8))) != 0;
	}

	void set(int i, int j, bool value) {
		int flag = (value ? 1 : 0) << (j % 8);

		needsAA[i][j / 8] = unsigned char((int(needsAA[i][j / 8]) & (~(1 << (j % 8)))) | flag);
	}
};

struct FragmentData {
	SceneProperties* sceneProperties;
	Rect fragment;
	Color data[MAX_BUCKET_SIZE][MAX_BUCKET_SIZE];
	AANeed needsAA;

	FragmentData() {}
	FragmentData(SceneProperties* sceneProperties, const Rect& fragment)
		: sceneProperties(sceneProperties), fragment(fragment) {}
};

// generate a list of buckets (image sub-rectangles) to be rendered, in a zigzag pattern
std::vector<Rect> getBucketsList(int W, int H);

extern bool disableAcceleratedStructures;
extern std::map<std::string, Bitmap*> filenameToBitmap;
//extern std::map<std::string, Heightfield::HighStruct*> filenameToHighStruct;
struct ObjFileProperties;
extern std::map<std::string, ObjFileProperties*> filenameToObjFileProperties;
extern Mutex cachedItemsMutex;
void freeCachedItems();

#endif // __SCENE_H__
