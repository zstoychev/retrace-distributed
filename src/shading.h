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
#ifndef __SHADING_H__
#define __SHADING_H__


#include "color.h"
#include "geometry.h"
#include "scene.h"

/**
 * @class BRDF - An abstract class, representing a Bidirectional Reflectance Distribution
 *               Function. The real BRDF (as per definition) is implemented in the eval() method.
 */
class BRDF {
public:
	virtual ~BRDF() {}
	/**
	 * @brief Calculate the BRDF, as by definition.
	 * @param x - the point in the scene under consideration
	 * @param w_in - incoming ray (w_in.dir points TOWARDS x)
	 * @param w_out - outgoing ray direction (w_out points OUTWARDS x)
	 * @return The BRDF value for the given parameter.
	 */
	virtual Color eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const;
	/**
	 * @brief Sample the BRDF.
	 * Produces an outgoing ray and calculate the BRDF for the given incoming and the generated
	 * outgoing rays.
	 * @param x - the point in the scene under consideration
	 * @param w_in - incoming ray (w_in.dir points TOWARDS x)
	 * @param w_out [out] - the generated outgoing ray will be written here. The BRDF may or may not
	 *                      use importance sampling. If it does, the generated rays follow the
	 *                      probabilistic distribution of the probability density function (PDF) for
	 *                      a reflection given the intersection point and the incoming direction.
	 *                      That is, more probable reflection rays would be more frequently generated.
	 * @param brdf [out] - the BRDF, evaluated at (x, w_in, w_out). This is mostly the same as the
	 *                     the result of eval(x, w_in, w_out.dir) after spawnRay() generated the w_out;
	 *                     Thus it is somewhat redundant - but it's here for speed and convenience.
	 * @param pdf [out] - the probability to generate this specific w_out ray; I.e., an evaluation
	 *                    of the PDF of reflectance at direction w_out.dir, given x and w_in.
	 */
	virtual void spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
		Color& brdf, float& pdf) const;
};

/// An abstract class, representing a shader in our scene.
class Shader: public SceneElement, public BRDF {
public:
	virtual ~Shader() {}
	
	ElementType getElementType() const { return ELEM_SHADER; }
	
	virtual Color shade(const Ray& ray, const IntersectionInfo& info) = 0;
};

/// An abstract class, representing a (2D) texture
class Texture: public SceneElement {
public:
	virtual ~Texture() {}
	ElementType getElementType() const { return ELEM_TEXTURE; }
	virtual Color getTexColor(const Ray& ray, const IntersectionInfo& info) = 0;
};

/// A checker texture
class Checker: public Texture {
	Color col1, col2; /// the colors of the alternating squares
	double squareSize; /// the size of a square side, in world units
public:
	Checker() { col1 = Color(1.0f, 1.0f, 1.0f); col2 = Color(0.0f, 0.0f, 0.0f); squareSize = 1.0; } // default BW checker
	Color getTexColor(const Ray& ray, const IntersectionInfo& info);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getColorProp("col1", &col1);
		pb.getColorProp("col2", &col2);
		pb.getDoubleProp("squareSize", &squareSize);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Checker");
		se.addProperty("col1", col1);
		se.addProperty("col2", col2);
		se.addProperty("squareSize", squareSize);
	}
};

class Bitmap;
class BitmapTexture: public Texture {
	Bitmap* map;
	double scaling;
	char fileName[256];
public:
	BitmapTexture() { map = NULL; scaling = 1; }
	~BitmapTexture();
	Color getTexColor(const Ray& ray, const IntersectionInfo& info);
	void fillProperties(ParsedBlock& pb)
	{
		pb.getBitmapFileProp("file", &map, fileName);
		pb.getDoubleProp("scaling", &scaling);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("BitmapTexture");
		se.addFilenameProperty("file", fileName);
		se.addProperty("scaling", scaling);
	}
};

/// A Lambert (flat) shader
class Lambert: public Shader {
	Color color; //!< This is the static color of the Lambert shader (to be used if a texture isn't present)
	Texture *texture; //!< a diffuse texture, if not NULL.
public:
	Lambert() { color = Color(1, 1, 1); texture = NULL; }
	Lambert(const Color& color) : color(color), texture(NULL) {}
	Color shade(const Ray& ray, const IntersectionInfo& info);

	void setColor(const Color& color) {
		this->color = color;
	}
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getColorProp("color", &color);
		pb.getTextureProp("texture", &texture);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Lambert");
		se.addProperty("color", color);
		if (texture != 0) {
			se.addProperty("texture", texture);
		}
	}
	Color eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const;
	void spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
		Color& brdf, float& pdf) const;

};

/// A Lambert (flat) shader
class Phong: public Shader {
	Color color; //!< This is the static color of the Lambert shader (to be used if a texture isn't present)
	Texture *texture; //!< a diffuse texture, if not NULL.
	double exponent;
	float specularMultiplier; //!< specular amount, defaults to 1
public:
	Phong() { color = Color(1, 1, 1); texture = NULL; exponent = 100.0; specularMultiplier = 1.0f; }
	Color shade(const Ray& ray, const IntersectionInfo& info);
	void fillProperties(ParsedBlock& pb)
	{
		pb.getColorProp("color", &color);
		pb.getTextureProp("texture", &texture);
		pb.getDoubleProp("exponent", &exponent, 0.5, 300000.0);
		pb.getFloatProp("specularMultiplier", &specularMultiplier);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Phong");
		se.addProperty("color", color);
		if (texture != 0) {
			se.addProperty("texture", texture);
		}
		se.addProperty("exponent", exponent);
		se.addProperty("specularMultiplier", specularMultiplier);
	}
	Color eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const;
	void spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
		Color& brdf, float& pdf) const;

};


class Reflection: public Shader {
	double glossiness, scaling;
	float lightMultiplier;
	int samples;
public:
	Reflection()
	{
		this->glossiness = 1;
		lightMultiplier = 0.8f;
		samples = 32;
	}
	Color shade(const Ray& ray, const IntersectionInfo& info);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("glossiness", &glossiness, 0, 1);
		scaling = tan((1 - glossiness) * PI /2);
		pb.getIntProp("samples", &samples, 1);
		pb.getFloatProp("lightMultiplier", &lightMultiplier, 0, 1);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Reflection");
		se.addProperty("glossiness", glossiness);
		se.addProperty("samples", samples);
		se.addProperty("lightMultiplier", lightMultiplier);
	}
	Color eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const;
	void spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
		Color& brdf, float& pdf) const;

};

class Refraction: public Shader {
	float ior;
	float lightMultiplier;
public:
	Refraction() { ior = 1.33f; lightMultiplier = 0.8f; }
	Color shade(const Ray& ray, const IntersectionInfo& info);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getFloatProp("ior", &ior, 1e-6f);
		pb.getFloatProp("lightMultiplier", &lightMultiplier, 0, 1);
	}
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Refraction");
		se.addProperty("ior", ior);
		se.addProperty("lightMultiplier", lightMultiplier);
	}
	Color eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const;
	void spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
		Color& brdf, float& pdf) const;

};

class Layered: public Shader {
	struct Layer {
		Shader * shader;
		Color blend;
		Texture* blendTex;
	};
	static const int MAX_LAYERS = 32;
	Layer layers[MAX_LAYERS];
	int nLayers;
public:
	Layered() { nLayers = 0; }
	void addLayer(Shader* shader, Color blend, Texture* blendTex = NULL);
	Color shade(const Ray& ray, const IntersectionInfo& info);
	void fillProperties(ParsedBlock& pb);
	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Layered");
		for (int i = 0; i < nLayers; i++) {
			char value[256];
			sprintf(value, "%s (%f, %f, %f) %s", layers[i].shader->name,
				layers[i].blend.r, layers[i].blend.g, layers[i].blend.b,
				layers[i].blendTex == 0 ? "" : layers[i].blendTex->name);
			se.addProperty("layer", value);
		}
	}
};

class Fresnel: public Texture {
	float ior;
public:
	Fresnel() { ior = 1.33f; }
	Color getTexColor(const Ray& ray, const IntersectionInfo& info);
	
	void fillProperties(ParsedBlock& pb)
	{
		pb.getFloatProp("ior", &ior, 1e-6f);
	}

	virtual void fillSceneEntry(SceneEntry& se) const {
		se.setType("Fresnel");
		se.addProperty("ior", ior);
	}
};

#endif // __SHADING_H__
