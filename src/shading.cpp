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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "shading.h"
#include "bitmap.h"
#include "lights.h"
#include "random_generator.h"

// "dummy" eval functions. Used where the BRDF for a specific Shader is not yet written.
// pass some special values to indicate an error.
Color BRDF::eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const
{
	return Color(1, 0, 0);
}
void BRDF::spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
	Color& brdf, float& pdf) const
{
	pdf = -1;
}


Color Checker::getTexColor(const Ray& ray, const IntersectionInfo& info)
{
	/*
	 * The checker texture works like that. Partition the whole 2D space
	 * in squares of squareSize. Use division and floor()ing to get the
	 * integral coordinates of the square, which our point happens to be. Then,
	 * use the parity of the sum of those coordinates to decide which color to return.
	 */
	int squareX = (int) floor(info.u / squareSize);
	int squareY = (int) floor(info.v / squareSize);
	if ((squareX + squareY) % 2 == 0) return col1;
	else return col2;
}

bool lightIsVisible(Vector p, Vector l, Scene& scene);

Color Lambert::shade(const Ray& ray, const IntersectionInfo& info)
{
	// fetch the material color. This is ether the solid color, or a color
	// from the texture, if it's set up.
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(ray, info) * color;
	else materialColor = color;
	
	Color lightMultiplier = scene->settings.ambientLight;
	
	for (int i = 0; i < (int) scene->lights.size(); i++) {
		Vector lightPos;
		Color lightCol;
		int n = scene->lights[i]->getNumSamples();
		Color accum(0.0f, 0.0f, 0.0f);
		for (int j = 0; j < n; j++) {
			scene->lights[i]->getNthSample(j, info.ip, lightPos, lightCol);
			Vector lightDir = lightPos - info.ip;
			double lightDist = lightDir.length();
			lightDir.normalize();
			// get the Lambertian cosine of the angle between the geometry's normal and
			// the direction to the light. This will scale the lighting:
			double normDotL = lightDir * info.norm;
			if (normDotL > 0 &&    // we require the light to be "above" the surface.
			lightIsVisible(info.ip, lightPos, *scene)) {// also check if our point (info.ip) is visible from the light...
				// ... in which case add the diffuse component, with the quadratic light attenuation law applied.
				accum += lightCol * float(normDotL / sqr(lightDist));
			}
		}
		lightMultiplier += accum / float(n);
	}
	return materialColor * lightMultiplier;
}

/* generate a random, uniformly-distributed point in the unit hemisphere,
 * such that dot(<generated point>, norm) >= 0. It can be also thought as, if
 * we have the surface normal at some scene point, generate a unit direction
 * out of that point, with uniform distribution. Used in Lambert and Phong BRDFs
 */
static inline Vector hemisphereSample(const Vector& norm)
{
	Random& rnd = getRandomGen();
	
	// first, generate a random sphere point (http://mathworld.wolfram.com/SpherePointPicking.html
	float u = rnd.randfloat();
	float v = rnd.randfloat();
	float theta = 2 * (float) PI * u;
	
	float cosPhi = 2 * v - 1;
	float sinPhi = sqrt(1 - sqr(cosPhi));
	
	Vector p;
	p.x = cos(theta) * sinPhi;
	p.y = cosPhi;
	p.z = sin(theta) * sinPhi;
	
	// Test if the generated point is in the right hemisphere; if not, just flip the signs.
	// This is guaranteed to make the dot product positive.
	if (dot(p, norm) < 0) p = -p;
	return p;
}


Color Lambert::eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const
{
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(w_in, x) * color;
	else materialColor = color;
	
	return materialColor / (float) PI;
}

void Lambert::spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
	Color& brdf, float& pdf) const
{
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(w_in, x) * color;
	else materialColor = color;
	
	w_out = w_in;
	w_out.depth++;
	w_out.start = x.ip + x.norm * 1e-6;
	w_out.dir = hemisphereSample(x.norm);
	w_out.flags |= RF_DIFFUSE;
	
	pdf = 1.0f / float(PI);
	brdf = materialColor * pdf;
}


Color Phong::shade(const Ray& ray, const IntersectionInfo& info)
{
	// fetch the material color. This is ether the solid color, or a color
	// from the texture, if it's set up.
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(ray, info) * color;
	else materialColor = color;
	
	Vector nrm = faceforward(info.norm, ray.dir);
	
	Color diffuse(0.0f, 0.0f, 0.0f);
	Color specular(0.0f, 0.0f, 0.0f);
	
	for (int i = 0; i < (int) scene->lights.size(); i++) {
		Color lightColor;
		Vector lightPos;
		Color accumDiffuse(0.0f, 0.0f, 0.0f);
		Color accumSpecular(0.0f, 0.0f, 0.0f);
		int n = scene->lights[i]->getNumSamples();
		for (int j = 0; j < n; j++) {
			scene->lights[i]->getNthSample(j, info.ip, lightPos, lightColor);
			if (!lightIsVisible(info.ip, lightPos, *scene)) continue; // skip this light - it isn't visible at all
		
			Vector lightDir = lightPos - info.ip;
			double lightDist = lightDir.length();
			lightDir.normalize();
			// get the Lambertian cosine of the angle between the geometry's normal and
			// the direction to the light. This will scale the lighting:
			double normDotL = lightDir * nrm;
			if (normDotL > 0) { // the light must be "above" the surface.
				accumDiffuse += lightColor * float(normDotL / sqr(lightDist));
				// calculate specular:
				Vector fromLight = info.ip - lightPos;
				fromLight.normalize();
			
				Vector r = reflect(fromLight, nrm);
				Vector toCamera = ray.start - info.ip;
				toCamera.normalize();
				double cosGamma = dot(toCamera, r);
				accumSpecular += lightColor * float(pow(cosGamma, exponent) / sqr(lightDist));
			}
		}
		float rN = 1.0f / n;
		diffuse += accumDiffuse * rN;
		specular += accumSpecular * (specularMultiplier * rN);
	}
	// phong shader: combine ambient, diffuse and specular:
	return (scene->settings.ambientLight + diffuse) * materialColor + specular;
}


Color Phong::eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const
{
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(w_in, x) * color;
	else materialColor = color;
	
	Vector nrm = faceforward(x.norm, w_in.dir);
	Vector r = reflect(-w_out, nrm);
	Vector toCamera = -w_in.dir;
	float cosGamma = (float) dot(toCamera, r);

	
	return (materialColor + Color(1, 1, 1) * specularMultiplier * pow(cosGamma, (float) exponent)) / (float) PI;
}

void Phong::spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
	Color& brdf, float& pdf) const
{
	Color materialColor;
	if (texture != NULL) materialColor = texture->getTexColor(w_in, x) * color;
	else materialColor = color;
	
	w_out = w_in;
	w_out.depth++;
	w_out.start = x.ip + x.norm * 1e-6;
	w_out.dir = hemisphereSample(x.norm);
	w_out.flags |= RF_DIFFUSE;
	
	Vector nrm = faceforward(x.norm, w_in.dir);
	Vector r = reflect(-w_out.dir, nrm);
	Vector toCamera = -w_in.dir;
	float cosGamma = (float) dot(toCamera, r);
	
	pdf = 1 / (float)PI;
	brdf = (materialColor + Color(1, 1, 1) * specularMultiplier * pow(cosGamma, (float) exponent)) * pdf;
}


BitmapTexture::~BitmapTexture(){ //if (map) delete map; }
}

Color BitmapTexture::getTexColor(const Ray& ray, const IntersectionInfo& info)
{
	if (!map) return Color(0, 0, 0);
	double u = info.u / scaling;
	double v = info.v / scaling;
	double fracu = u - floor(u);
	double fracv = v - floor(v);
	float tx = (float) fracu * map->getWidth();
	float ty = (float) fracv * map->getHeight();
	return map->getFilteredPixel(tx, ty);
}

Color Reflection::shade(const Ray& ray, const IntersectionInfo& info)
{
	extern Color raytrace(const Ray& ray, Scene& scene);
	
	Vector n = faceforward(info.norm, ray.dir);
	if (glossiness == 1) {
	 	// simple reflection, launch a single ray:
		Vector r = reflect(ray.dir, n);
		Ray newRay = ray;
		newRay.dir = r;
		newRay.start = info.ip + n * 1e-6;
		newRay.depth = ray.depth + 1;
		return raytrace(newRay, *scene) * lightMultiplier;
	} else {
		Random& rnd = getRandomGen();
		Vector dx, dy;
		// generate an orthonormed system; the new vectors dx and dy will be orthogonal
		// to each other, and to n, in the same time.
		orthonormedSystem(n, dx, dy);
		Color sum(0, 0, 0);
		// If we're considering a primary ray, use the full-quality sampling to achieve
		// best result. If not (e.g., we came here via a reflection or refraction), then
		// use five samples only, to avoid a potential combinatorial explosion caused by
		// inter-reflecting surfaces with glossiness < 1:
		int actual_samples = (ray.flags & RF_GLOSSY) ? 5 : this->samples;
		int samples_done = 0; //< count of samples actually traced
		while (samples_done < actual_samples) {
			double x, y;
			// get a random point on the unit disc and scale it:
			rnd.unitDiscSample(x, y);
			x *= scaling;
			y *= scaling;
			// modify the normal according to the random offset:
			Vector n_new = n + x * dx + y * dy;
			n_new.normalize();
			// reflect the incoming ray around the new normal:
			Vector r = reflect(ray.dir, n_new);
			
			// the resulting reflection might "dig in" below the actual geometry surface;
			// if this is the case, just ignore that ray, and only trace "real" reflection rays:
			if (dot(r, n) > 0) {
				Ray newRay = ray;
				newRay.dir = r;
				newRay.start = info.ip + n * 1e-6;
				newRay.depth = ray.depth + 1;
				newRay.flags |= RF_GLOSSY;
				sum += raytrace(newRay, *scene);
				// we added a new sample, increment the count
				samples_done++;
			}
		}
		return sum * (lightMultiplier / actual_samples);
	}
}

Color Reflection::eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const
{
	return Color(0, 0, 0);
}
void Reflection::spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
	Color& brdf, float& pdf) const
{
	Vector n = faceforward(x.norm, w_in.dir);
	// simple reflection, launch a single ray:
	Vector r = reflect(w_in.dir, n);
	w_out = w_in;
	w_out.depth++;
	w_out.dir = r;
	w_out.start = x.ip + n * 1e-6;
	w_out.flags &= ~RF_DIFFUSE;
	pdf = 1e10;
	float f = lightMultiplier * pdf;
	brdf = Color(f, f, f);
}


Color Refraction::shade(const Ray& ray, const IntersectionInfo& info)
{
	extern Color raytrace(const Ray& ray, Scene& scene);
	
	float eta = ior;
	if (info.norm * ray.dir < 0)
		eta = 1.0f / eta;
	Vector n = faceforward(info.norm, ray.dir);
	Vector refr = refract(ray.dir, n, eta);
	if (refr.length() == 0)
		return Color(0, 0, 0);
	Ray newRay = ray;
	newRay.dir = refr;
	newRay.start = info.ip + ray.dir * 1e-6;
	newRay.depth = ray.depth + 1;
	return raytrace(newRay, *scene) * lightMultiplier;
}

Color Refraction::eval(const IntersectionInfo& x, const Ray& w_in, const Vector& w_out) const
{
	return Color(0, 0, 0);
}
void Refraction::spawnRay(const IntersectionInfo& x, const Ray& w_in, Ray& w_out, 
	Color& brdf, float& pdf) const
{
	float eta = ior;
	if (x.norm * w_in.dir < 0)
		eta = 1.0f / eta;
	Vector n = faceforward(x.norm, w_in.dir);
	Vector refr = refract(w_in.dir, n, eta);
	if (refr.length() == 0) {
		pdf = 0;
		return;
	}
	w_out = w_in;
	w_out.depth++;
	w_out.dir = refr;
	w_out.start = x.ip + w_in.dir * 1e-6;
	w_out.flags &= ~RF_DIFFUSE;

	pdf = 1e10;
	float f = lightMultiplier * pdf;
	brdf = Color(f, f, f);
}

void Layered::addLayer(Shader* shader, Color blend, Texture* blendTex)
{
	layers[nLayers].shader = shader;
	layers[nLayers].blend = blend;
	layers[nLayers].blendTex = blendTex;
	nLayers++;
}

void Layered::fillProperties(ParsedBlock& pb)
{
	char name[128];
	char value[256];
	int srcLine;
	for (int i = 0; i < pb.getBlockLines(); i++) {
		// fetch and parse all lines like "layer <shader>, <color>[, <texture>]"
		pb.getBlockLine(i, srcLine, name, value);
		if (!strcmp(name, "layer")) {
			char shaderName[200];
			char textureName[200] = "";
			bool err = false;
			if (!getFrontToken(value, shaderName)) {
				err = true;
			} else {
				stripPunctuation(shaderName);
			}
			if (!strlen(value)) err = true;
			if (!err && value[strlen(value) - 1] != ')') {
				if (!getLastToken(value, textureName)) {
					err = true;
				} else {
					stripPunctuation(textureName);
				}
			}
			if (!err && !strcmp(textureName, "NULL")) strcpy(textureName, "");
			Shader* shader = NULL;
			Texture* texture = NULL;
			if (!err) {
				shader = pb.getParser().findShaderByName(shaderName);
				err = (shader == NULL);
			}
			if (!err && strlen(textureName)) {
				texture = pb.getParser().findTextureByName(textureName);
				err = (texture == NULL);
			}
			if (err) throw SyntaxError(srcLine, "Expected a line like `layer <shader>, <color>[, <texture>]'");
			double x, y, z;
			get3Doubles(srcLine, value, x, y, z);
			addLayer(shader, Color((float) x, (float) y, (float) z), texture);
		}
	}
}

Color Layered::shade(const Ray& ray, const IntersectionInfo& info)
{
	Color result(0, 0, 0);
	for (int i = 0; i < nLayers; i++) {
		Layer& layer = layers[i];
		Color blend = layer.blendTex ? layer.blendTex->getTexColor(ray, info) : layer.blend;
		Color color = layer.shader->shade(ray, info);
		Color bottom = Color(1, 1, 1) - blend;
		result = result * bottom + color * blend;
	}
	return result;
}

static float fresnel(const Vector& i, const Vector& n, float eta)
{
	float f = sqr((1 - eta) / (1 + eta));
	float NdotI = fabs((float) dot(i, n));
	return f + (1 - f) * pow(1 - NdotI, 5.0f);
}

Color Fresnel::getTexColor(const Ray& ray, const IntersectionInfo& info)
{
	float eta = ior;
	if (ray.dir * info.norm > 0) eta = 1/eta;
	Vector n = faceforward(info.norm, ray.dir);
	float fr = fresnel(ray.dir, n, eta);
	return Color(fr, fr, fr);
}
