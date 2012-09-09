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

#include "camera.h"
#include "matrix.h"
#include "sdl.h"
#include "random_generator.h"

Camera::Camera()
{
	yaw = pitch = roll = 0;
	aspect = 4.0/3.0;
	fov = 90;
	pos.makeZero();
	stereoSeparation = 0;
	leftMask = Color(1.0f, 0.0f, 0.0f);
	rightMask = Color(0.0f, 1.0f, 1.0f);
	samples = 0;
}

void Camera::fillProperties(ParsedBlock& pb)
{
	dof = false;
	autofocus = true;
	pb.getBoolProp("dof", &dof);
	pb.getBoolProp("autofocus", &autofocus);
	pb.getDoubleProp("stereoSeparation", &stereoSeparation);
	pb.getColorProp("leftMask", &leftMask);
	pb.getColorProp("rightMask", &rightMask);
	pb.getDoubleProp("fNumber", &fNumber);
	pb.getIntProp("samples", &samples);
	pb.getDoubleProp("focusDist", &focusDist);
	pb.getDoubleProp("aspect", &aspect, 1e-6);
	pb.getDoubleProp("fov", &fov, 0.0001, 179);
	if (!pb.getVectorProp("pos", &pos))
		pb.requiredProp("pos");
	bool hasYaw = pb.getDoubleProp("yaw", &yaw);
	bool hasPitch = pb.getDoubleProp("pitch", &pitch, -90, 90);
	bool hasRoll = pb.getDoubleProp("roll", &roll);
	bool yawPitchRoll = hasYaw || hasPitch || hasRoll;
	// check if lookAt / upDir are present in the parsed block:
	Vector lookAt, upDir = Vector(0, 1, 0);
	if (pb.getVectorProp("upDir", &upDir)) {
		if (!pb.getVectorProp("lookAt", &lookAt))
			pb.requiredProp("lookAt");
	}
	if (pb.getVectorProp("lookAt", &lookAt)) {
		if (yawPitchRoll) {
			// user has specified yaw/pitch/roll and lookAt at the same time:
			pb.signalError("The yaw/pitch/roll and the lookAt/upDir systems cannot be used in the same time!");
		} else {
			// check if the two given directions (the lookAt direction and the upDir direction) aren't colinear:
			if (((lookAt - pos) ^ upDir).length() < 1e-6) pb.signalError("Gimbal lock: upDir is collinear to lookAt direction");
			setPosition(pos, lookAt, upDir);
		}
	}
}

void Camera::beginFrame(void)
{
	double x, y;
	x = -aspect;
	y = 1;
	
	double curLength = sqrt(x*x + y*y);
	double L = tan(toRadians(fov / 2));
	double factor = L / curLength;
	x *= factor;
	y *= factor;
	
	upLeft = Vector(x, y, 1);
	upRight = Vector(-x, y, 1);
	downLeft = Vector(x, -y, 1);
	
	Matrix mroll = rotationAroundZ(toRadians(roll));
	Matrix mpitch = rotationAroundX(toRadians(pitch));
	Matrix myaw = rotationAroundY(toRadians(yaw));
	
	Matrix all = mroll * mpitch * myaw;
	
	upLeft *= all;
	upRight *= all;
	downLeft *= all;
	
	frontDir = Vector(0, 0, 1);
	upDir = Vector(0, 1, 0);
	rightDir = Vector(1, 0, 0);
	
	frontDir *= all;
	upDir *= all;
	rightDir *= all;
	
	radius = 1.0 / fNumber;
	
	if (autofocus && dof) {
		extern double getClosestGeom(Ray ray, Scene& scene);
		focusDist = getClosestGeom(Ray(this->pos, frontDir), *scene);
		printf("autofocus: focus distance is %.1lf units\n", focusDist);
	}
	
	upLeft = upLeft + pos;
	upRight = upRight + pos;
	downLeft = downLeft + pos;
}

Ray Camera::getScreenRay(double x, double y, int offset)
{
	Vector dest = upLeft + (upRight - upLeft) * (x / scene->settings.frameWidth) 
	                     + (downLeft - upLeft) * (y / scene->settings.frameHeight);
	Ray result;
	result.start = pos;
	result.dir = dest - pos;
	result.dir.normalize();
	if (offset != 0) {
		// move the ray start, for stereoscopic rendering:
		result.start += rightDir * (stereoSeparation * 0.5 * offset);
	}
	if (!dof) return result;
	//
	double cosTheta = frontDir * result.dir;
	double M = focusDist / cosTheta;
	Vector T = this->pos + result.dir * M;
	Random& rnd = getRandomGen();
	double u, v;
	rnd.unitDiscSample(u, v);
	u *= radius;
	v *= radius;
	result.start = this->pos + u * rightDir + v * upDir;
	result.dir = T - result.start;
	result.dir.normalize();
	return result;
}

Ray Camera::getCenterRay() const {
	Vector dest = upLeft + (upRight - upLeft) * 0.5
	                     + (downLeft - upLeft) * 0.5;
	Ray result;
	result.start = pos;
	result.dir = dest - pos;
	result.dir.normalize();
	return result;
}

void Camera::setPosition(Vector pos, Vector lookAt, Vector upDir)
{
	this->pos = pos;
	Vector frontDir = lookAt - pos;
	frontDir.normalize();
	upDir.normalize();
	//
	// create a version of upDir, which is "canonic" - i.e., truly orthogonal to frontDir
	Vector rightDir = upDir ^ frontDir;
	upDir = frontDir ^ rightDir;
	// calculate yaw and pitch from frontDir alone
	yaw = toDegrees(atan2(frontDir.z, frontDir.x)) - 90; // yaw = 0 -> frontDir.z = 1; we compensate for this angle here.
	pitch = toDegrees(asin(frontDir.y));
	// get the canonic upDir (positive Y) and apply pitch and yaw to it. We'll get a vector, which lies
	// in the plane, determined by upDir and rightDir, and coincides with the upDir when pitch is 0.
	// Then express this vector in terms of rightDir and upDir. The polar angle of that 2D point is our
	// roll angle.
	Vector upAfterYawPitch = Vector(0, 1, 0) * rotationAroundX(toRadians(pitch)) * rotationAroundY(toRadians(yaw));
	double fx = upAfterYawPitch * rightDir;
	double fy = upAfterYawPitch * upDir;
	roll = toDegrees(atan2(fy, fx)) - 90; // the canonic updir would have fx = 0, fy = 1, thus atan() will be 90 degrees.
}

void Camera::move(double ff, double ss)
{
	pos += rightDir * ss + frontDir * ff;
}

void Camera::rotate(double xx, double yy)
{
	yaw -= xx / 3.14 * 180;
	pitch += yy / 3.14 * 180;
	if (pitch > 90) pitch = 90;
	if (pitch < -90) pitch = -90;
}

