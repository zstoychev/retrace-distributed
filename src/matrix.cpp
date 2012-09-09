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
#include "matrix.h"
#include "util.h"

Matrix rotationAroundX(double angle)
{
	double S = sin(angle);
	double C = cos(angle);
	Matrix a(1.0);
	a.m[1][1] = C;
	a.m[2][1] = S;
	a.m[1][2] = -S;
	a.m[2][2] = C;
	return a;
}

Matrix rotationAroundY(double angle)
{
	double S = sin(angle);
	double C = cos(angle);
	Matrix a(1.0);
	a.m[0][0] = C;
	a.m[2][0] = -S;
	a.m[0][2] = S;
	a.m[2][2] = C;
	return a;
}

Matrix rotationAroundZ(double angle)
{
	double S = sin(angle);
	double C = cos(angle);
	Matrix a(1.0);
	a.m[0][0] = C;
	a.m[1][0] = S;
	a.m[0][1] = -S;
	a.m[1][1] = C;
	return a;
}

Matrix operator * (const Matrix& a, const Matrix& b)
{
	Matrix c(0.0);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				c.m[i][j] += a.m[i][k] * b.m[k][j];
	return c;
}

double determinant(const Matrix& a)
{
	return a.m[0][0] * a.m[1][1] * a.m[2][2]
	     - a.m[0][0] * a.m[1][2] * a.m[2][1]
	     - a.m[0][1] * a.m[1][0] * a.m[2][2]
	     + a.m[0][1] * a.m[1][2] * a.m[2][0]
	     + a.m[0][2] * a.m[1][0] * a.m[2][1]
	     - a.m[0][2] * a.m[1][1] * a.m[2][0];
}

double cofactor(const Matrix& m, int ii, int jj)
{
	int rows[2], rc = 0, cols[2], cc = 0;
	for (int i = 0; i < 3; i++)
		if (i != ii) rows[rc++] = i;
	for (int j = 0; j < 3; j++)
		if (j != jj) cols[cc++] = j;
	double t = m.m[rows[0]][cols[0]] * m.m[rows[1]][cols[1]] - m.m[rows[1]][cols[0]] * m.m[rows[0]][cols[1]];
	if ((ii + jj) % 2) t = -t;
	return t;
}

Matrix inverseMatrix(const Matrix& m)
{
	double D = determinant(m);
	if (fabs(D) < 1e-12) return m; // an error; matrix is not invertible
	double rD = 1.0 / D;
	Matrix result;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			result.m[i][j] = rD * cofactor(m, j, i);
	return result;
}

Matrix transposeMatrix(const Matrix& a)
{
	Matrix result;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			result.m[i][j] = a.m[j][i];
	return result;
}

Transform::Transform()
{
	resetTransform();
}
void Transform::resetTransform(void)
{
	reversed = false;
	offset.makeZero();
	m = Matrix(1);
}

void Transform::scale(double sx, double sy, double sz)
{
	Matrix scaleM(1);
	scaleM.m[0][0] = sx;
	scaleM.m[1][1] = sy;
	scaleM.m[2][2] = sz;
	// sx  0  0
	//  0 sy  0
	//  0  0 sz
	m = m * scaleM;
}
void Transform::rotate(double yaw, double pitch, double roll)
{
	m = m * rotationAroundZ(toRadians(roll))
	      * rotationAroundX(toRadians(pitch))
	      * rotationAroundY(toRadians(yaw));
}
void Transform::translate(const Vector& t)
{
	offset += t;
}

Transform Transform::inverseTransform(void) const
{
	Transform result(*this);
	result.reversed = !reversed;
	result.m = inverseMatrix(m);
	result.offset = -offset;
	return result;
}

Vector Transform::transformPoint(const Vector& p) const
{
	if (!reversed)
		return (p * m) + offset;
	else
		return (p + offset) * m;
}

Vector Transform::transformDir(const Vector& dir) const
{
	return dir * m;
}

Ray Transform::transformRay(Ray ray) const
{
	ray.start = transformPoint(ray.start);
	ray.dir = transformDir(ray.dir);
	ray.dir.normalize();
	return ray;
}

