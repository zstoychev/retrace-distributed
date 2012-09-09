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
#ifndef __BITMAP_H__
#define __BITMAP_H__

#include "color.h"
#include <stdio.h>

/// @brief a class that represents a bitmap (2d array of colors), e.g. a image
/// supports loading/saving to BMP
class Bitmap {
	int width, height;
	Color* data;
	bool differentiated;
	
	bool tryLoadPFM(FILE* f);
	
public:
	Bitmap(); //!< Generates an empty bitmap
	~Bitmap();
	void freeMem(void); //!< Deletes the memory, associated with the bitmap
	int getWidth(void) const; //!< Gets the width of the image (X-dimension)
	int getHeight(void) const; //!< Gets the height of the image (Y-dimension)
	bool isOK(void) const; //!< Returns true if the bitmap is valid
	void generateEmptyImage(int width, int height); //!< Creates an empty image with the given dimensions
	Color getPixel(int x, int y) const; //!< Gets the pixel at coordinates (x, y). Returns black if (x, y) is outside of the image
	Color getFilteredPixel(float x, float y) const; //!< Gets the filtered pixel at coordinates (x, y). Returns black if (x, y) is outside of the image
	void setPixel(int x, int y, const Color& col); //!< Sets the pixel at coordinates (x, y).
	
	bool loadBMP(const char* filename); //!< Loads an image from a BMP file. Returns false in the case of an error
	
	/// Saves the image to a BMP file (with clamping, etc). Returns false in the case of an error (e.g. read-only media)
	bool saveBMP(const char* filename, float outputGamma = 1.0f);
	
	bool loadPFM(const char* filename); //!< Loads an image from a .PFM file (a simplistic file format with floating-point data)
	
	bool loadImage(const char* filename); //!< Loads an image (type is detected by its extension)
	
	void differentiate(void);
};

#endif // __BITMAP_H__
