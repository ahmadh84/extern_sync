/*
Copyright (C) 2014 Jerome Revaud

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
*/
#ifndef ___IO_H___
#define ___IO_H___
#include <stdlib.h>

#include "image.h"
#include "deep_matching.h"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"


// output correspondences to a file or on the stdout
void output_correspondences( FILE* const fd, const corres_t* corres, 
					         int nb, float fx, float fy, const bool binary=false );

/* load a color image from a file */
color_image_t *color_image_load(const char *fname);

/* convert opencv mat to color_image_t */
color_image_t *cvmat_to_color_im(const cv::Mat& im_mat);

/* convert opencv mat to grayscale image_t (similar to image_gray_from_color()) 
   If a non-NULL im_t is given, it is used for the output grayscale image, rather 
   than allocating new memory for the image */
image_t *cvmat_to_gray_im(const cv::Mat& im_mat, image_t* im_t=NULL);

#endif
