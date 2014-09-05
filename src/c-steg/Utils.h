//
//  Utils.h
//  c-steg
//

#ifndef CSTEG_UTILS_H
#define CSTEG_UTILS_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "../libjpeg/jpeglib.h"
#include "../libjpeg/jerror.h"

#ifdef __cplusplus
}
#endif

class DecodeResponse {
public:
	const unsigned int getLength() const {
		return length;
	}

	const unsigned char* const getData() const {
		return data;
	}

	~DecodeResponse() {
		free(data);
		data = NULL;
	}

private:
	unsigned int length;
	unsigned char* data;

	friend class StegUtils;
	friend class LossyDataFormatter;
};

class JpegCoefficients {
public:
    JBLOCKARRAY* coefBuffers;
    JBLOCKARRAY* rowPtrs;
    struct jpeg_decompress_struct* info;
    jvirt_barray_ptr *coefArrays;
    int numComponents;
    size_t* blockRowSize;
    int* widthInBlocks;
    int* heightInBlocks;

    /**
     * destructor for a coefficient_store class call it to clean up and deallocate memory
     */
    ~JpegCoefficients() {
    	jpeg_finish_decompress(info);
		jpeg_destroy_decompress(info);

		free(coefBuffers);
		free(rowPtrs);
		free(info);
		free(blockRowSize);
		free(widthInBlocks);
		free(heightInBlocks);

		coefBuffers = NULL;
		rowPtrs = NULL;
		info = NULL;
		coefArrays = NULL;
		blockRowSize = NULL;
		widthInBlocks = NULL;
		heightInBlocks = NULL;
    }
};

class Image {
public:
    unsigned char* rawImage;
    unsigned int width;
    unsigned int height;
    unsigned int numberOfComponents;
    J_COLOR_SPACE colorSpace;

    /**
     * constructor function for an image struct, this
     * function just zeroes out memory
     *
     * @param image
     *     the image struct to setup
     */
    Image() {
    	width = 0;
    	height = 0;
    	numberOfComponents = 0;
    	colorSpace = JCS_RGB;
    	rawImage = NULL;
    }

    /**
     * this function initializes and image struct with the parameters
     * from the jpeg_decompress_struct
     *
     * @param image
     *     the image struct to initialize
     *
     * @param info
     *     the jpge_decompress_struct to copy parameters from
     */
    void initImage(struct jpeg_decompress_struct* info) {
    	width = info->image_width;
    	height = info->image_height;
    	numberOfComponents = info->num_components;
    	colorSpace = info->jpeg_color_space;

    	if (rawImage == NULL) {
    		// allocate memory to hold the uncompressed image
    		rawImage = (unsigned char*) malloc(info->image_width * info->image_height * info->num_components);
    	}
    }

    /**
     * a destructor for an image class, call it to clean up and deallocate memory
     *
     * @param image
     *     the image struct to destroy
     */
    ~Image() {
    	free(rawImage);
    	rawImage = NULL;
    }
};

class Utils {
public:
	static void getFileNameWithoutExtension(const char* file, char* out_buffer);
	static short getLosslessEncodedShort(short coef, short coef_bits, short data_bits);
};

#endif
