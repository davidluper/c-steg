//
//  jpeg_tools.cpp
//	c-steg
//
//  this is a layer on top of the The Independent JPEG Group's JPEG software stack in the jpeg folder
//  i did not write The Independent JPEG Group's JPEG software stack which i think is also referred to
//  as libjpeg i am merely interfacing into it.  thanks for providing it guys!

#include "JpegUtils.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <stdio.h>

/**
 * this function reads a jpeg image from disk and stores it in memory
 * in the provided image class
 *
 * @param image
 *     the image class to store the uncompressed image in
 * @param imageFileName
 *     the path to an image on disk
 *
 * returns true everything went ok, false if some error occured
 */
int JpegUtils::readJpegFile(Image& image, const char* imageFileName) {
	//these are standard libjpeg structures for decompression
	struct jpeg_decompress_struct info;
	struct jpeg_error_mgr jerr;
	// libjpeg data structure for storing one row, that is, scanline of an image
	JSAMPROW row_pointer[1];

	FILE* image_file = fopen(imageFileName, "rb");
	int i = 0;

	if (!image_file) {
		printf("error opening jpeg file %s\n!", imageFileName);
		return 0;
	}

	// here we set up the standard libjpeg error handler
	info.err = jpeg_std_error(&jerr);
	// setup decompression process and source, then read JPEG header
	jpeg_create_decompress(&info);
	// this makes the library read from image_file
	jpeg_stdio_src(&info, image_file);
	// reading the image header which contains image information
	jpeg_read_header(&info, TRUE);

	jpeg_start_decompress(&info);

	image.initImage(&info);

	// now actually read the jpeg into the raw buffer
	row_pointer[0] = (unsigned char *) malloc(info.image_width * info.num_components);

	unsigned long location = 0;
	// read one scan line at a time
	while (info.output_scanline < info.image_height) {
		jpeg_read_scanlines(&info, row_pointer, 1);
		for (i = 0; i < info.image_width * info.num_components; ++i) {
			image.rawImage[location++] = row_pointer[0][i];
		}
	}

	// wrap up decompression, destroy objects, free pointers and close open files
	jpeg_finish_decompress(&info);
	jpeg_destroy_decompress(&info);
	free(row_pointer[0]);
	fclose(image_file);

	return 1;
}

/**
 * this function writes an uncompressed image from an image class to disk as a
 * compressed jpeg file
 *
 * @param image
 *     the image class that holds the uncompressed image
 * @param outputFileName
 *     the path to save the image to on disk
 * @param quality
 *     the quality to save the compressed jpeg with (0 to 100)
 *
 * returns true if everything went as planned, false if an error occurred.
 */
int JpegUtils::writeJpegFile(Image& image, const char* outputFileName, int quality) {
	if (quality < 1) {
		quality = 1;
	}

	if (quality > 100) {
		quality = 100;
	}

	struct jpeg_compress_struct info;
	struct jpeg_error_mgr jerr;

	// this is a pointer to one row of image data
	JSAMPROW row_pointer[1];
	FILE* image_file = fopen(outputFileName, "wb");

	if (!image_file) {
		printf("error opening output jpeg file %s\n!", outputFileName);
		return 0;
	}

	info.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&info);
	jpeg_stdio_dest(&info, image_file);

	// Setting the parameters of the output file here
	info.image_width = image.width;
	info.image_height = image.height;
	info.input_components = image.numberOfComponents;
	info.in_color_space = JCS_RGB;

	// default compression parameters, we shouldn't be worried about these
	jpeg_set_defaults(&info);

	info.num_components = 3;
	//info.data_precision = 4;
	info.dct_method = JDCT_ISLOW;
	//info.dct_method = JDCT_IFAST;
	jpeg_set_quality(&info, quality, TRUE);

	// Now do the compression ..
	jpeg_start_compress(&info, TRUE);
	// like reading a file, this time write one row at a time
	while (info.next_scanline < info.image_height) {
		row_pointer[0] = &(image.rawImage)[info.next_scanline
				* info.image_width * info.input_components];

		jpeg_write_scanlines(&info, row_pointer, 1);
	}

	// similar to read file, clean up after we're done compressing
	jpeg_finish_compress(&info);
	jpeg_destroy_compress(&info);
	fclose(image_file);

	return 1;
}

/**
 * this function writes a jpeg to disk from the supplied jpeg coefficient object
 *
 * @param coefficients
 *     the image class that holds the uncompressed image
 * @param outputFileName
 *     the path to save the image to on disk
 * @param quality
 *     the quality to save the compressed jpeg with (0 to 100)
 *
 * returns true if everything went as planned, false if an error occurred.
 */
bool JpegUtils::writeJpegCoefficients(JpegCoefficients& coefficients, const char* outputFileName, int quality) {
	FILE* output_file = fopen(outputFileName, "wb");
	if (!output_file) {
		printf("error opening jpeg file %s\n!", outputFileName);
		return false;
	}

	// these are standard libjpeg structures for reading(decompression)
	struct jpeg_compress_struct output_info;
	struct jpeg_error_mgr jerr;

	output_info.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&output_info);
	jpeg_stdio_dest(&output_info, output_file);

	// the coefficients are more stable with SLOW method
	// meaning when you uncompress then recompress they
	// seems to change less
	output_info.dct_method = JDCT_ISLOW;
	//info.dct_method = JDCT_IFAST;
	jpeg_set_quality(&output_info, quality, TRUE);
	// copy compression parameters from the input file to the output file
	jpeg_copy_critical_parameters(coefficients.info, &output_info);

	// output the new DCT coeffs to a JPEG file
	for (JDIMENSION compnum = 0; compnum < coefficients.numComponents; ++compnum) {
		for (JDIMENSION rownum = 0; rownum < coefficients.heightInBlocks[compnum]; ++rownum) {
			coefficients.rowPtrs[compnum] = ((&output_info)->mem->access_virt_barray)(
					(j_common_ptr) &output_info, coefficients.coefArrays[compnum], rownum, (JDIMENSION) 1, TRUE);
			memcpy(coefficients.rowPtrs[compnum][0][0], coefficients.coefBuffers[compnum][rownum][0], coefficients.blockRowSize[compnum]);
		}
	}

	// write to the output file
	jpeg_write_coefficients(&output_info, coefficients.coefArrays);

	// finish compression and release memory
	jpeg_finish_compress(&output_info);
	jpeg_destroy_compress(&output_info);

	// close files
	fclose(output_file);

	return true;
}

/**
 * this function loads the coefficients from a jpeg image into a coefficient store struct
 *
 * @param inputFileName
 *     the full path/name of the jpeg file on disk for which to load coefficients from
 * @param coefficients
 * 	   the jpeg coefficients object to store the read coefficients in
 *
 * returns true on success false otherwise
 */
bool JpegUtils::readJpegCoefficients(const char *inputFileName, JpegCoefficients& coefficients) {
	struct jpeg_error_mgr jerr;
	JDIMENSION i, rownum, blocknum;
	JBLOCKARRAY* coef_buffers = new JBLOCKARRAY[MAX_COMPONENTS];
	JBLOCKARRAY* row_ptrs = new JBLOCKARRAY[MAX_COMPONENTS];

	/* Open the input and output files */
	FILE* jpeg_file = fopen(inputFileName, "rb");
	if (!jpeg_file) {
		printf("error opening jpeg file %s\n!", inputFileName);
		return false;
	}

	// initialize the JPEG compression and decompression objects with default error handling
	struct jpeg_decompress_struct* info =
			(struct jpeg_decompress_struct*) malloc(sizeof(struct jpeg_decompress_struct));
	(*info).err = jpeg_std_error(&jerr);
	jpeg_create_decompress(info);

	// specify data source for decompression
	jpeg_stdio_src(info, jpeg_file);

	// read file header
	(void) jpeg_read_header(info, TRUE);

	// allocate memory for reading out DCT coeffs
	for (JDIMENSION compnum = 0; compnum < (*info).num_components; ++compnum) {
		coef_buffers[compnum] = (info->mem->alloc_barray)((j_common_ptr) info,
				JPOOL_IMAGE, (*info).comp_info[compnum].width_in_blocks,
				(*info).comp_info[compnum].height_in_blocks);
	}

	// read input file as DCT coeffs
	jvirt_barray_ptr *coef_arrays = jpeg_read_coefficients(info);

	// copy DCT coeffs to a new array
	int num_components = (*info).num_components;

	size_t* block_row_size = new size_t[num_components];
	int* width_in_blocks = new int[num_components];
	int* height_in_blocks = new int[num_components];
	for (JDIMENSION compnum = 0; compnum < num_components; ++compnum) {
		height_in_blocks[compnum] = (*info).comp_info[compnum].height_in_blocks;
		width_in_blocks[compnum] = (*info).comp_info[compnum].width_in_blocks;
		block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
		for (rownum = 0; rownum < height_in_blocks[compnum]; rownum++) {
			row_ptrs[compnum] = (info->mem->access_virt_barray)(
					(j_common_ptr) &info, coef_arrays[compnum], rownum,
					(JDIMENSION) 1, FALSE);
			for (blocknum = 0; blocknum < width_in_blocks[compnum]; blocknum++) {
				for (i = 0; i < DCTSIZE2; i++) {
					coef_buffers[compnum][rownum][blocknum][i] =
							row_ptrs[compnum][0][blocknum][i];
				}
			}
		}
	}

	fclose(jpeg_file);

	coefficients.coefBuffers = coef_buffers;
	coefficients.rowPtrs = row_ptrs;
	coefficients.info = info;
	coefficients.numComponents = num_components;
	coefficients.blockRowSize = block_row_size;
	coefficients.widthInBlocks = width_in_blocks;
	coefficients.heightInBlocks = height_in_blocks;
	coefficients.coefArrays = coef_arrays;

	return true;
}
