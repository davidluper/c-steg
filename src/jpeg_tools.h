//
//  jpeg_tools.h
//  c-steg
//
//

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "lossy_data_formatter.h"
#ifdef __cplusplus
extern "C" {
#endif

#include "jpeg/jpeglib.h"
#include "jpeg/jerror.h"

#ifdef __cplusplus
}
#endif

struct coefficient_store {
    JBLOCKARRAY* coefficient_buffer;
    struct jpeg_decompress_struct* info;
};

struct image {
    unsigned char* raw_image;
    unsigned int width;
    unsigned int height;
    unsigned int number_of_components;
    J_COLOR_SPACE color_space;
};

#define CSTEG_JPEG_VERSION 1
#define CSTEG_JPEG_VERSION_1_HEADER_LENGTH 51 //16 + 32 + 3
#define NOISE_MODULO_VALUE 5
#define VERSION_BIT_LENGTH 16
#define DATA_LENGTH_BIT_LENGTH 32
#define BITS_TO_STEAL_BIT_LENGTH 3
#define MAGIC_XOR 1563965293
#define MAGIC_XOR_EVEN 139
#define MAGIC_XOR_ODD 242
#define COEFFICIENT_MAX 32765 //give some room to spare (actual max should be 32767)
#define COEFFICIENT_MIN -32765 //give some room to spare (actual in should be -32767)

unsigned int lossless_jpeg_file_get_number_of_coefficients(const char *input_file_name);
void get_without_extension(const char* file, char* out_buffer);
void construct_image(struct image* image);
void init_image(struct image* image, struct jpeg_decompress_struct* info);
void destruct_image(struct image* image);
void destruct_coefficient_store(struct coefficient_store* coefficient_store);
int read_jpeg_file(struct image* image, const char* image_file_name);
int write_jpeg_file(struct image* image, const char *image_file_name, int quality);
int encode_lossless_jpeg_file(const char *input_file_name, const char* output_file_name,
                               unsigned char* data, int data_length,
                               unsigned short bits_to_steal, int quality);
int decode_lossless_jpeg_file(const char *input_file_name, unsigned char** data_out, unsigned int* data_out_length);
int encode_lossy_jpeg_file(const char *input_file_name, const char* output_file_name, struct packet_collection* packet_collection, int quality);
int decode_lossy_jpeg_file(const char *input_file_name, struct packet_collection* packet_collection);
int diff_coefficients(struct coefficient_store* _1, struct coefficient_store* _2);
struct coefficient_store* load_coefficients(char *input_file_name_1);
void do_compressions(int n, const char* file_path, int quality);
int do_coefficient_diff(const char* file_path_1, const char* file_path_2);
void do_compression_series_coefficient_diffs(int n, const char* file_path);
void purge(const char* jpeg_file_path, int limit);
void purge(const char* directory, const char* jpeg_file_name, int limit);
void move_quality_result(int limit, int quality, const char* directory, const char* jpeg_file_name);
void make_quality_panel(int bottom, int top, int number_of_compressions, const char* directory, const char* jpeg_file_name);

#endif
