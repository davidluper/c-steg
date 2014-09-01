//  c_steg_api.h
//  c-steg
//
//  this header file defines what is meant to be the externally visible api
//  for the c-steg library
//

#ifndef C_STEG_API_H
#define C_STEG_API_H

#ifdef __cplusplus
extern "C" {
#endif

int c_steg_encode_lossy_jpeg(const char* input_file, const char* output_file, const char* data,
                              unsigned int data_length, int quality, int number_of_compressions);

int c_steg_decode_lossy_jpeg(const char* input_file, unsigned char** data_out, unsigned int* data_out_length);

int c_steg_encode_lossless_jpeg(const char *input_file_name, const char* output_file_name,
                               unsigned char* data, int data_length,
                               unsigned short bits_to_steal, int quality);

int c_steg_decode_lossless_jpeg(const char *input_file_name, unsigned char** data_out, unsigned int* data_out_length);

unsigned int c_steg_lossless_jpeg_file_get_number_of_coefficients(const char* input_file_name);

#ifdef __cplusplus
}
#endif

#endif
