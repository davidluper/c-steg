//
//  c_steg_api.cpp
//  c-steg
//
//  this file contains what is meant to be the external api for
//  the c-steg library
//  

#include "c_steg_api.h"
#include "jpeg_tools.h"
#include "lossy_data_formatter.h"

/**
 * this function exposes the ability to encode a lossy byte payload into a jpeg image.  the
 * lossy payload is meant to survive recompressions
 *
 * @param input_file
 *     the full path/name to the jpeg to encode information into
 * @param output_file
 *     the full path/name to the output image that contains the encoded information
 * @param data
 *     a pointer to the data to encode into the image
 * @param data_length
 *     the length of the data in the data buffer of information to encode
 * @param quality
 *     the quality number to use when saving the image with the encoded information
 * @param number_of_compressions
 *     the number of times to uncompress and recompress the image to stabilize the DCT coefficients
 *
 * returns 1 if successful 0 otherwise
 */
int c_steg_encode_lossy_jpeg(const char* input_file, const char* output_file, const char* data,
                              unsigned int data_length, int quality, int number_of_compressions) {
    do_compressions(number_of_compressions, input_file, quality);
    purge(input_file, number_of_compressions - 1);

    char file_name_buffer[256];
    memset(file_name_buffer, 0, 256);
    get_without_extension(input_file, file_name_buffer);

    char file_path_buffer[256];
    memset(file_path_buffer, 0, 256);
    int i = number_of_compressions - 1;
    sprintf(file_path_buffer, "%s-%d.jpeg", file_name_buffer, i);

    struct packet_collection* to_packet_data = to_lossy(data, data_length);
    int rtn = encode_lossy_jpeg_file(file_path_buffer, output_file, to_packet_data, quality);

    remove(file_path_buffer);
    free(to_packet_data->data);
    free(to_packet_data);

    return rtn;
}

/*
 * this function exposes the ability to decode a lossy byte payload from a jpeg that was encoded
 * with lossy formatted information
 *
 * @param input_file
 *     the full path/name to the jpeg image to decode information from
 * @param data_out
 *     a pointer that if successful will point to the decoded information
 * @param data_out_length
 *     a pointer to an unsigned int that if successful will be set to the 
 *     length of the decoded information in data_out
 *
 * returns 1 if successful 0 otherwise
 */
int c_steg_decode_lossy_jpeg(const char* input_file, unsigned char** data_out, unsigned int* data_out_length) {
	char* data = NULL;
    unsigned int data_length = 0;
    struct packet_collection from_packet_data;
    int rtn = decode_lossy_jpeg_file(input_file, &from_packet_data);

    if (rtn) {
        *data_out = NULL;
        from_lossy(&from_packet_data, data_out, data_out_length);
    }
    
    free(from_packet_data.data);
    
    return rtn;
}

/**
 * this function exposes the ability to encode lossless information in a jpeg.  its assumed 
 * when decoding the image it has not been recompressed.
 *
 * @param input_file
 *     the full path/name to the jpeg to encode information into
 * @param output_file
 *     the full path/name to the output image that contains the encoded information
 * @param data
 *     a pointer to the data to encode into the image
 * @param data_length
 *     the length of the data in the data buffer of information to encode
 * @param bits_to_steal
 *     the number of bits to steal per coefficient when embedding information
 * @param quality
 *     the quality number to use when saving the image with the encoded information
 * 
 * returns 1 if successful 0 otherwise
 */
int c_steg_encode_lossless_jpeg(const char *input_file_name, const char* output_file_name,
                               unsigned char* data, int data_length,
                               unsigned short bits_to_steal, int quality) {
    return encode_lossless_jpeg_file(input_file_name, output_file_name, data, data_length, bits_to_steal, quality);
}

/**
 * this function exposes the ability to decode information that was encoded with the lossless 
 * format.  it is assumed the image has not been recompressed.
 *
 * @param input_file
 *     the full path/name to the jpeg image to decode information from
 * @param data_out
 *     a pointer that if successful will point to the decoded information
 * @param data_out_length
 *     a pointer to an unsigned int that if successful will be set to the 
 *     length of the decoded information in data_out
 *
 * returns 1 if successful 0 otherwise
 */
int c_steg_decode_lossless_jpeg(const char *input_file_name, unsigned char** data_out, unsigned int* data_out_length) {
    return decode_lossless_jpeg_file(input_file_name, data_out, data_out_length);
}

/*
 * this function calculates the number of useable coefficients for encoding 
 * information into an image.
 *
 * @param input_file_name
 *     the full name/path to the jpeg to calculate the number of useable coefficient for
 *
 * returns the number of useable coefficients for embedding information
 */
unsigned int c_steg_lossless_jpeg_file_get_number_of_coefficients(const char *input_file_name) {
    return lossless_jpeg_file_get_number_of_coefficients(input_file_name);
}
