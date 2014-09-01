//
//  jpeg_tools.cpp
//  c-steg
//

#include "jpeg_tools.h"
#include <math.h>
#include <time.h>

/**
 * this function removes a file extension from a path and stores the 
 * result in out_buffer
 *
 * @param file
 *     path to a file to remove the extension from
 * @param out_buffer
 *     char buffer to store the copy of the path to the file without
 *     the extension
 */
void get_without_extension(const char* file, char* out_buffer) {
    int file_name_len = strlen(file);
    memset(out_buffer, 0, file_name_len + 1);
    strcat(out_buffer, file);
    for (int i = file_name_len - 1; i > -1; --i) {
        char c = out_buffer[i];
        out_buffer[i] = 0x00;
        if (c == '.') {
            break;
        }
    }
}

/**
 * constructor function for an image struct, this 
 * function just zeroes out memory
 *
 * @param image
 *     the image struct to setup 
 */
void construct_image(struct image* image) {
    image->width = 0;
    image->height = 0;
    image->number_of_components = 0;
    image->color_space = JCS_RGB;
    image->raw_image = NULL;
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
void init_image(struct image* image, struct jpeg_decompress_struct* info) {
    image->width = info->image_width;
    image->height = info->image_height;
    image->number_of_components = info->num_components;
    image->color_space = info->jpeg_color_space;

    if (image->raw_image == NULL) {
        // allocate memory to hold the uncompressed image
        image->raw_image = (unsigned char*)
            malloc(
                info->image_width *
                info->image_height *
                info->num_components
            );
    }
}

/**
 * this function acts as a destructor for an image struct
 * call it to clean up and deallocate memory
 *
 * @param image
 *     the image struct to destroy
 */
void destruct_image(struct image* image) {
    free(image->raw_image);
    image->raw_image = NULL;
}

/**
 * this function acts as a destructor for a coefficient_store struct
 * call it to clean up and deallocate memory
 *
 * @param coefficient_store
 *     the coefficient_store struct to destroy
 */
void destruct_coefficient_store(struct coefficient_store* coefficient_store) {
    jpeg_finish_decompress(coefficient_store->info);
    jpeg_destroy_decompress(coefficient_store->info);

    free(coefficient_store->coefficient_buffer);
    free(coefficient_store->info);
}

/**
 * this function reads a jpeg image from disk and stores it in memory
 * in the provided image struct
 *
 * @param image
 *     the image struct to store the uncompressed image in
 * @param image_file_name
 *     the path to an image on disk
 *
 * returns true everything went ok, false if some error occured
 */
int read_jpeg_file(struct image* image, const char* image_file_name) {
    //these are standard libjpeg structures for decompression
    struct jpeg_decompress_struct info;
    struct jpeg_error_mgr jerr;
    // libjpeg data structure for storing one row, that is, scanline of an image
    JSAMPROW row_pointer[1];
    
    FILE* image_file = fopen(image_file_name, "rb");
    int i = 0;
    
    if (!image_file) {
        printf("error opening jpeg file %s\n!", image_file_name);
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

    init_image(image, &info);
    
    // now actually read the jpeg into the raw buffer
    row_pointer[0] =
    	(unsigned char *) malloc(
    		info.image_width *
            info.num_components
        );

    unsigned long location = 0;
    // read one scan line at a time
    while(info.output_scanline < info.image_height) {
        jpeg_read_scanlines(&info, row_pointer, 1);
        for(i = 0; i < info.image_width * info.num_components; ++i) {
            image->raw_image[location++] = row_pointer[0][i];
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
 * this function writes an uncompressed image from an image struct to disk as a 
 * compressed jpeg file
 *
 * @param image
 *     the image struct that holds the uncompressed image
 * @param file_name
 *     the path to save the iamge to on disk
 * @param quality
 *     the quality to save the compressed jpeg with (0 to 100)
 *
 * returns true if everything went as planned, false if an error occurred.
 */
int write_jpeg_file(struct image* image, const char* file_name, int quality) {
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
    FILE* image_file = fopen(file_name, "wb");
    
    if (!image_file) {
        printf("error opening output jpeg file %s\n!", file_name);
        return 0;
    }

    info.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&info);
    jpeg_stdio_dest(&info, image_file);
    
    // Setting the parameters of the output file here
    info.image_width = image->width;
    info.image_height = image->height;
    info.input_components = image->number_of_components;
    info.in_color_space = JCS_RGB;

    // default compression parameters, we shouldn't be worried about these    
    jpeg_set_defaults(&info);

    info.num_components = 3;
    //info.data_precision = 4;
    info.dct_method = JDCT_ISLOW;
    //info.dct_method = JDCT_IFAST;
    jpeg_set_quality(&info, quality, TRUE);

    // Now do the compression ..
    jpeg_start_compress( &info, TRUE );
    // like reading a file, this time write one row at a time
    while(info.next_scanline < info.image_height) {
        row_pointer[0] = &(image->raw_image)[
            info.next_scanline * 
            info.image_width * 
            info.input_components
        ];

        jpeg_write_scanlines(&info, row_pointer, 1);
    }

    // similar to read file, clean up after we're done compressing
    jpeg_finish_compress(&info);
    jpeg_destroy_compress(&info);
    fclose(image_file);

    return 1;
}

/**
 * this function is used in the lossless encoding function to modify the coefficient
 * to hold a single bit from the data. this is used for encoding the header of the lossless
 * encoder
 *
 * @param coef
 *     the coefficient that need to hold the information from the data
 * @ param coef_bits
 *     the bits from the coefficient to compare against the data_bits to
 *     see if the coefficient needs to be modified
 * @param data_bits
 *     the bits from the data that need to be in the coefficient
 *
 * returns the coefficient modified (if necessary) to hold the information 
 * from the data
 */
short get_lossless_encode_bit(short coef, short coef_bits, short data_bits) {
    // this will always be a single bit, if the bit from data and coef already
    // equal each other then nothing is necessary
    if (coef_bits != data_bits) {
        if (coef >= COEFFICIENT_MAX) {
            //be careful to never over flow
            --coef;
        } else if (coef <= COEFFICIENT_MIN) {
            ++coef;
        } else {
            //if the bits are different then modify the coefficient bit
            //to hold the information from the data
            switch (coef) {
            case 1:
                //never let a coefficient that was not already zero turn into a zero
                ++coef;
                break;
            case -1:
                //never let a coefficient that was not already zero turn into a zero
                --coef;
                break;
            default:
                //randomly either go up one or down one to modifiy the coefficient bit
                if (rand() % 2 == 1) {
                    ++coef;
                } else {
                    --coef;
                }
            }
        }
    }

    return coef;
}

/**
 * this function uses LSB encoding to embed a data payload into the DCT coefficients of 
 * a jpeg image.  this functions assumes the decode lossless function will be run on
 * the exact image produced as ouput of this routine.  if the output image from this
 * routine is altered in any way decoding will fail.  if the image is expected to be
 * altered either by cropping or recompression you should use the lossy encode/decode
 * functions
 *
 * @param input_file_name 
 *     the complete file path/name to the image to embed information into
 * @param output_file_name
 *     the complete file path/name to the output image that will contain
 *     the image with embedded information
 * @param data
 *     the data (presented as an unsigned char array) to embed into the coefficient
 *     of the specified jpeg
 * @param data_length
 *     the number of unsigned char bytes to encode from data (the unsigned 
 *     char buffer input)
 * @param bit_to_steal
 *     the number of least significant bits to steal from a coefficient and replace 
 *     with data bits
 * @param quality
 *     the quality number to use when saving the jpeg with embeded information
 *
 * returns true if everything went as planned, false if an error occurred
 */
int encode_lossless_jpeg_file(const char *input_file_name, const char* output_file_name,
                               unsigned char* data, int data_length, 
                               unsigned short bits_to_steal, int quality) {
    for (int i = 0; i < data_length; ++i) {
        if (i % 2 == 1) {
            //odd
            data[i] = data[i] ^ MAGIC_XOR_ODD;
        } else {
            //even
            data[i] = data[i] ^ MAGIC_XOR_EVEN;
        }
    }
    
    // these are standard libjpeg structures for reading(decompression)
    struct jpeg_decompress_struct input_info;
    struct jpeg_compress_struct output_info;
    struct jpeg_error_mgr jerr;
    JBLOCKARRAY coef_buffers[MAX_COMPONENTS];
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    // libjpeg data structure for storing one row, that is, scanline of an image
    
    // open the input and output files
    FILE* input_file = fopen(input_file_name, "rb");
    if ( !input_file ) {
        printf("error opening jpeg file %s\n!", input_file_name);
        return 0;
    }
    
    FILE* output_file = fopen(output_file_name, "wb");
    if (!output_file) {
        printf("error opening jpeg file %s\n!", output_file_name);
        return 0;
    }

    // initialize the JPEG compression and decompression objects with default error handling.
    input_info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&input_info);
    output_info.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&output_info);
    
    // specify data source for decompression and recompression
    jpeg_stdio_src(&input_info, input_file);
    jpeg_stdio_dest(&output_info, output_file);

    // default compression parameters, we shouldn't be worried about these    
    //jpeg_set_defaults(&output_info);
    //output_info.num_components = 3;
    //info.data_precision = 4;

    // make sure to use JDCT_ISLOW
    // it helps with facebook compression series
    // this method only requires 5 or 6 
    // compressions/decompressions for the coefficients
    // to become stable
    output_info.dct_method = JDCT_ISLOW;
    //info.dct_method = JDCT_IFAST;
    jpeg_set_quality(&output_info, quality, TRUE);

    // read file header
    (void) jpeg_read_header(&input_info, TRUE);
    
    jvirt_barray_ptr *coef_arrays;
    
    // allocate memory for reading out DCT coeffs
    for (JDIMENSION compnum = 0; compnum < input_info.num_components; ++compnum) {
        coef_buffers[compnum] = ((&input_info)->mem->alloc_barray)
            ((j_common_ptr) &input_info, JPOOL_IMAGE,
             input_info.comp_info[compnum].width_in_blocks,
             input_info.comp_info[compnum].height_in_blocks);
    }

    // read input file as DCT coeffs 
    coef_arrays = jpeg_read_coefficients(&input_info);
    
    // copy compression parameters from the input file to the output file
    jpeg_copy_critical_parameters(&input_info, &output_info);
    
    // copy DCT coeffs to a new array
    int num_components = input_info.num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        height_in_blocks[compnum] = input_info.comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = input_info.comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&input_info)->mem->access_virt_barray)
            ((j_common_ptr) &input_info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; ++blocknum)
            {
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef_buffers[compnum][rownum][blocknum][i] = row_ptrs[compnum][0][blocknum][i];
                }
            }
        }
    }
    
    // -------------------------------------------------------------------------------------------------
    // begin LSB methodology
    // loop over the coefficients and replace LSBs with data bits
    // -------------------------------------------------------------------------------------------------
    int number_of_coefficients = 0;
    int data_byte_index = 0;
    int data_bit_index = 0;
    int coef = 0;
    int coef_bits = 0;
    int data_bits = 0;
    int data_mask = (int)(pow((float)2, bits_to_steal) - 1);
    int coefficient_mask = (255 >> bits_to_steal) << bits_to_steal;
    int length_payload_bit_count = 0;
    int version_payload_bit_count = 0;
    int bits_to_steal_payload_count = 0;
    int length_payload = data_length ^ MAGIC_XOR;
    int data_length_minus_one = data_length - 1;
    int coef_max = COEFFICIENT_MAX - 256;
    int coef_min = -coef_max;
    unsigned short short_buf;
    unsigned short version = 0;
    int header_complete = 0;
    int encode_bit_count = 0;
    int stop = 0;
    time_t t;
    srand((unsigned) time(&t));
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        if (stop) {
            break;
        }
        
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            if (stop) {
                break;
            }
                
            for (JDIMENSION blocknum = 0; blocknum<width_in_blocks[compnum]; ++blocknum)
            {
                if (stop) {
                    break;
                }

                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef = coef_buffers[compnum][rownum][blocknum][i];
                    //only steal coefficients that don't have a zero value
                    if (coef != 0) {
                        if (header_complete) {
                            //if we have not encoded the entire payload
                            if (data_byte_index < data_length) {
                                // use a short buf here to pack the current
                                // byte and the next byte of the data.  this
                                // makes it easier for values for bits to steal like
                                // 3 which can spill over between bytes
                                short_buf = data[data_byte_index];
                                if (data_byte_index < data_length_minus_one) {
                                    short_buf = (short_buf << 8) | data[data_byte_index + 1];
                                } else {
                                    short_buf = short_buf << 8;;
                                }
                                
                                data_bits = (short_buf >> (16 - (data_bit_index + bits_to_steal))) & data_mask;
                                coef = (coef >> bits_to_steal << bits_to_steal) | data_bits;
                                if (coef == 0) {
                                    coef = coef | (1 << bits_to_steal);
                                }
                                
                                //add the bits we stole to the bit count and byte count
                                data_bit_index += bits_to_steal;
                                if (data_bit_index >= 8) {
                                    data_bit_index %= 8;
                                    ++data_byte_index;
                                }

                                coef_buffers[compnum][rownum][blocknum][i] = coef;
                            } else {
                                stop = 1;
                                break;
                            } 
                        } else {
                            // the header contains the version number, payload length,
                            // and bits to steal encoded by stealing 1 bit
                            // the payload will use the number of bits to steal passed
                            // into the function
                            // noise bit every 5 values
                            if (encode_bit_count % NOISE_MODULO_VALUE == 0) {
                                coef_bits = coef & 1;
                                data_bits = (rand() % 2) & 1;
                                    
                                coef_buffers[compnum][rownum][blocknum][i] = get_lossless_encode_bit(coef, coef_bits, data_bits);
                            } else if (version_payload_bit_count < VERSION_BIT_LENGTH) {
                                coef_bits = coef & 1;
                                data_bits = (CSTEG_JPEG_VERSION >> version_payload_bit_count) & 1;
                                
                                coef_buffers[compnum][rownum][blocknum][i] = get_lossless_encode_bit(coef, coef_bits, data_bits);
                                ++version_payload_bit_count;
                            } else if (length_payload_bit_count < DATA_LENGTH_BIT_LENGTH) {
                                coef_bits = coef & 1;
                                data_bits = (length_payload >> length_payload_bit_count) & 1;
                                
                                coef_buffers[compnum][rownum][blocknum][i] = get_lossless_encode_bit(coef, coef_bits, data_bits);
                                ++length_payload_bit_count;
                            } else if (bits_to_steal_payload_count < BITS_TO_STEAL_BIT_LENGTH) {
                                //encode the number of bits we stole when encoding this image
                                coef_bits = coef & 1;
                                data_bits = (bits_to_steal >> bits_to_steal_payload_count) & 1;
                                
                                coef_buffers[compnum][rownum][blocknum][i] = get_lossless_encode_bit(coef, coef_bits, data_bits); 
                                ++bits_to_steal_payload_count;
                            } else {
                                header_complete = 1;
                            }
                            
                            ++encode_bit_count;
                        } 

                        ++number_of_coefficients;
                    }
                }
            }
        }
    }
    
    // -----------------------------------------------------------------------------------------------------------------------------

    // output the new DCT coeffs to a JPEG file 
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&output_info)->mem->access_virt_barray)
                ((j_common_ptr) &output_info, coef_arrays[compnum],
                 rownum, (JDIMENSION) 1, TRUE);
            memcpy(row_ptrs[compnum][0][0],
                coef_buffers[compnum][rownum][0],
                 block_row_size[compnum]);
        }
    }
    
    // write to the output file 
    jpeg_write_coefficients(&output_info, coef_arrays);
    
    // finish compression and release memory 
    jpeg_finish_compress(&output_info);
    jpeg_destroy_compress(&output_info);
    jpeg_finish_decompress(&input_info);
    jpeg_destroy_decompress(&input_info);
    
    // close files
    fclose(input_file);
    fclose(output_file);
    
    return 1;
}

/**
 * this functino decodeds the LSB payload from a jpeg image that was embedded through the encode lossless
 * function.  it is assumed the image being decoded is exactly the same (byte for byte) as the output 
 * image from the encode lossless function.  if the encoded image was cropped, recompressed, or manipulated
 * in any way, the lossy encode/decode functions should be used.
 *
 * @param input_file_name
 *     the full file path/name to the image to be decoded
 * @param data_out
 *     a double unsigned char pointer.  this function will figure out what length array the supplied pointer should
 *     point to an reassign the point to point to it.  THIS FUNCTION WILL ALLOCATE MEMORY THAT IS POINTED TO BY
 *     THIS.  the memory allocated by this function should be freed elsewhere!!!!!! resposibly deallocate please!!!
 * @param data_out_length
 *     this unsigned int will be assigned a value equal to the length of the decoded data_out buffer
 *
 * returns true if everything went as expected, false otherwise.  memory for data_out will only be allocated if
 *     this function returned true.
 */
int decode_lossless_jpeg_file(const char *input_file_name, unsigned char** data_out, unsigned int* data_out_length) {
    // these are standard libjpeg structures for reading(decompression)
    struct jpeg_decompress_struct input_info;
    struct jpeg_error_mgr jerr;
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    // libjpeg data structure for storing one row, that is, scanline of an image
    
    // open the input and output files
    FILE* input_file = fopen(input_file_name, "rb");
    if (!input_file) {
        printf("error opening jpeg file %s\n!", input_file_name);
        return 0;
    }
    
    // initialize the JPEG compression and decompression objects with default error handling.
    input_info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&input_info);
    
    // specify data source for decompression and recompression
    jpeg_stdio_src(&input_info, input_file);
    
    // read file header
    (void) jpeg_read_header(&input_info, TRUE);
    
    // read input file as DCT coeffs 
    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&input_info);
    
    // -------------------------------------------------------------------------------------------
    // begin LSB decoding
    // -------------------------------------------------------------------------------------------
    int num_components = input_info.num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    int data_byte_index = 0;
    int data_bit_index = 0;
    int coef;
    int buffer_length = 0;
    unsigned char* buffer = NULL;
    unsigned int decode_bit_count = 0;
    unsigned char coef_bit = 0;
    int coefficient_mask = 0;
    int buff_buff = 1;
    int buff_buff_cnt = 0;
    int data_bits = 0;
    int stop = 0;
    int bits_to_steal = 0;
    int header_complete = 0;
    int version_payload_bit_count = 0;
    int bits_to_steal_payload_count = 0;
    int length_payload_bit_count = 0;
    int version = 0;
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        if (stop) {
            break;
        }
        height_in_blocks[compnum] = input_info.comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = input_info.comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            if (stop) {
                break;
            }
            row_ptrs[compnum] = ((&input_info)->mem->access_virt_barray)
            ((j_common_ptr) &input_info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; ++blocknum)
            {
                if (stop) {
                    break;
                }
                
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef = row_ptrs[compnum][0][blocknum][i];

                    if (coef != 0) {
                        if (header_complete) {
                            if (data_byte_index < buffer_length) {
                                data_bits = coef & coefficient_mask;
                                buff_buff = (buff_buff << bits_to_steal) | data_bits;
                                buff_buff_cnt += bits_to_steal;
                                if (buff_buff_cnt >= 8) {
                                    buff_buff_cnt -= 8;
                                    buffer[data_byte_index] = ((buff_buff >> buff_buff_cnt) & 255) ^ (data_byte_index % 2 == 1 ? MAGIC_XOR_ODD : MAGIC_XOR_EVEN);
                                    ++data_byte_index;
                                }
                            } else {
                                stop = true;
                                break;
                            }
                        } else {
                            if (decode_bit_count % NOISE_MODULO_VALUE == 0) {
                            } else if (version_payload_bit_count < VERSION_BIT_LENGTH) {
                                coef_bit = coef & 1;
                                
                                version = version | (coef_bit << version_payload_bit_count);
                                ++version_payload_bit_count;
                            } else if (length_payload_bit_count < DATA_LENGTH_BIT_LENGTH) {
                                coef_bit = coef & 1;
                                
                                buffer_length = buffer_length | (coef_bit << length_payload_bit_count);
                                ++length_payload_bit_count;
                            } else if (bits_to_steal_payload_count < BITS_TO_STEAL_BIT_LENGTH) {
                                //encode the number of bits we stole when encoding this image
                                coef_bit = coef & 1;

                                bits_to_steal = bits_to_steal | (coef_bit << bits_to_steal_payload_count);                            
                                ++bits_to_steal_payload_count;
                            } else { 
                                header_complete = 1;

                                buffer_length = buffer_length ^ MAGIC_XOR;
                                buffer = (unsigned char*) malloc(buffer_length);
                                memset(buffer, 0, buffer_length);
                                coefficient_mask = (int)pow((float)2, bits_to_steal) - 1;
                                
                                //debug info
                                //printf("length: %d\n", buffer_length);
                                //printf("mask: %d\n", coefficient_mask);
                                //printf("bits: %d\n", bits_to_steal);
                                //printf("version: %d\n", version);
                            }
                            
                            ++decode_bit_count;
                        }
                    }
                }
            }
        }
    }
    
    // finish compression and release memory 
    jpeg_finish_decompress(&input_info);
    jpeg_destroy_decompress(&input_info);
    
    // close file
    fclose(input_file);

    // set the supplied pointer to point to this memory
    // DEALLOCATING THIS MEMORY IS DONEN BY THE CALLER
    // yeah i know that is dangerous
    *data_out = buffer;
    *data_out_length = buffer_length;

    return 1;
}

/**
 * this file returns the number of usable coefficients a lossless encode can expect when 
 * encoding information.  this can be used to check whether payloads will fit into the image
 *
 * @param input_file_name
 *     the path to a jpeg image file to use
 *
 * returns the number of non-zero coefficients data can be embed into.  this number can 
 *     be used to see whether data will fit into the image
 */
unsigned int lossless_jpeg_file_get_number_of_coefficients(const char *input_file_name) {
    // these are standard libjpeg structures for reading(decompression)
    struct jpeg_decompress_struct input_info;
    struct jpeg_error_mgr jerr;
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    // libjpeg data structure for storing one row, that is, scanline of an image
    
    // open the input and output files
    FILE* input_file = fopen(input_file_name, "rb");
    if (!input_file) {
        printf("error opening jpeg file %s\n!", input_file_name);
        return 0;
    }
    
    // initialize the JPEG compression and decompression objects with default error handling.
    input_info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&input_info);
    
    // specify data source for decompression and recompression
    jpeg_stdio_src(&input_info, input_file);
    
    // read file header
    (void) jpeg_read_header(&input_info, TRUE);

    // read input file as DCT coeffs 
    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&input_info);
    
    // copy DCT coeffs to a new array
    int num_components = input_info.num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    int number_of_coefficients = 0;
    short coef = 0;

    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        height_in_blocks[compnum] = input_info.comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = input_info.comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&input_info)->mem->access_virt_barray)
            ((j_common_ptr) &input_info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; ++blocknum)
            {
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef = row_ptrs[compnum][0][blocknum][i];
                    
                    // we only use non-zero coefficients
                    if (coef != 0) {
                       ++ number_of_coefficients;
                    }
                }
            }
        }
    }
    
    // finish compression and release memory 
    jpeg_finish_decompress(&input_info);
    jpeg_destroy_decompress(&input_info);
    
    // close file
    fclose(input_file);

    return number_of_coefficients - CSTEG_JPEG_VERSION_1_HEADER_LENGTH;
}

/*
 * this function encodes shorter amounts of information redundantly into a jpeg image.  this information
 * can survive manipulation such as recompression or cropping.  this function take a specially formated input
 * that is stored into a packet_collection.  this function puts packeted into the single LSB of non-zero 
 * coefficients.  those packets have information beyond just the original data (checksum, rolling sum, packet
 * index, etc.) that allows the decode lossy function to attempt a decode.  the encode and decode methods here
 * are not responsible to formatting the packets of reassembling the original data.  that is done elsewhere
 * this method only encodes supplied data 1 bit at a time into "more reliable/stable" DCT coefficients.  there
 * is a layer above this that does the packet formatting and unpacking.
 *
 * @param input_file_name
 *     the full path/name to a jpeg image on disk to encode
 * @param output_file_name
 *     the full path/name to store the encoded jpeg image on disk
 * @param packet_collection
 *     the packet_collection object that holds the formatted packets to encode into the single LSB of
 *     the non-zero coefficients of the input image
 * @param quality
 *     the quality to save the output image with (it would be advisable to make this equal to the 
 *     quality of the input image)
 *
 * returns true if everything went as planned, false if an error occurred
 */
int encode_lossy_jpeg_file(const char *input_file_name, const char* output_file_name,
                            struct packet_collection* packet_collection, int quality) {
    // these are standard libjpeg structures for reading(decompression)
    struct jpeg_decompress_struct input_info;
    struct jpeg_compress_struct output_info;
    struct jpeg_error_mgr jerr;
    JBLOCKARRAY coef_buffers[MAX_COMPONENTS];
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    // libjpeg data structure for storing one row, that is, scanline of an image
    
    // open the input and output files
    FILE* input_file = fopen(input_file_name, "rb");
    if ( !input_file ) {
        printf("error opening jpeg file %s\n!", input_file_name);
        return 0;
    }
    
    FILE* output_file = fopen(output_file_name, "wb");
    if (!output_file) {
        printf("error opening jpeg file %s\n!", output_file_name);
        return 0;
    }

    // initialize the JPEG compression and decompression objects with default error handling.
    input_info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&input_info);
    output_info.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&output_info);
    
    // specify data source for decompression and recompression
    jpeg_stdio_src(&input_info, input_file);
    jpeg_stdio_dest(&output_info, output_file);

    // default compression parameters, we shouldn't be worried about these    
    //jpeg_set_defaults(&output_info);

    //output_info.num_components = 3;
    //info.data_precision = 4;
    output_info.dct_method = JDCT_ISLOW;
    //info.dct_method = JDCT_IFAST;
    jpeg_set_quality(&output_info, quality, TRUE);

    
    // read file header
    (void) jpeg_read_header(&input_info, TRUE);
    
    jvirt_barray_ptr *coef_arrays;
    
    // allocate memory for reading out DCT coeffs
    for (JDIMENSION compnum = 0; compnum < input_info.num_components; ++compnum) {
        coef_buffers[compnum] = ((&input_info)->mem->alloc_barray)
            ((j_common_ptr) &input_info, JPOOL_IMAGE,
             input_info.comp_info[compnum].width_in_blocks,
             input_info.comp_info[compnum].height_in_blocks);
    }

    // read input file as DCT coeffs 
    coef_arrays = jpeg_read_coefficients(&input_info);
    
    // copy compression parameters from the input file to the output file
    jpeg_copy_critical_parameters(&input_info, &output_info);
    
    // copy DCT coeffs to a new array
    int num_components = input_info.num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        height_in_blocks[compnum] = input_info.comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = input_info.comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&input_info)->mem->access_virt_barray)
            ((j_common_ptr) &input_info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; ++blocknum)
            {
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef_buffers[compnum][rownum][blocknum][i] = row_ptrs[compnum][0][blocknum][i];
                }
            }
        }
    }
    
    int number_of_coefficients = 0;
    int data_byte_index = 0;
    int data_bit_index = 0;
    int coef = 0;
    int coef_bit = 0;
    int data_bit = 0;
    int x = 0;
    int y = 0;
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            for (JDIMENSION blocknum = 0; blocknum<width_in_blocks[compnum]; ++blocknum)
            {
                x = 0;
                y = 0;
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef = coef_buffers[compnum][rownum][blocknum][i];
                    x = (x + 1) % DCTSIZE;
                    if (x == 0) {
                        ++y;
                    }

                    if (coef != 0 && x > 1 && y > 1) {
                        coef_bit = coef & 1;
                        data_bit = (packet_collection->data[data_byte_index] >> data_bit_index) & 1;

                        ++data_bit_index;
                        data_bit_index %= 8;
                        if (data_bit_index == 0) {
                            ++data_byte_index;
                            data_byte_index %= packet_collection->length;
                        }
                        
                        if (data_bit != coef_bit) {
                            switch(coef) {
                            case 1:
                                ++coef;
                                break;
                            case -1:
                                --coef;
                                break;
                            default:
                                if (coef >= COEFFICIENT_MAX) {
                                    --coef;
                                } else if (coef <= COEFFICIENT_MIN) {
                                    ++coef;
                                } else {
                                    if ((coef ^ number_of_coefficients) == 1) {
                                        ++coef;
                                    } else {
                                        --coef;
                                    }
                                }
                            }
                        }

                        coef_buffers[compnum][rownum][blocknum][i] = coef;
                        ++number_of_coefficients;
                    }
                }
            }
        }
    }
    
    // output the new DCT coeffs to a JPEG file 
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&output_info)->mem->access_virt_barray)
                ((j_common_ptr) &output_info, coef_arrays[compnum],
                 rownum, (JDIMENSION) 1, TRUE);
            memcpy(row_ptrs[compnum][0][0],
                coef_buffers[compnum][rownum][0],
                 block_row_size[compnum]);
        }
    }
    
    // write to the output file 
    jpeg_write_coefficients(&output_info, coef_arrays);
    
    // finish compression and release memory 
    jpeg_finish_compress(&output_info);
    jpeg_destroy_compress(&output_info);
    jpeg_finish_decompress(&input_info);
    jpeg_destroy_decompress(&input_info);
    
    // close files
    fclose(input_file);
    fclose(output_file);
    
    return 1;
}

/**
 * this function attempts to decode information encoded into a jpeg image using the lossy data encoding function.
 * this function puts the extracted bit stream into a packet collection. it should be assumed by layers above 
 * this one that the packet stream extracted from this image is corrupted.  bits could be flipped or dropped, 
 * which means layers above this should not even assume the data is byte aligned!  the encode and decode methods 
 * here are not responsible to formatting the packets of reassembling the original data.  that is done elsewhere
 * this method only encodes supplied data 1 bit at a time into "more reliable/stable" DCT coefficients.  there
 * is a layer above this that does the packet formatting and unpacking.
 *
 * @param input_file_name
 *     the full path/name to the input jpeg file to decode
 * @param packet_collection
 *     the packet collection object to store the extracted bit stream from the single LSB of non-zero DCT
 *     coefficients
 *
 * returns true if the information was extracted, false only if there was an error like file open error
 */
int decode_lossy_jpeg_file(const char *input_file_name, struct packet_collection* packet_collection) {
    // these are standard libjpeg structures for reading(decompression)
    struct jpeg_decompress_struct input_info;
    struct jpeg_error_mgr jerr;
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    // libjpeg data structure for storing one row, that is, scanline of an image
    
    // open the input and output files
    FILE* input_file = fopen(input_file_name, "rb");
    if (!input_file) {
        printf("error opening jpeg file %s\n!", input_file_name);
        return 0;
    }
    
    // initialize the JPEG compression and decompression objects with default error handling.
    input_info.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&input_info);
    
    // specify data source for decompression and recompression
    jpeg_stdio_src(&input_info, input_file);
    
    // read file header
    (void) jpeg_read_header(&input_info, TRUE);
    
    // read input file as DCT coeffs 
    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&input_info);

    // ----------------------------------------------------------------------------
    // begin LSB decoding 1 bit only for lossy encoding/decoding
    // ----------------------------------------------------------------------------
    int num_components = input_info.num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    int number_of_coefficients = 0;
    int data_byte_index = 0;
    int data_bit_index = 0;
    int coef;
    int buffer_length = 
        (sizeof(unsigned char) * 
        input_info.image_width *
        input_info.image_height *
         input_info.num_components) / 8;

    unsigned char* buffer = (unsigned char*) malloc(buffer_length);
    memset(buffer, 0, buffer_length);
    
    int x = 0;
    int y = 0;
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        height_in_blocks[compnum] = input_info.comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = input_info.comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            row_ptrs[compnum] = ((&input_info)->mem->access_virt_barray)
            ((j_common_ptr) &input_info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; ++blocknum)
            {
                x = 0;
                y = 0;
                for (JDIMENSION i = 0; i < DCTSIZE2; ++i)
                {
                    coef = row_ptrs[compnum][0][blocknum][i];
                    x = (x + 1) % DCTSIZE;
                    if (x == 0) {
                        ++y;
                    }

                    if (coef != 0 && x > 1 && y > 1) {
                        if ((coef & 1) == 1) {
                            buffer[data_byte_index] = buffer[data_byte_index] | (1 << data_bit_index);
                        }

                        ++data_bit_index;
                        data_bit_index %= 8;
                        if (data_bit_index == 0) {
                            ++data_byte_index;
                        }

                        ++number_of_coefficients;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------
    
    // finish compression and release memory 
    jpeg_finish_decompress(&input_info);
    jpeg_destroy_decompress(&input_info);
    
    // close file
    fclose(input_file);

    packet_collection->data = buffer;
    packet_collection->length = number_of_coefficients / 8;

    return 1;
}

/*
 * this function loops through two coefficient store structs and returns the number of coefficients
 * between the two that differ.  this is a debug/tool function that probably wont be used in production
 *
 * @param coefficient_store_1
 *     the coefficients from the 1st image in the diff
 * @param coefficient_store_2
 *     the coefficients from the 2nd image in the diff
 *
 * returns the number of coefficients that differ between the two sets of coefficients
 */
int diff_coefficients(struct coefficient_store* coefficient_store_1, struct coefficient_store* coefficient_store_2) {
    int diff = 0;
    int count = 0;
    int num_components = coefficient_store_1->info->num_components;
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    int x = 0;
    int y = 0;
    for (JDIMENSION compnum = 0; compnum < coefficient_store_1->info->num_components; ++compnum)
    {
        height_in_blocks[compnum] = coefficient_store_1->info->comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = coefficient_store_1->info->comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof(JCOEF) * DCTSIZE2 * width_in_blocks[compnum];
        for (JDIMENSION rownum=0; rownum < height_in_blocks[compnum]; ++rownum)
        {
            for (JDIMENSION blocknum=0; blocknum<width_in_blocks[compnum]; blocknum++)
            {
                x = 0;
                y = 0;
                for (JDIMENSION i = 0; i < DCTSIZE2; i++)
                {
                    x = (x + 1) % DCTSIZE;
                    if (x == 0) {
                        ++y;
                    }

                    int val1 = coefficient_store_1->coefficient_buffer[compnum][rownum][blocknum][i];
                    int val2 = coefficient_store_2->coefficient_buffer[compnum][rownum][blocknum][i];

                    if (val1 != val2) {
                        printf("%d != %d, for %d, %d at rowblock %d, %d\n", val1, val2, x, y, rownum, blocknum);
                        ++diff;
                    }
                }
            }
        }
    }

    printf("number of diffs = %d\n", diff);

    return diff;
}

/**
 * this function loads the coefficients from a jpeg image into a coefficient store struct
 *
 * @param jpeg_file_name
 *     the full path/name of the jpeg file on disk for which to load coefficients from
 *
 * returns the coefficients in the jpeg file inside a coefficient_store struct
 */
struct coefficient_store* load_coefficients(const char *jpeg_file_name)
{
    struct jpeg_error_mgr jerr;
    JDIMENSION i, rownum, blocknum;
    JBLOCKARRAY* coef_buffers = (JBLOCKARRAY*)malloc(sizeof(JBLOCKARRAY) * MAX_COMPONENTS);
    JBLOCKARRAY row_ptrs[MAX_COMPONENTS];
    
    /* Open the input and output files */
    FILE* jpeg_file = fopen(jpeg_file_name, "rb");
    if (!jpeg_file) {
        printf("error opening jpeg file %s\n!", jpeg_file_name);
        return NULL;
    }
    
    // initialize the JPEG compression and decompression objects with default error handling
    struct jpeg_decompress_struct* info = (struct jpeg_decompress_struct*) malloc(sizeof(struct jpeg_decompress_struct));
    (*info).err = jpeg_std_error(&jerr);
    jpeg_create_decompress(info);
    
    // specify data source for decompression 
    jpeg_stdio_src(info, jpeg_file);
    
    // read file header
    (void) jpeg_read_header(info, TRUE);
    
    // allocate memory for reading out DCT coeffs 
    for (JDIMENSION compnum = 0; compnum < (*info).num_components; ++compnum) {
        coef_buffers[compnum] =
            (info->mem->alloc_barray)
                ((j_common_ptr) info, JPOOL_IMAGE,
                 (*info).comp_info[compnum].width_in_blocks,
                 (*info).comp_info[compnum].height_in_blocks);
    }
    
    // read input file as DCT coeffs 
    jvirt_barray_ptr *coef_arrays = jpeg_read_coefficients(info);
    
    // copy DCT coeffs to a new array 
    int num_components = (*info).num_components;
    
    size_t block_row_size[num_components];
    int width_in_blocks[num_components];
    int height_in_blocks[num_components];
    for (JDIMENSION compnum = 0; compnum < num_components; ++compnum)
    {
        height_in_blocks[compnum] = (*info).comp_info[compnum].height_in_blocks;
        width_in_blocks[compnum] = (*info).comp_info[compnum].width_in_blocks;
        block_row_size[compnum] = (size_t) sizeof( JCOEF ) * DCTSIZE2 * width_in_blocks[compnum];
        for (rownum=0; rownum<height_in_blocks[compnum]; rownum++)
        {
            row_ptrs[compnum] = (info->mem->access_virt_barray)
            ((j_common_ptr) info, coef_arrays[compnum],
             rownum, (JDIMENSION) 1, FALSE);
            for (blocknum=0; blocknum<width_in_blocks[compnum]; blocknum++)
            {
                for (i=0; i<DCTSIZE2; i++)
                {
                    coef_buffers[compnum][rownum][blocknum][i] = row_ptrs[compnum][0][blocknum][i];
                }
            }
        }
    }
    
    fclose(jpeg_file);
    
    struct coefficient_store* response = 
        (struct coefficient_store*)malloc(sizeof(struct coefficient_store));
    response->coefficient_buffer = coef_buffers;
    response->info = info;
    
    return response;
}

/**
 * this function takes an image and runs a series of compressions and then recompressions n times.
 * the series of images is left on disk, use purge to clean up the images.  uncompressing and then 
 * recompressing an image makes the coefficients more stable for subsequent recompressions.  this
 * means that when facebook recompresses an image it will destroy less of our information.  the 
 * images on disk will be names using the input image name.  so /some/path/my.jpeg will be turned
 * into a sequence of /some/path/my-0.jpeg, /some/path/my-1.jpeg, etc.  the final image in the 
 * series will be n - 1.
 *
 * @param n
 *     the number of times to perform uncompression/recompression
 * @param file_path
 *     the full path/name to the input image on disk
 * @param quality
 *     the quality to use when compressing the image
 */
void do_compressions(int n, const char* file_path, int quality) {
    char file_path_buffer[256];
    get_without_extension(file_path, file_path_buffer);

    char file_name_buffer_1[256];
    char file_name_buffer_2[256];
    memset(file_name_buffer_1, 0, 256);
    
    char* ptr_1 = file_name_buffer_1;
    char* ptr_2 = file_name_buffer_2;
        
    strcat(ptr_1, file_path);
    struct image image;
    construct_image(&image);
    
    for (int i = 0; i < n; ++i) {
        memset(ptr_2, 0, 256);
        sprintf(ptr_2, "%s-%d.jpeg", file_path_buffer, i);

        read_jpeg_file(&image, ptr_1);
        write_jpeg_file(&image, ptr_2, quality);
        
        char* temp = ptr_1;
        ptr_1 = ptr_2;
        ptr_2 = temp;
    }

    free(image.raw_image);
}

/**
 * this function does a coefficient diff between to jpeg image files on disk
 *
 * @param file_path_1
 *     the full path/name to the first image on disk
 * @param file_path_2
 *     the full path/name to the second image on disk
 */
int do_coefficient_diff(const char* file_path_1, const char* file_path_2) {
    struct coefficient_store* coefficient_store_1 = load_coefficients(file_path_1);
    struct coefficient_store* coefficient_store_2 = load_coefficients(file_path_2);

    int rtn = diff_coefficients(coefficient_store_1, coefficient_store_2);

    destruct_coefficient_store(coefficient_store_1);
    destruct_coefficient_store(coefficient_store_2);

    free(coefficient_store_1);
    free(coefficient_store_2);
    
    return rtn;
}

/**
 * this funtion runs over a compression series and diffs n and n - 1.
 *
 * @param n
 *     the n used for the compression series (or a lesser n)
 * @param file_path
 *     the name of the original image input into the compression series
 */
void do_compression_series_coefficient_diffs(int n, const char* file_path) {
    char file_path_buffer[256];
    get_without_extension(file_path, file_path_buffer);

    char buffer_1[256];
    char buffer_2[256];
    memset(buffer_1, 0, 256);
    
    char* ptr_1 = buffer_1;
    char* ptr_2 = buffer_2;
    
    strcat(ptr_1, file_path);
    for (int i = 0; i < n; ++i) {
        memset(ptr_2, 0, 256);
        sprintf(ptr_2, "%s-%d.jpeg", file_path_buffer, i);
        
        if (i > 0) {
            struct coefficient_store* coefficient_store_1 = load_coefficients((const char*)ptr_1);
            struct coefficient_store* coefficient_store_2 = load_coefficients((const char*)ptr_2);
            
            printf("\n\n------------------------------\n");
            printf("%s\n%s\n", ptr_1, ptr_2);
            diff_coefficients(coefficient_store_1, coefficient_store_2);

            destruct_coefficient_store(coefficient_store_1);
            destruct_coefficient_store(coefficient_store_2);

            free(coefficient_store_1);
            free(coefficient_store_2);
        }

        char* temp = ptr_1;
        ptr_1 = ptr_2;
        ptr_2 = temp;
    }
}

/**
 * this function deletes the intermediate files in a compression series
 *
 * @param jpeg_file_path
 *     the original input path/name for the compression series (this file 
 *     will not be deleted)
 * @param limit
 *     every file in the compression series (except the original iamge) 
 *     whose n value is strictly less than limit will be deleted
 */
void purge(const char* jpeg_file_path, int limit) {
    char jpeg_path_buffer[256];
    get_without_extension(jpeg_file_path, jpeg_path_buffer);

    char file_path_buffer[256];
    for (int i = 0; i < limit; ++i) {
        memset(file_path_buffer, 0, 256);
        sprintf(file_path_buffer, "%s-%d.jpeg", jpeg_path_buffer, i);
        remove(file_path_buffer);
    }
}

/**
 * this function deletes the intermediate files in a compression series
 *
 * @param directory
 *     the directory the original jpeg resided in
 * @param jpeg_file_name
 *     the original input file name for the compression series (this file 
 *     will not be deleted)
 * @param limit
 *     every file in the compression series (except the original iamge) 
 *     whose n value is strictly less than limit will be deleted
 */
void purge(const char* directory, const char* jpeg_file_name, int limit) {
    char jpeg_name_buffer[256];
    get_without_extension(jpeg_file_name, jpeg_name_buffer);

    char file_path_buffer[256];
    for (int i = 0; i < limit; ++i) {
        memset(file_path_buffer, 0, 256);
        sprintf(file_path_buffer, "%s/%s-%d.jpeg", directory, jpeg_name_buffer, i);
        remove(file_path_buffer);
    }
}

/**
 * this function moves the final image in a compression series to a file named -n-final.jpeg
 *
 * @param limit
 *     the n used in the compression series
 * @param quality
 *     the quality number used in the compression series
 * @param directory
 *      the directory the used for the compressions series
 * @param jpeg_file_name
 *     the name of the original jpeg file used in the compression series
 */
void move_quality_result(int limit, int quality, const char* directory, const char* jpeg_file_name) {
    char jpeg_name_buffer[256];
    get_without_extension(jpeg_file_name, jpeg_name_buffer);

    char file_name_buffer_1[256];
    char file_name_buffer_2[256];
    memset(file_name_buffer_1, 0, 256);
    memset(file_name_buffer_2, 0, 256);
    
    --limit;
    sprintf(file_name_buffer_1, "%s/%s-%d.jpeg", directory, jpeg_name_buffer, limit);
    sprintf(file_name_buffer_2, "%s/%s-%d-final.jpeg", directory, jpeg_name_buffer, quality);
    
    rename(file_name_buffer_1, file_name_buffer_2);
}

/**
 * this function does a series of compressions series and stores the result of each
 * the final image in the series is renamed to my-jpeg-[qualitynumber]-final.jpeg
 * 
 * @param bottom
 *     the lowest quality number to use in the series of compression series
 * @param top
 *      the highest quality number to use in the series of compression series
 * @param number_of_compressions
 *     the number of compressions to run for each quality number
 * @param directory
 *     the diretory the original input image resides in
 * @param jpeg_file_name
 *     the name of the original jpeg image
 */
void make_quality_panel(int bottom, int top, int number_of_compressions, const char* directory, const char* jpeg_file_name) {
    char path_to_jpeg_file[256];
    memset(path_to_jpeg_file, 0, 256);
    strcat(path_to_jpeg_file, directory);
    strcat(path_to_jpeg_file, "/");
    strcat(path_to_jpeg_file, jpeg_file_name);

    for (int i = bottom; i <= top; ++i) {
        printf("generating image for quality: %d\n", i);
        do_compressions(number_of_compressions, path_to_jpeg_file, i);
        purge(directory, jpeg_file_name, number_of_compressions - 1);
        move_quality_result(number_of_compressions, i, directory, jpeg_file_name);
    }
}
