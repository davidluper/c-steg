//
//  main.cpp
//  steg
//

#include "jpeg_tools.h"
#include "lossy_data_formatter.h"
#include "c_steg_api.h"
#include <string.h>

void print_help() {
    printf("this is the c-steg command line utility\n\n");

    printf("\tto encode a lossless jpeg image with text\n");
    printf("\t\t--encode-lossless\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\t\t--output-file\t\t<required string (path to output image)>\n");
    printf("\t\t--text\t\t\t<required string>\n");
    printf("\t\t--quality\t\t<optional int defaults to 71>\n");
    printf("\t\t--bits-to-steal\t\t<optional int defaults to 2>\n");
    printf("\n");

    printf("\tto decode a lossless jpeg\n");
    printf("\t\t--decode-lossless\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\n");

    printf("\tto encode a lossy jpeg image with text\n");
    printf("\t\t--encode-lossy\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\t\t--output-file\t\t<required string (path to output image)>\n");
    printf("\t\t--text\t\t\t<required string>\n");
    printf("\t\t--quality\t\t<optional int defaults to 71>\n");
    printf("\t\t--compressions\t\t<optional int defaults to 5 (number of times to uncompress/recompress the image)>\n");
    printf("\n");
    
    printf("\tto decode a lossy jpeg\n");
    printf("\t\t--decode-lossy\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\n");

    printf("\tto get the number of usable coefficients in a file\n");
    printf("\t\t--coefficients\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\n");

    printf("\tto do a coefficient diff for two jpeg images\n");
    printf("\t\t--diff\n");
    printf("\t\t--input-file-1\t\t<required string (path to an image)>\n");
    printf("\t\t--input-file-2\t\t<required string (path to an image)>\n");
    printf("\n");
    
    printf("\tto do a compression series/diff\n");
    printf("\t\t--diff-series\n");
    printf("\t\t--input-file\t\t<required string (path to input image)>\n");
    printf("\t\t--quality\t\t<optional int defaults to 71>\n");
    printf("\t\t--compressions\t\t<optional int defaults to 5 (number of times to uncompress/recompress the image)>\n");
    printf("\n");

    printf("\n");
}

int main(int argc, const char* argv[])
{
    const char* input_file = NULL;
    const char* output_file = NULL;
    const char* input_file_1 = NULL;
    const char* input_file_2 = NULL;
    const char* text = NULL;
    bool encode_lossy = false;
    bool encode_lossless = false;
    bool decode_lossless = false;
    bool decode_lossy = false;
    bool number_of_coefficients = false;
    bool diff = false;
    bool diff_series = false;
    int quality = 71;
    int number_of_compressions = 5;
    int bits_to_steal = 2;
    
    for (int i = 0; i < argc; ++i) {
        if (strcmp("--help", argv[i]) == 0) {
            print_help();
            return 0;
        } else if (strcmp("--input-file", argv[i]) == 0) {
            input_file = argv[++i];
        } else if (strcmp("--input-file-1", argv[i]) == 0) {
            input_file_1 = argv[++i];
        } else if (strcmp("--input-file-2", argv[i]) == 0) {
            input_file_2 = argv[++i];
        } else if (strcmp("--output-file", argv[i]) == 0) {
            output_file = argv[++i];
        } else if (strcmp("--text", argv[i]) == 0) {
            text = argv[++i];
        } else if (strcmp("--encode-lossy", argv[i]) == 0) {
            encode_lossy = true;
        } else if (strcmp("--coefficients", argv[i]) == 0) {
            number_of_coefficients = true;
        } else if (strcmp("--encode-lossless", argv[i]) == 0) {
            encode_lossless = true;
        } else if (strcmp("--decode-lossless", argv[i]) == 0) {
            decode_lossless = true;
        } else if (strcmp("--quality", argv[i]) == 0) {
            sscanf(argv[++i],"%d", &quality);
        } else if (strcmp("--bits-to-steal", argv[i]) == 0) {
            sscanf(argv[++i],"%d", &bits_to_steal);
        } else if (strcmp("--compressions", argv[i]) == 0) {
            sscanf(argv[++i],"%d", &number_of_compressions);
        } else if (strcmp("--decode-lossy", argv[i]) == 0) {
            decode_lossy = true;
        } else if (strcmp("--diff", argv[i]) == 0) {
            diff = true;
        } else if (strcmp("--diff-series", argv[i]) == 0) {
            diff_series = true;
        }
    }

    if (encode_lossy) {
        unsigned int length = strlen(text);
        c_steg_encode_lossy_jpeg(input_file, output_file, text, length,  quality, number_of_compressions);
    } else if (decode_lossy) {
        unsigned char* data_out = NULL;
        unsigned int data_out_length = 0;
        c_steg_decode_lossy_jpeg(input_file, &data_out, &data_out_length);
        
        printf("data_out_length: %d\n", data_out_length);        
        for (int i = 0; i < data_out_length; ++i) {
            printf("%c", data_out[i]);
        }

        printf("\n");

        free(data_out);
    } else if (diff) {
        do_coefficient_diff(input_file_1, input_file_2);
    } else if (diff_series) {
        do_compressions(number_of_compressions, input_file, quality);
        do_compression_series_coefficient_diffs(number_of_compressions, input_file);
        purge(input_file, number_of_compressions);
    } else if (encode_lossless) {
        unsigned int length = strlen(text);
        c_steg_encode_lossless_jpeg(input_file, output_file, (unsigned char*)text, length, bits_to_steal, quality);
        printf("\n");
    } else if (decode_lossless) {
        unsigned char* data_out = NULL;
        unsigned int data_out_length = 0;
        
        c_steg_decode_lossless_jpeg(input_file, &data_out, &data_out_length);

        printf("data_out_length: %d\n", data_out_length);
        for (int i = 0; i < data_out_length; ++i) {
            printf("%c", data_out[i]);
        }
        printf("\n");

        free(data_out);
    } else if (number_of_coefficients) {
        printf("number of coefficients in file: %d\n", c_steg_lossless_jpeg_file_get_number_of_coefficients(input_file));
    }

    return 0;
}

