//
//  main.cpp
//  c-steg
//

#include "CStegApi.h"
#include <string.h>
#include <iostream>

void print_help() {
    std::cout << "this is the c-steg command line utility" << std::endl << std::endl;

    std::cout << "\tto encode a lossless jpeg image with text" << std::endl;
    std::cout << "\t\t--encode-lossless" << std::endl;
    std::cout << "\t\t--input-file\t\t<required string (path to input image)>" << std::endl;
    std::cout << "\t\t--output-file\t\t<required string (path to output image)>" << std::endl;
    std::cout << "\t\t--text\t\t\t<required string>" << std::endl;
    std::cout << "\t\t--quality\t\t<optional int defaults to 71>" << std::endl;
    std::cout << "\t\t--bits-to-steal\t\t<optional int defaults to 2>" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "\tto decode a lossless jpeg" << std::endl;
    std::cout << "\t\t--decode-lossless" << std::endl;
    std::cout << "\t\t--input-file\t\t<required string (path to input image)>" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "\tto encode a lossy jpeg image with text" << std::endl;
    std::cout << "\t\t--encode-lossy" << std::endl;
    std::cout << "\t\t--input-file\t\t<required string (path to input image)>" << std::endl;
    std::cout << "\t\t--output-file\t\t<required string (path to output image)>" << std::endl;
    std::cout << "\t\t--text\t\t\t<required string>" << std::endl;
    std::cout << "\t\t--quality\t\t<optional int defaults to 71>" << std::endl;
    std::cout << "\t\t--compressions\t\t<optional int defaults to 5 (number of times " << std::endl
    		<< "\t\t\t\t\t to uncompress/recompress the image before " << std::endl
    		<< "\t\t\t\t\t encoding the payload)>" << std::endl;
    std::cout << "" << std::endl;
    
    std::cout << "\tto decode a lossy jpeg" << std::endl;
    std::cout << "\t\t--decode-lossy" << std::endl;
    std::cout << "\t\t--input-file\t\t<required string (path to input image)>" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "\tto get the number of usable coefficients in a file" << std::endl;
    std::cout << "\t\t--coefficients" << std::endl;
    std::cout << "\t\t--input-file\t\t<required string (path to input image)>" << std::endl;
    std::cout << "" << std::endl;

    std::cout << "\tto do a coefficient diff for two jpeg images" << std::endl;
    std::cout << "\t\t--diff" << std::endl;
    std::cout << "\t\t--input-file-1\t\t<required string (path to an image)>" << std::endl;
    std::cout << "\t\t--input-file-2\t\t<required string (path to an image)>" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
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
    int quality = 71;
    int number_of_compressions = 5;
    int bits_to_steal = 2;
    
    if (argc <= 1) {
    	std::cout << "***********************" << std::endl;
    	std::cout << "you need to pass in some command line flags, here is the help documentation" << std::endl;
    	std::cout << "***********************" << std::endl << std::endl;
    	print_help();
    }

    for (int i = 1; i < argc; ++i) {
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
        } else {
        	throw "invalid command line arg, use --help to see what commands are valid";
        }
    }

    if (encode_lossy) {
        unsigned int length = strlen(text);
        CSteg::encodeLossyJpeg(input_file, output_file, text, length,  quality, number_of_compressions);
    } else if (decode_lossy) {
        DecodeResponse decodeResponse;
        CSteg::decodeLossyJpeg(input_file, decodeResponse);
        std::cout << "data out length: " << decodeResponse.getLength() << std::endl;
        
        std::string printString;
        printString.assign(decodeResponse.getData(), decodeResponse.getData() + decodeResponse.getLength());
        std::cout << printString << std::endl;
    } else if (diff) {
        CSteg::doCoefficientDiff(input_file_1, input_file_2);
    } else if (encode_lossless) {
        unsigned int length = strlen(text);
        CSteg::encodeLosslessJpeg(input_file, output_file, (unsigned char*)text, length, bits_to_steal, quality);
        std::cout << std::endl;
    } else if (decode_lossless) {
        DecodeResponse decodeResponse;
        CSteg::decodeLosslessJpeg(input_file, decodeResponse);
        std::cout << "data out length: " << decodeResponse.getLength() << std::endl;

        std::string printString;
        printString.assign(decodeResponse.getData(), decodeResponse.getData() + decodeResponse.getLength());
        std::cout << printString << std::endl;
    } else if (number_of_coefficients) {
        std::cout << "number of coefficients in file: " << CSteg::getNumberOfCoefficients(input_file) << std::endl;
    }

    return 0;
}

