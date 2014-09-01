//
//  lossy_data_formatter.h
//  c-steg
//

#ifndef DATA_FORMATTER_H
#define DATA_FORMATTER_H

#define DATA_FORMATTER_PACKET_LENGTH 4
#define NUMBER_OF_NECESSARY_VALID_PACKETS 6
#define NUMBER_OF_PACKETS_FOR_BYTE_ALIGNMENT 2
#define NUMBER_OF_BITS_TO_REWIND 4
#define TERMINATE -1
#define BAD_PACKET 1
#define FIXED_IT 2
#define BYTE_ALIGNMENT_PROBLEM -2

struct packet {
    unsigned short data_length;
    unsigned char bytes[DATA_FORMATTER_PACKET_LENGTH];
    unsigned char index;
    unsigned char sum;
    unsigned char check_sum;
};

struct packet_collection {
    unsigned int length;
    unsigned char* data;
};

struct packet_collection_enumerator {
    int current_bit;
    int current_byte;
    unsigned int packet_collection_length;
};

unsigned short get_number_of_packets(unsigned short data_length);
struct packet_collection* to_lossy(const char* data, unsigned short data_length);
void from_lossy(struct packet_collection* packet_collection, unsigned char** data_out, unsigned int* data_out_length);

#endif
