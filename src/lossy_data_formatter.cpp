//
//  lossy_data_formatter.cpp
//  c-steg
//
//  this file contains functions that format and retrieve data into
//  packets that can survive bits flipping or dissapearing
//

#include "lossy_data_formatter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/**
 * this function returns the number of packets needed for transmitting
 * a data payload.
 *
 * @param data_length
 *     the length of the data being wrapped in packets
 */
unsigned short get_number_of_packets(unsigned short data_length) {
    int number_of_packets = (data_length / DATA_FORMATTER_PACKET_LENGTH);
    int number_of_packets_remainder = data_length % DATA_FORMATTER_PACKET_LENGTH;

    if (number_of_packets_remainder != 0) {
        ++number_of_packets;
    }

    return number_of_packets;
}

/**
 * this function wraps the input data into packets and puts them in a
 * packet collection.  the packet collection is more resiliant to the data
 * payload getting corrupted during transimssion.
 *
 * @param data
 *     the data to wrap into packets
 * @param data_length
 *     the length of the data being put into packets
 */
packet_collection* to_lossy(const char* data, unsigned short data_length) {
    time_t t;   
    srand((unsigned) time(&t));

    unsigned short number_of_packets = get_number_of_packets(data_length);

    // each byte in data needs to be included in the legnth, plus
    // 1 check sum per packet, plus the space for storing the packet 
    // length multiple time and the xors for the length (1 per time
    // length is represented)
    unsigned int lossy_length = ((number_of_packets) * sizeof(struct packet));

    struct packet_collection* packet_collection = 
        (struct packet_collection*) malloc(sizeof(struct packet_collection));
    packet_collection->length = lossy_length * sizeof(unsigned char);
    packet_collection->data =
        (unsigned char*) malloc(packet_collection->length * sizeof(unsigned char));
     
    memset(packet_collection->data, 0, packet_collection->length * sizeof(unsigned char));
    
    int index = 0;
    struct packet* packets = (struct packet*)&packet_collection->data[0];
    unsigned char packet_xor;
    unsigned char packet_sum;
    for (int i = 0; i < number_of_packets; ++i) {
        packet_xor = 0;
        packet_sum = 0;
        for (int j = 0; j < DATA_FORMATTER_PACKET_LENGTH; ++j) {
            if (index < data_length) {
                packets[i].bytes[j] = data[index++];
            } else {
                packets[i].bytes[j] = 0;
            }

            packet_xor ^= packets[i].bytes[j];
            packet_sum += packets[i].bytes[j];
        }
        
        packets[i].data_length = data_length;
        packets[i].index = (i % 256);
        packets[i].sum = (packet_sum + packets[i].index) % 256;
        packets[i].check_sum = (packet_xor ^ packets[i].index);
    }
    
    return packet_collection;
}

/**
 * this function increments the bit stream enumerator for a packet collection 1 bit
 *
 * @param enu
 *     the bit stream enumerator to operate on
 *
 * returns true if the enumerator is left in a valid state (meaning it has more bits) 
 */
int increment_enu_1_bit(struct packet_collection_enumerator* enu) {
    if (enu->current_byte < enu->packet_collection_length) {
        ++enu->current_bit;
        if (enu->current_bit % 8 == 0) {
            enu->current_bit = 0;
            ++enu->current_byte;
        }
        
        return enu->current_byte < enu->packet_collection_length;
    } else {
        return 0;
    }
}

/*
 * this function decrements the bit stream enumerator for a packet collection 1 bit
 *
 * @param enu
 *     the bit stream enumerator to operate on
 *
 * returns true if the enumerator was left in a valid state (meaning the byte index
 *     is greater than -1). the enumerator can be positioned before the first bit
 *     and that is considered an invalid position.
 */
int decrement_enu_1_bit(struct packet_collection_enumerator* enu) {
    if (enu->current_byte > -1) {
        --enu->current_bit;
        if (enu->current_bit < 0) {
            enu->current_bit = 7;
            --enu->current_byte;
        }
        
        return enu->current_byte > -1;
    } else {
        return 0;
    }
}

/*
 * this function increments the bit stream enumerator for a packet collection n bytes
 *
 * @param enu
 *     the bit stream enumerator to operate on
 * @param n
 *     the number of bytes to increment the bit stream operator
 *
 * returns true if there are more bytes after the increment operation completed
 */
int increment_enu_n_bytes(struct packet_collection_enumerator* enu, unsigned int n) {
    enu->current_byte += n;
    if (enu->current_byte >= enu->packet_collection_length) {
        enu->current_bit = 0;
        enu->current_byte = enu->packet_collection_length;
        return 0;
    } else {
        return 1;
    }
}

/*
 * this function decrements the bit stream enumerator for a pack collection n bytes
 *
 * @param enu
 *     the bit stream enumerator to operate on
 * @param n
 *     the number of bytes to decrement the bit stream operator
 *
 * returns true if the enumerator byte index is greater than or equal to 0 after the 
 *     decrement operation
 */
int decrement_enu_n_bytes(struct packet_collection_enumerator* enu, unsigned int n) {
    enu->current_byte -= n;
    if (enu->current_byte < 0) {
        enu->current_bit = 7;
        enu->current_byte = -1;
        return 0;
    } else {
        return 1;
    }
}

/*
 * get the current byte out of the packet collection pointed to by a bit stream enumerator
 *
 * @param packet_collection
 *     the packet collection from which to retrieve a byte
 * @param enu
 *     the bit stream enumerator pointing to a byte to retrieve
 *
 * returns the unsigned char from the packet collection pointed to be the bit stream enumerator
 */
unsigned char get_current_byte(struct packet_collection* packet_collection, struct packet_collection_enumerator* enu) {
    return 
        (packet_collection->data[enu->current_byte] >> enu->current_bit) |
        ((packet_collection->data[enu->current_byte + 1] << (8 - enu->current_bit)) & 255);
}

/*
 * this function can be used to see if a bit stream enumerator has more bits
 *
 * @param enu
 *     the bit stream enumerator to check
 *
 * returns true if the bit stream enumerator has more bits
 */
int has_next_bit(struct packet_collection_enumerator* enu) {
    return enu->current_byte < enu->packet_collection_length;
}

/* 
 * this function can be used to see if a bit stream enumerator has more bytes
 *
 * @param enu
 *     the bit stream enumerator to check
 *
 * returns true if the bit stream enumerator has more bytes
 */
int has_next_byte(struct packet_collection_enumerator* enu) {
    return 
        enu->current_byte < enu->packet_collection_length - 1 ||
        (enu->current_byte == enu->packet_collection_length - 1 && enu->current_bit == 0);
}

/**
 * this struct is used in the linked lists for each packet index to 
 * keep track of the the majority of what has come in for a particular
 * packet index.  in the event there is some corruption taking the 
 * majority vote for a packet helps fight the corruption
 */
struct node {
    unsigned char bytes[DATA_FORMATTER_PACKET_LENGTH];
    struct node* next;
};

/**
 * this function inserts a packet into the node buffer.  the node
 * buffer is an array of linked lists.  each packet gets inserted into 
 * linked list corresponding to its packet index
 *
 * @param packet
 *     the extracted packet to insert into the node_buffer
 * @param node_buffer
 *     the array of node linked lists to insert a packet into
 */
void insert_packet(struct packet* packet, struct node** node_buffer) {
    struct node* new_node = (struct node*) malloc(sizeof(struct node));
    memcpy(new_node->bytes, packet->bytes, DATA_FORMATTER_PACKET_LENGTH);    
    //packet->index is an unsigned char and 
    //node_buffer is of lenght 256 so no bounds 
    //checking necessary
    new_node->next = node_buffer[packet->index];
    node_buffer[packet->index] = new_node;
}

/**
 * this function calculates and returns the remaining bytes in a bit 
 * stream enumerator
 *
 * @param enu
 *     the bit stream enumerator to use
 *
 * returns the number of bytes remaining in the bit stream enumerator
 */
int get_remaining_bytes(struct packet_collection_enumerator* enu) {
    int rtn = enu->current_byte;
    if (enu->current_bit != 0) {
        --rtn;
    }

    rtn = enu->packet_collection_length - rtn;

    if (rtn < 0) { 
        return 0; 
    } else { 
        return rtn;
    }
}

/*
 * this function reads a specified number of packets from the packet collection and
 * places them in the buffer.
 *
 * @param packet_collection
 *     the packet collection containing the payload of packets to unpack
 * @param enu
 *     the bit stream enumerator being used for the from lossy operation
 * @param buffer
 *     the buffer to store the the read packets into
 * @param number_of_packets
 *     the number of packets to read from the packet collection
 *
 * returns true if the packets were able to be read, false if the
 *     requested number of packets exceeds the number of remaining 
 *     byets in the packet collection
 */
int read(struct packet_collection* packet_collection, struct packet_collection_enumerator* enu,
          unsigned char* buffer, int number_of_packets) {
    if (number_of_packets * sizeof(struct packet) > get_remaining_bytes(enu)) {
        return 0;
    }

    //this should not happen but if the is some bug 
    //somewhere that lets enu.current_byte be less than
    //0 this will stop from reading invalid memory
    if (enu->current_byte < 0) {
        enu->current_byte = 0;
    }

    for (int h = 0; h < number_of_packets; ++h) {
        for (int i = 0; i < sizeof(struct packet); ++i) {
            if (!has_next_byte(enu)) {
                return 0;
            }
            
            buffer[(h * sizeof(struct packet)) + i] = get_current_byte(packet_collection, enu);
            increment_enu_n_bytes(enu, 1);
        }
    }

    return 1;
}

/**
 * this function validates a packet. if a bit stream becomes corrupted by dropping
 * bits or flipping bits packets can become invalid.  the checks performed are an
 * xor comparison against a check sum and a rolling sum of the byte payload
 *
 * @param packet
 *     the packet to check for validity
 *
 * returns true if the packet passes validity checks, false otherwise
 */
int validate_packet(struct packet* packet) {
    unsigned char packet_xor = 0;
    unsigned char packet_sum = 0;
    for (int i = 0; i < DATA_FORMATTER_PACKET_LENGTH; ++i) {
        packet_xor ^= packet->bytes[i];
        packet_sum += packet->bytes[i];        
    }

    return 
        (packet->sum == (packet_sum + packet->index) % 256) &&
        (packet->check_sum ^ packet_xor ^ packet->index == 0);
}

/**
 * this method checks a packet buffer containing multiple packets to see if
 * enough of the packets checkout that it would be prudent to assume no bits
 * were dropped from the packet buffer and it is byte algined.  a bit being 
 * dropped should throw off a large number of packets as opposed a bit being
 * flipped which would only throw off one packet.
 *
 * @param buffer
 *     the buffer containing packets (a byte buffer so is stores the indivdual 
 *     bytes of all the packets sequentially)
 * @param number_of_packets
 *     the number of packets the buffer contains
 *
 * returns true if the buffer should be considered byte aligned, false otherwise
 */
int is_byte_aligned(unsigned char* buffer, int number_of_packets) {
    struct packet* packet = (struct packet*)&buffer[0];
    int number_of_valid_packets = 0;
    int number_of_matching_indexes = 0;
    int packet_index = packet[0].index;
    for (int i = 0; i < number_of_packets; ++i) {
        if (validate_packet(&packet[i])) {
            ++number_of_valid_packets;
        }

        if (i > 0) {
            if (packet[i].index - packet[i - 1].index == 1) {
                ++number_of_matching_indexes;
            }
        } else {
            if (packet[i + 1].index - packet[i].index == 1) {
                ++number_of_matching_indexes;
            }
        }
    }

    if (number_of_valid_packets >= NUMBER_OF_PACKETS_FOR_BYTE_ALIGNMENT) {
        return 1;
    } else if (number_of_packets >= 2 && number_of_valid_packets >= (number_of_packets / 2)) {
        return 1;
    } else if (number_of_packets == 1 && number_of_valid_packets == 1) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * this method handles the situations that arise when a packet is invalid.  a packet can be invalid
 * due to a bit flipping or a bit dropping.  this function tries to figure out which one of those cases
 * happened. if the buffer is not considered to be byte algined this function will shift backwards 1 bit
 * at a time to try to become byte aligned.
 *
 * @param packet_collection
 *     the packet collection holding the payload being unpacked
 * @param enu
 *     the bit stream enumerator being used for the unpacking of the packet collection
 *
 * returns some constant from the h file if byte alignment can be reached FIXED_IT will be returned,
 *     if there are not enough packets left in the packet collection TERMINATE is returned, or if 
 *     byte alignement can not be reached BYTE_ALIGNMENT_PROBLEM is returned
 */
int handle_cluster_fuck(struct packet_collection* packet_collection, struct packet_collection_enumerator* enu) {
    packet_collection_enumerator enu_hold;
    packet_collection_enumerator enu_start = *enu;

    int remaining_packets = get_remaining_bytes(enu) / sizeof(struct packet);
    if (remaining_packets <= 0) {
        return TERMINATE;
    }

    int number_of_packets = 
        remaining_packets < NUMBER_OF_NECESSARY_VALID_PACKETS ? 
        remaining_packets :
        NUMBER_OF_NECESSARY_VALID_PACKETS;

    unsigned char buffer[sizeof(struct packet) * number_of_packets];
     
    // this for loop will read a number of packets, then 
    // check them for byte alignement. if they are not aligned
    // the bit stream enumerator will be decremented one bit
    // NUMBER_OF_BITS_TO_REWIND controls how many bits we 
    // will try to rewind, which corollates to how many
    // bit we think could be dropped in a single packet
    for (int i = 0; i < NUMBER_OF_BITS_TO_REWIND; ++i) {
        // remember where we are before the read, the read is only to
        // determine byte alignement, not to actually advance the enumerator
        enu_hold = *enu;
        
        // if read is unsuccessfull that means we are at the end of the packet collection
        // alert the calling function to stop
        if (!read(packet_collection, enu, buffer, NUMBER_OF_NECESSARY_VALID_PACKETS)) {
            return TERMINATE;
        }
        
        // put the enumerator back where it was prior to the read
        *enu = enu_hold;

        // if we are byte aligned then return, the problem was
        // just a flipped bit
        if (is_byte_aligned(buffer, number_of_packets)) {
            return FIXED_IT;
        } 
        
        // if byte alignement could not be determined, move the 
        // bit stream enumerator back one to accomadate a dropped 
        // bit
        decrement_enu_1_bit(enu);
    }

    *enu = enu_start;
    return BYTE_ALIGNMENT_PROBLEM;
}

/**
 * the data lenght is a short and is embedded in every packet.  this way
 * if the length is corrupted in one packet i can make it through in another.  
 * the assumption here is than it will make it through more often than it will 
 * be corrupted.  the when unpacking the length that occurs the most times is
 * considered the actual length.  this function looks for the length that occurred
 * the most.
 *
 * @param data_lengths
 *     the unsigned short array of lengths that came through the packets that
 *     were unpacked
 * @param data_lengths_size
 *     the number of length values in the data_length array
 *
 * returns the length that occurred the most times in the data_lenghts array
 */
unsigned short calculate_data_length(unsigned short* data_lengths, int data_lengths_size) {
    if (data_lengths_size <= 0) {
        return 0;
    }

    unsigned short unique_values[256];
    memset(unique_values, 0, 256 * sizeof(unsigned short));
    int unique_counts[256];
    memset(unique_counts, 0, 256 * sizeof(int));
    int unique_values_size = 0;

    for (int i = 0; i < data_lengths_size; ++i) {
        int contains = 0;
        for (int j = 0; j < unique_values_size; ++j) {
            if (unique_values[j] == data_lengths[i]) {
                contains = 1;
                ++unique_counts[j];
                break;
            }
        }

        // we only look at the first 256 uniquely occuring
        // length values.  that is way more than should
        // ever occur but if for some reason more than that
        // occur (or even anything remotely close to that occur)
        // something has gone wrong.  in practice there should be a 
        // relatively small number of lengths that occur (less than 10
        // probably)    
        if (!contains && unique_values_size < 256) {
            unique_values[unique_values_size] = data_lengths[i];
            unique_counts[++unique_values_size] = 1;
        }
    }

    // look through the unique value histogram and find
    // the value that occurred the most
    unsigned short max_index = 0;
    int max_value = unique_counts[0];
    for (int i = 1; i < unique_values_size; ++i) {
        if (unique_counts[i] > max_value) {
            max_index = i;
            max_value = unique_counts[i];
        }
    }

    // return the unique value whose count was maximal
    return unique_values[max_index];
}

/**
 * this function reclaims memory allocated in the linked lists of the node_buffer 
 * array of linked lists
 *
 * @param node_buffer
 *     a node buffer used when unpacking a lossy data stream
 */
void free_nodes(struct node** node_buffer) {
    for (int i = 0; i < 256; ++i) {
        struct node* node = node_buffer[i];
        struct node* hold;
        while (node != NULL) {
            hold = node->next;
            free(node);
            node = hold;
        }
    }
}

/**
 * this function looks over the supplied linked list and logs how many unique
 * values were occurred in the data stream for the character at a position of the
 * data stream.  return the character that occurs the most in the linked list
 *
 * @param node
 *     a linked list of packets that occurred while unpacking a packet collection
 * @param bytes_offset
 *     the bytes position to look at in the packet linked list, this corresponds to
 *     the character being saught for the data stream being unpacked.
 */
unsigned char calculate_value(struct node* node, int bytes_offset) {
    if (node == NULL) {
        return 0;
    }

    unsigned char unique_values[256];
    int unique_values_size = 0;
    int unique_counts[256];
    
    // loop through the values that occurred for this character position 
    // in the data stream and keep a histogram of how many times a value 
    // occurred.  ultimately we want to return the value that occured the most
    while (node != NULL) {
        int contains = 0;
        for (int j = 0; j < unique_values_size; ++j) {
            if (unique_values[j] == node->bytes[bytes_offset]) {
                contains = 1;
                ++unique_counts[j];
                break;
            }
        }

        // only look at the 1st 256 unique values for this character in the 
        // data stream.  that should be way more than enough
        if (!contains && unique_values_size < 256) {
            unique_values[unique_values_size] = node->bytes[bytes_offset];
            unique_counts[++unique_values_size] = 1;
        }

        node = node->next;
    }

    // look through the unique values histogram and find which value occured
    // the most number of times
    unsigned short max_index = 0;
    unsigned char max_value = unique_counts[0];
    for (int i = 1; i < unique_values_size; ++i) {
        if (unique_counts[i] > max_value) {
            max_index = i;
            max_value = unique_counts[i];
        }
    }
    
    //return the value that occurred the most number of times
    return unique_values[max_index];
}

/**
 * this function iterates through the node buffer and gets the character that occurred the most
 * at every position in the output data stream while unpacking the packet collection
 *
 * @param node_buffer
 *     the node buffer the packets were stored in
 * @param data_out
 *     the buffer to store the output data stream to
 * @param data_out_length
 *     the length of the data_out buffer
 */
void calculate_data_out(struct node** node_buffer, unsigned char* data_out, unsigned int data_out_length) {
    unsigned char number_of_packets = get_number_of_packets(data_out_length);
    for (int i = 0; i < number_of_packets; ++i) {
        for (int j = 0; j < DATA_FORMATTER_PACKET_LENGTH; ++j) {
            unsigned short data_index = (i * DATA_FORMATTER_PACKET_LENGTH) + j;
            if (data_index < data_out_length && node_buffer[i] != NULL) {
                data_out[data_index] = calculate_value(node_buffer[i], j);
            }
        }
    }
}

/**
 * this function unpacks a packet_collection and produces the original input data that was packed into it. 
 * it is assumed the packet collection could have become corrupted and this function will try to not let any
 * corruption stop it from reconstructing the original data.
 *
 * @param packet_collection
 *     a packet collection that was deserliazed from a jpeg image
 * @param data_out
 *     a pointer to set to point at the resulting, unpacked data from the packet collection
 * @param data_out_length
 *     this variable will be set to the length of the data_out
 */
void from_lossy(struct packet_collection* packet_collection, unsigned char** data_out, unsigned int* data_out_length) {
    struct node* node_buffer[256];
    memset(node_buffer, 0, 256 * sizeof(struct node*));

    unsigned short data_lengths[256];
    memset(data_lengths, 0, 256 * sizeof(unsigned short));    
    int data_lengths_size = 0;

    unsigned char* payload;
    struct packet_collection_enumerator enu;
    enu.current_bit = 0;
    enu.current_byte = 0;
    enu.packet_collection_length = packet_collection->length;

    unsigned char buffer[sizeof(struct packet)];
    struct packet* packet = (struct packet*)&buffer[0];

    while (read(packet_collection, &enu, buffer, 1)) {
        int is_valid_packet = validate_packet(packet);
        
        // if packet is not valid attempt to figure out why
        // if a bit was dropped handle_cluster_fuck will position
        // the enu to account for the dropped bit
        if (!is_valid_packet) {
            int c_fuck_return = handle_cluster_fuck(packet_collection, &enu);
            if (c_fuck_return == TERMINATE) {
                break;
            }
        } else {
            // if the packet is valid put it in the node buffer at the index
            // of the packet (packet->index). 
            insert_packet(packet, node_buffer);
            // keep track of the data lengths we get from the packets
            // data lengths is the size in chars of the original encoded
            // data
            if (data_lengths_size < 256) {
                data_lengths[++data_lengths_size] = packet->data_length;
            }
        }
    }

    // if we did not get any packets that were valid
    if (data_lengths_size < 1) {
        *data_out = NULL;
        *data_out_length = 0;
        free_nodes(node_buffer);
        return;
    }

    unsigned short data_length = calculate_data_length(data_lengths, data_lengths_size);    
    // this check prevent allocating rediculous amounts of memory as a result of an
    // improperly unpacked packet collection.  we should not be encoding anything in 
    // the lossy format that exceeds this amount of memory
    if (data_length > 4096) {
        *data_out = NULL;
        *data_out_length = 0;
        free_nodes(node_buffer);
        return;
    }
    
    // allocate the memory for the returned original data
    unsigned char* data = (unsigned char*) malloc(data_length * sizeof(unsigned char));
    memset(data, 0, data_length * sizeof(unsigned char));

    // reconstruct the original data from the unpacked packets
    calculate_data_out(node_buffer, data, data_length);
    
    //free memory, be responsible
    free_nodes(node_buffer);

    *data_out = data;
    *data_out_length = data_length;
}
