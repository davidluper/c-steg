//
//  Utils.cpp
//  c-steg
//

#include "utils.h"
#include <string.h>

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
void Utils::getFileNameWithoutExtension(const char* file, char* out_buffer) {
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

