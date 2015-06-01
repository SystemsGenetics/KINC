#include "error.h"

/**
 * Prints an error message to the terminal throws and terminates the program.
 *
 * @param char * message
 *   The message to print to the terminal.
 */
void handle_error(char * message) {
  printf("ERROR: %s\nTerminating.", message);
  exit(-1);
}
/**
 * Prints a warning message to the terminal but allows the program to continue.
 *
 * @param char * message
 *   The message to print to the terminal.
 */
void handle_warning(char * message) {
//  printf("WARNING: %s\n.", message);
}
