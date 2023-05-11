#include <stdlib.h>
#include <stdio.h>
#include "create_array.h"

void *create_array(int size) {
  void *x = malloc(size);
  if (x == NULL) {
	printf("Not enough memory");
	exit(1);
  }
  return x;
}


