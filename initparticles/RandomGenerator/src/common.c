#include <stdio.h>
#include <stdlib.h>
#include "common.h"


/* ! calloc with error handling ! 8 ! */
void *common_calloc(size_t count, size_t size){
  void *ptr = calloc(count, size);
  if(ptr == NULL){
    fprintf(stderr, "memory allocation error\n");
    exit(EXIT_FAILURE);
  }
  return ptr;
}

/* ! free with error handling ! 4 ! */
void common_free(void *ptr){
  free(ptr);
  ptr = NULL;
}

