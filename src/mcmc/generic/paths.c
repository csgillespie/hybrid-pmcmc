#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "paths.h"

char *createPath(char *dir, char *file) 
{
  char *path;
  path = (char *)malloc((strlen(dir) + strlen(file) + 1)*sizeof(char));

  strcpy(path, dir);
  strcat(path, file);
  return(path);
}

char *addInputPath(char *dir) {
  return(createPath("../input/", dir));
}
