#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fim.h"
#include "fim_tiff.h"

void usage(int argc, char ** argv)
{
  printf("Make a maximum intensity projection of uint16_t multi-page tiff file\n");
  printf("Usage: %s input.tif [output.tif]\n", argv[0]);
}

int main(int argc, char ** argv)
{
  if(argc < 2)
  {
    usage(argc, argv);
    exit(1);
  }

  char * inFile = argv[1];
  char * outFile = NULL;
  if(argc > 2)
  {
    outFile = argv[2];
  } else {
    outFile = malloc(strlen(inFile) + 10);
    sprintf(outFile, "max_%s", inFile);
  }

  printf("Input: %s\n", inFile);
  printf("Output: %s\n", outFile);

  fim_tiff_maxproj(inFile, outFile);

}

