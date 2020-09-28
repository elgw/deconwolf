#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fim.h"
#include "fim_tiff.h"
#include "dw_version.h"

void usage(int argc, char ** argv)
{
  printf("Make a maximum intensity projection of uint16_t multi-page tiff file\n");
  printf("Usage: %s input1.tif input2.tif ... \n", argv[0]);
  printf("Example: %s dw*tif\n", argv[0]);
  printf("Part of deconwolf %s", deconwolf_version);
}

int file_exist(char * fname)
{
  if( access( fname, F_OK ) != -1 ) {
    return 1; // File exist
  } else {
    return 0;
  }
}

int main(int argc, char ** argv)
{
  if(argc < 2)
  {
    usage(argc, argv);
    exit(1);
  }

  char * inFile = NULL;
  char * outFile = NULL;

  for(int ff = 1; ff<argc; ff++)
  {

    inFile = argv[ff];
    outFile = malloc(strlen(inFile) + 10);
    sprintf(outFile, "max_%s", inFile);
    if(file_exist(outFile))
    {
      printf("%s exists, skipping.\n", outFile);
    } else {
      printf("%s -> %s\n", inFile, outFile);
      fim_tiff_maxproj(inFile, outFile);
    }
    free(outFile);
  }

}

