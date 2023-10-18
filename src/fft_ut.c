#include <stdlib.h>
#include <stdio.h>
#include "fft.h"

int main(int argc, char ** argv)
{
    if(argc != 1)
    {
        printf("%s does not use any command line arguments\n", argv[0]);
    }
  fft_ut();
  return 0;
}
