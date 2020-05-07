#include <stdlib.h>
#include <stdio.h>
#include "fft.h"
#include "fft.c"

int main(int argc, char ** argv)
{
  /* Try this when $HOME/.config/deconwolf/ does not exist
   * and when it does ... 
   * could also test it when that dir isn't writeable.
   * What to do on OSX?
   * */

  int nThreads = 4;
  char * swf = get_swf_file_name(nThreads);
  printf("swf = '%s'\n", swf);
  free(swf);

  return 0;
}
