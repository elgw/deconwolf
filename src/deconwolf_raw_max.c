#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage(int argc, char ** argv)
{
  printf("Usage: %s inputfile M N P\n", argv[0]);
}

int main(int argc, char ** argv)
{
  if(argc < 5)
  {
    usage(argc, argv);
  }

  size_t M = atol(argv[2]);
  size_t N = atol(argv[3]);
  size_t P = atol(argv[4]);

  char * ofilename = malloc(strlen(argv[1]) + 10);
  sprintf(ofilename, "%s.maxz", argv[1]);

  FILE * fin = fopen(argv[1], "r");
  FILE * fout = fopen(ofilename, "w");

  float * maxline = malloc(M*sizeof(float));
  float * buf = malloc(M*sizeof(float));

  // Process one plane at a time
  // Never using more than one column of memory
  printf("Assuming size %zu %zu %zu\n", M, N, P);
  printf("Writing to %s\n", ofilename);
  printf("Processing ...\n");
  for(size_t nn = 0; nn<N; nn++)
  {
    if(nn%10 == 0){
    printf("\r %zu / %zu", nn, N);
    }
    memset(maxline, 0, M*sizeof(float));
    for(size_t pp = 0; pp<P; pp++)
    {
      size_t colpos = (nn*M + pp*M*N)*sizeof(float);
      fseek(fin, colpos, SEEK_SET);
      size_t ok = fread(buf, M*sizeof(float), 1, fin);
      if(ok != 1)
      {
        printf("An error occured while reading! Check the image dimensions.\n");
        exit(1);
      }
      for(size_t kk = 0; kk<M; kk++)
      {
        if(buf[kk] > maxline[kk])
        {
          maxline[kk] = buf[kk];
        }
      }
      fseek(fout, nn*M*sizeof(float), SEEK_SET);
      fwrite(maxline, M*sizeof(float), 1, fout);
    }
  }
  printf("\nDone\n");
  
  free(buf);
  free(maxline);
  fclose(fin);
  fclose(fout);

  printf("%% In MATLAB\n:");
  printf("f = fopen('%s', 'rb')\n", ofilename);
  printf("I = fread(f, inf, 'float32'); fclose(f);\n");
  printf("I = reshape(I, [%zu, %zu]);\n", M, N);
  printf("imagesc(I), colormap gray, axis image\n");

  return 0;
}

