#include <stdlib.h>
#include <stdio.h>
#include "fim.h"
#include "dw_util.h"

int main(int argc, char ** argv)
{
    if(argc > 1)
    {
        printf("WARNING: %s does not use any command line arguments\n", argv[0]);
    }
    fim_ut();
    return EXIT_SUCCESS;
}
