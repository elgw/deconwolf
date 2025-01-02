/* To test that the library is installed and can be used
 *
 * Compile with:
 * gcc minimal_example.c -ltrafo -o minimal_example -Wall -Wextra --pedantic
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <trafo.h>

int main (int argc, char ** argv)
{
    if(argc != 1)
    {
        printf("%s does not take any command line arguments\n", argv[0]);
    }

    /* This information is in the header file */
    printf("Using trafo v.%d.%d.%d\n",
           TRAFO_VERSION_MAJOR, TRAFO_VERSION_MINOR, TRAFO_VERSION_PATCH);

    printf("Trying to fit using invalid settings "
           "-- expecting an error message\n");

    trafo_settings C = {0};

    /* This function is in the library */
    trf * T = trafo_fit(&C);
    assert(T == NULL);

    return EXIT_SUCCESS;
}
