/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Command line interface for deconwolf */

/* Extra modules can be enabled by un-commenting in
 * the header file. */
#include "dw.h"

int main(int argc, char ** argv)
{
    // setlocale(LC_ALL, NULL);
    if(argc > 1)
    {
        if(strcmp(argv[1], "maxproj") == 0)
        {
            return dw_tiff_max(argc-1, argv+1);
        }
        if(strcmp(argv[1], "imshift") == 0)
        {
            return dw_imshift(argc-1, argv+1);
        }
        if(strcmp(argv[1], "merge") == 0)
        {
            return dw_tiff_merge(argc-1, argv+1);
        }

        if(strcmp(argv[1], "nuclei") == 0)
        {
#ifdef dw_module_nuclei
            return dw_nuclei(argc-1, argv+1);
#else
            fprintf(stderr, "dw was built without the 'nuclei' module\n");
            exit(EXIT_FAILURE);
#endif
        }

        if(strcmp(argv[1], "dots") == 0)
        {
#ifdef dw_module_dots
            return dw_dots(argc-1, argv+1);
#else
            fprintf(stderr, "dw was built without the 'dots' module\n");
            exit(EXIT_FAILURE);
#endif
        }

        if(strcmp(argv[1], "psf") == 0)
        {
#ifdef dw_module_psf
            return dw_psf_cli(argc-1, argv+1);
#else
            fprintf(stderr, "dw was built without the 'psf' module\n");
            exit(EXIT_FAILURE);
#endif
        }

        if(strcmp(argv[1], "psf-STED") == 0)
        {
#ifdef dw_module_psf_sted
            return dw_psf_sted_cli(argc-1, argv+1);
#else
            fprintf(stderr, "dw was built without the 'psf-STED' module\n");
            exit(EXIT_FAILURE);
#endif
        }

        if(strcmp(argv[1], "noise1") == 0)
        {
            return sparse_preprocess_cli(argc-1, argv+1);
        }

    }

    dw_opts * s = dw_opts_new(); /* Load default settings and initialize */
    dw_argparsing(argc, argv, s); /* Parse command line */
    return dw_run(s); /* And do the job */
}
