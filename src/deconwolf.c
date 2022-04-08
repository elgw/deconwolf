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
        #ifdef __dw_otsu_h__
        if(strcmp(argv[1], "otsu") == 0)
        {
            return dw_otsu(argc-1, argv+1);
        }
        #endif
        #ifdef __dw_dots_h__
        if(strcmp(argv[1], "dots") == 0)
        {
            return dw_dots(argc-1, argv+1);
        }
        #endif
    }

    dw_opts * s = dw_opts_new(); /* Load default settings and initialize */
    dw_argparsing(argc, argv, s); /* Parse command line */
    return dw_run(s); /* And do the job */
}
