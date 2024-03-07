#pragma once

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cairo.h>
#include <cairo-svg.h>

#include "fim.h"
#include "ftab.h"
#include "fim_tiff.h"
#include "dw_version.h"
#include "dw_util.h"

int dw_render(int argc, char ** argv);
