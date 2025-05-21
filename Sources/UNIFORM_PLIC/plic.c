#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uniform.h"
#include "nrutil.h"

void plic(void)
{
    // printf("PLIC\n");

    // printf("PLIC X\n");
    clear_vof(temp_vof);
    plicx(vof, temp_vof);
    correct_vof(temp_vof);

    // printf("PLIC Y\n");
    clear_vof(vof);
    plicy(temp_vof, vof);
    correct_vof(vof);

    // printf("PLIC done\n");
}