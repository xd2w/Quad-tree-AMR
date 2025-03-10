#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

#define pi 3.14159265

// void computeVelocity(void)
// {
// }

Real computeVX(Real x, Real y)
{
    // return 1;
    return -y;
    // return -2 * sin(pi * x) * sin(pi * x) * sin(pi * y) * cos(pi * y);
}

Real computeVY(Real x, Real y)
{
    // return 1;
    return x;
    // return -2 * cos(pi * x) * sin(pi * x) * cos(pi * y) * cos(pi * y);
}
