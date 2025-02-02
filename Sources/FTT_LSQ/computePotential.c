#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "pfplib.h"
#include "ftt.h"
#include "nrutil.h"
#include "box.h"
#include "pfplib.h"

#define pi 3.14159265

Real computeVX(Real x, Real y)
{
    if (time > 62.51)
    {
        return -2 * sin(pi * x) * sin(pi * x) * sin(pi * y) * cos(pi * y);
    }
    return 2 * sin(pi * x) * sin(pi * x) * sin(pi * y) * cos(pi * y);
}

Real computeVY(Real x, Real y)
{
    if (time > 62.51)
    {
        return -2 * cos(pi * x) * sin(pi * x) * cos(pi * y) * cos(pi * y);
    }
    return 2 * cos(pi * x) * sin(pi * x) * cos(pi * y) * cos(pi * y);
}

// Real computeVX(Real x, Real y)
// {
//     return x;
// }

// Real computeVY(Real x, Real y)
// {
//     return -y;
// }