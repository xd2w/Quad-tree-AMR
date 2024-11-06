#include <stdio.h>
#include <stdlib.h>
#include "ftt.h"
#include "pfplib.h"

// balance the FTT tree at all levels
void octPropagation(void)
{
  int iprop, it, nbPropagation;

  nbPropagation = 1;
  ifetch("nbPropagation", &nbPropagation);

  for (iprop = 0; iprop < nbPropagation; iprop++)
  {
    for (it = maxLevel - 1; it > 0; it--)
    {
      balanceFTT(it, 0);
    }
  }
}
