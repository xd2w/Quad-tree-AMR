#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

int redin()
{
	return (!isatty(fileno(stdin)));
}
