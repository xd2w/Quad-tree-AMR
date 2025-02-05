#include <stdio.h>

int redin()
{
	return (!isatty(fileno(stdin)));
}
