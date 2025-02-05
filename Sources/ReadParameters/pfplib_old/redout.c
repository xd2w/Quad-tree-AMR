#include <stdio.h>

int redout()
{
	return (!isatty(fileno(stdout)));
}
