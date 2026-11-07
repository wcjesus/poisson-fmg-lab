#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"
#include "problema.h"

double soma (double x)
{
	double s=0.0;
	double freq = cfg_K;

	for (int j=1; j<=EXP; j++)
	{
		s += sin(freq*M_PI*x);
		freq*=cfg_K;
	}
	return s;
}

double soma2 (double x)
{
	double s=0.0;
	double freq =cfg_K;

	for (int i=1; i<=EXP; i++)
	{
		s+= (double)freq*freq*sin(freq*M_PI*x);
		freq *=cfg_K;
	}
	return s;
}

double uex(double x, double y)
{
    return (soma(x) * soma(y));
}

double func(double x, double y)
{
    return M_PI*M_PI*(soma2(x)*soma(y)+soma(x)*soma2(y));
}

double contorno(double x, double y) { return uex(x, y); }

