#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"

//CONFIGURACAO GLOBAL DO LABORATORIO 
int    cfg_A = 2;
int    cfg_B = 1;
int    cfg_C = 1;

double cfg_K       = K;
double cfg_omega_A = OMEGAGS;
double cfg_omega_B = OMEGA;
double cfg_omega_C = OMEGA;
int    cfg_NU1     = NU1;
int    cfg_NU2     = NU2;
double cfg_tol_global = T;
double cfg_tol_B   = TOLERANCIA_B;
double cfg_tol_C   = TOLERANCIA_C;
int    cfg_N       = NUM;
const char *nome_A(int a)
{
    switch (a)
    { 
	    case 1: return "Jacobi"; 
	    case 2: return "Gauss-Seidel"; 
    }
    return "?";
}
const char *nome_B(int b)
{
    switch (b) 
    {
        case 1: return "LU";
        case 2: return "LDL";
        case 3: return "Gauss-Jordan";
        case 4: return "Jacobi";
        case 5: return "Gauss-Seidel";
    }
    return "?";
}

const char *nome_C(int c)
{
    switch (c) 
    {
        case 1: return "Jacobi";
        case 2: return "Gauss-Seidel";
        case 3: return "LU";
        case 4: return "LDL";
        case 5: return "Gauss-Jordan";
    }
    return "?";
}



