#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "configuracao.h"
#include "suavizacao.h"

/* 
 (A) Gauss-Seidel com peso omega.                                                                                                                 
 IMPORTANTE: GS *nao* e Jacobi. GS usa os valores ja atualizados nesta   
 varredura (u[i-1][j] e u[i][j-1] do estado novo). Por isso a funcao     
 precisa ser chamada escrevendo em u *in-place* -- nao em u_novo. Quem   
 faz a chamada correta e smooth_step().                                  
                                                                         
 Formula: u_i = (1-omega) u_i^velho + omega u_i^GS                          
 */



double gaussseidel(double **u, double **f, int i, int j, double h, double omega)
{
    // u[i-1][j] e u[i][j-1] JA ESTAO ATUALIZADOS (varredura crescente) 
    double u_anterior = u[i][j];
    double u_gs = 0.25 * (h*h*f[i][j] + u[i-1][j] + u[i+1][j]
                                       + u[i][j-1] + u[i][j+1]);
    return (1-omega)*u_anterior+omega*u_gs;
}

double jacobi_ponderado(double **u, double **f, int i, int j, double h, double omega)
{
    double u_jacobi = 0.25 * (h*h*f[i][j] + u[i-1][j] + u[i+1][j]
                                           + u[i][j-1] + u[i][j+1]);
    return (1.0 - omega) * u[i][j] + omega * u_jacobi;
}


/* 
  smooth_step: aplica 1 passo de Jacobi (method=1) ou GS (method=2).      
                                                                          
 - Jacobi: le de u, escreve em u_novo (buffer separado).                  
 - GS:     escreve in-place em u, e ao final copia para u_novo,           
           deixando o resultado em ambos (preserva a interface).          
                                                                          
 Apos a chamada, u_novo sempre contem o resultado do passo.               
  */
double **smooth_step(double **u, double **u_novo, double **f,
                     int nx, int ny, double h, int method, double omega)
{
    int i, j;

    if (method == 1) {
        // Jacobi: le u, escreve interior de u_novo 
        for (i = 1; i < nx; i++)
            for (j = 1; j < ny; j++)
                u_novo[i][j] = jacobi_ponderado(u, f, i, j, h, omega);

        // Bordas de u_novo (o buffer devolvido precisa ter contorno valido) 
        for (i = 0; i <= nx; i++) { u_novo[i][0] = u[i][0]; u_novo[i][ny] = u[i][ny]; }
        for (j = 0; j <= ny; j++) { u_novo[0][j] = u[0][j]; u_novo[nx][j] = u[nx][j]; }

        return u_novo;
    } else {
        // GS (SOR): in-place em u; bordas de u ja estao corretas 
        for (i = 1; i < nx; i++)
            for (j = 1; j < ny; j++)
                u[i][j] = gaussseidel(u, f, i, j, h, omega);

        return u;
    }
}


//  para o smoother A (chamado por suaviza/suavizacoes).          
void suaviza(double **u, double **u_novo, double **f, int nx, int ny, double h)
{
    //A 
    smooth_step(u, u_novo, f, nx, ny, h, cfg_A, cfg_omega_A);
}

void suavizacoes(double **u, double **u_novo, double **f,
                 int nx, int ny, double h, int numero)
{
    double **u_original = u;

    for (int k=0;k<numero;k++)
    {
	    double **res = smooth_step(u, u_novo,f,nx,ny,h,cfg_A,cfg_omega_A);

	    if (res==u_novo)
	    {
		    double **temp=u;
		    u=u_novo;
		    u_novo=temp;
	    }
    }
    if (u!=u_original)
    {
	    for (int i=0;i<=nx;i++)
		    for (int j=0;j<=ny;j++)
			    u_original[i][j]=u[i][j];
    }
}

