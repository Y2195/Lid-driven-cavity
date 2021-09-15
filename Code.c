#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main()
{
	int i, j, k = 0;
	int m, n, l;
	double dx, dy, beta, Err_psi, Err_w;
	double L = 1.0;    			  // Length of square cavity
	double Re = 400.0;			  // Reynolds number
	m = 129;
	n = 129;
	dx = L/(m-1);
	dy = L/(n-1);
	l = (m-2)*(n-2);
	beta = dx/dy;
	double u[m][n], v[m][n];      // Velocity
	double psi[m][n]; 			  // Stream function
	double w[m][n];				  // Vorticity
	double psi_old[m][n], w_old[m][n];
	double x, y;
	// Initialization
	for(j=0; j<n; j++)
	{
		for(i=0; i<m; i++)
		{
			if(j==0)
			{
				u[i][j]=0;
				v[i][j]=0;
				psi[i][j]=0;
			}
			else if(j==(n-1))
			{
				u[i][j]=1;
				v[i][j]=0;
				psi[i][j]=0;
			}
			else if(i==0)
			{
				u[i][j]=0;
				v[i][j]=0;
				psi[i][j]=0;
			}
			else if(i==(m-1))
			{
				u[i][j]=0;
				v[i][j]=0;
				psi[i][j]=0;
			}
			else
			{
				u[i][j]=0;
				v[i][j]=0;
				psi[i][j]=0;
			}
		}
	}
	for(j=0; j<n; j++)
	{
		for(i=0; i<m; i++)
		{
			if(j==0)
			{
				w[i][j]=(2/pow(dy,2))*(psi[i][j]-psi[i][j+1]);
			}
			else if(j==(n-1))
			{
				w[i][j]=((2/pow(dy,2))*(psi[i][j]-psi[i][j-1]))-((2*u[i][j])/dy);
			}
			else if(i==0)
			{
				w[i][j]=(2/pow(dx,2))*(psi[i][j]-psi[i+1][j]);
			}
			else if(i==(m-1))
			{
				w[i][j]=(2/pow(dx,2))*(psi[i][j]-psi[i-1][j]);
			}
			else
			{
				w[i][j] = 0.0;
			}
		}
	}
	// Gauss Siedel Method
	do
	{
		for(j=0; j<n; j++)
		{
			for(i=0; i<m; i++)
			{
				psi_old[i][j] = psi[i][j];
				w_old[i][j] = w[i][j];
			}
		}
		// Solving for stream function and vorticity
		for(j=1; j<(n-1); j++)
		{
			for(i=1; i<(m-1); i++)
			{
				psi[i][j] = (0.5/(1+pow(beta,2)))*((psi[i+1][j]+psi[i-1][j])+(pow(beta,2)*(psi[i][j+1]+psi[i][j-1]))+(pow(dx,2)*w[i][j])); 
			}
		}
		for(j=1; j<(n-1); j++)
		{
			for(i=1; i<(m-1); i++)
			{
			//	w[i][j] = (0.5/(1+pow(beta,2)))*(w[i+1][j]+w[i-1][j]+(pow(beta,2)*(w[i][j+1]+w[i][j-1]))-((0.25*beta*Re)*(w[i+1][j]-w[i-1][j])*(psi[i][j+1]-psi[i][j-1]))+((0.25*beta*Re)*(w[i][j+1]-w[i][j-1])*(psi[i+1][j]-psi[i-1][j])));
				w[i][j]=(0.5/(1+pow(beta,2)))*(((1-(0.25*Re*beta*(psi[i][j+1]-psi[i][j-1])))*w[i+1][j])+((1+(0.25*Re*beta*(psi[i][j+1]-psi[i][j-1])))*w[i-1][j])+(pow(beta,2)*(1+(0.25*(Re/beta)*(psi[i+1][j]-psi[i-1][j])))*w[i][j+1])+(pow(beta,2)*((1-(0.25*(Re/beta)*(psi[i+1][j]-psi[i-1][j])))*w[i][j-1])));
			}
		}
		// Updating vorticity boundary condition
		for(j=0; j<n; j++)
		{
			for(i=0; i<m; i++)
			{
				if(j==0)
				{
					w[i][j]=(2/pow(dy,2))*(psi[i][j] - psi[i][j+1]);
				}
				else if(j == (n-1))
				{
					w[i][j] = (2/pow(dy,2))*(psi[i][j] - psi[i][j-1]) - ((2*u[i][j])/dy);
				}
				else if(i == 0)
				{
					w[i][j] = (2/pow(dx,2))*(psi[i][j] - psi[i+1][j]);
				}
				else if(i == (m-1))
				{
					w[i][j] = (2/pow(dx,2))*(psi[i][j] - psi[i-1][j]);
				}
			}
		}
		// Error calculation
		Err_psi = 0;
		Err_w = 0;
		for(j=1; j<(n-1); j++)
		{
			for(i=1; i<(m-1); i++)
			{
				Err_psi = Err_psi + pow((psi[i][j]-psi_old[i][j]),2);
				Err_w = Err_w + pow((w[i][j]-w_old[i][j]),2);
			}
		}
		Err_psi = sqrt(Err_psi/l);
		Err_w = sqrt(Err_w/l);
		//printf("%d \t %lf \t %lf\n", k, Err_psi, Err_w);
		k++;
	}
	while(Err_psi > pow(10,-4) || Err_w > pow(10,-4));
	// update velocities
	for(j=1; j<(n-1); j++)
	{
		for(i=1; i<(m-1); i++)
		{
			u[i][j] = (0.5/dy)*(psi[i][j+1]-psi[i][j-1]);
			v[i][j] = (-0.5/dx)*(psi[i+1][j]-psi[i-1][j]);
		}
	}
	FILE *fp, *fp1, *fp2;
	fp = fopen("lid_driven_cavity.txt","w");
	fp1 = fopen("u_along_vertical.txt","w");
	fp2 = fopen("v_along_horizontal.txt","w");
	printf("Node-x\t\t Node-y \t Stream function \t  vorticity \t u \t v\n");
	for(j=0; j < n; j++)
	{
		for(i=0; i < m; i++)
		{
			x=i*dx;
			y=j*dy;
			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, u[i][j], v[i][j], psi[i][j], w[i][j]);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x, y, psi[i][j], w[i][j],u[i][j], v[i][j]);
		}
	}
	printf("\nTotal iteration = %d\n", k);	
	for(j=0; j<n; j++)
	{
		y = j*dy;
		fprintf(fp1,"%lf \t %lf\n", y, u[64][j]);
	}
	for(i=0; i<m; i++)
	{
		x = i*dx;
		fprintf(fp2,"%lf \t %lf\n", x, v[i][64]);
	}
	fclose(fp);
	return 0;
}
