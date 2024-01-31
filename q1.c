#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main()
{
// Creating File
    FILE *u1;
    FILE *v1;
    u1=fopen("u_1.plt","w");
	v1=fopen("v_1.plt","w");
// declaration of variables
	int i,j,m,n,k=0;
	double sum1,sum2,delx,dely,e1=1.0,e2=1.0,Tn,Re;
	printf("Enter the values For m, n and Re \n");
	scanf("%d%d%lf",&m,&n,&Re);  //taking input
	delx=1.0/(m-1);
	dely=1.0/(n-1);
	double fe,fn,de,dn,ff,c[m][m];
	double ae[m][n],aw[m][n],as[m][n],an[m][n],u[m][n],v[m][n],un[m][n],vn[m][n],p[m][n],pn[m][n],be[m][n],bw[m][n],bs[m][n],bn[m][n],a[m][n],b[m][n];
	double ce[m][n],cw[m][n],cs[m][n],cn[m][n];
    de=dely/(delx*Re);
    dn=delx/(dely*Re);
    fe=dely/2;
    fn=delx/2;

//initializing the interior points
	for(j=1;j<n-2;j++)
    {
        for(i=1;i<m-2;i++)
        {
          pn[i][j]=0;
        }
    }
    for(j=1;j<n-2;j++)
    {
        for(i=1;i<m-1;i++)
        {
          u[i][j]=0;
        }
    }
    for(j=1;j<n-1;j++)
    {
        for(i=1;i<m-2;i++)
        {
          v[i][j]=0;
        }
    }
//boundary conditions
	for(i=0;i<m-1;i++)
    {
        pn[i][0]=0.0;
        pn[i][n-1]=0.0;
        v[i][0]=vn[i][0]=0;
        v[i][n-1]=vn[i][n-1]=0;
    }
    for(i=0;i<m;i++)
    {
        u[i][0]=un[i][0]=0;
        u[i][n-1]=un[i][n-1]=1;
    }

    for(j=0;j<n-1;j++)
    {
        pn[0][j]=0.0;
        pn[m-1][j]=0.0;
        u[0][j]=un[0][j]=0;
        u[m-1][j]=un[m-1][j]=0;
    }
    for(j=0;j<n;j++)
    {
        v[0][j]=vn[0][j]=0;
        v[m-1][j]=vn[m-1][j]=0;
    }
   for(j=0;j<n-1;j++)
    {
        for(i=0;i<m-1;i++)
        {
          p[i][j]=0;
        }
    }

     do
    {
        sum1=0;
        sum2=0;
        for(j=0;j<n-1;j++)
   {
       a[0][j]=a[m-1][j]=1.0;
   }
    for(i=0;i<m-1;i++)
   {
       b[i][0]=b[i][n-1]=1.0;
   }

        for(j=1;j<n-1;j++)  //storing the values to calculate error for vorticity and stream function
   {
        for(i=1;i<m-1;i++)
        {
            ae[i][j]=de-(u[i][j]+u[i+1][j])*fe/2;
            an[i][j]=dn-(v[i][j]+v[i+1][j])*fn/2;
            aw[i][j]=de+(u[i][j]+u[i-1][j])*fe/2;
            as[i][j]=dn+(v[i+1][j-1]+v[i][j-1])*fn/2;
            a[i][j]=ae[i][j]+an[i][j]+aw[i][j]+as[i][j];
            be[i][j]=de-(u[i][j]+u[i][j+1])*fe/2;
            bn[i][j]=dn-(v[i][j]+v[i][j+1])*fn/2;
            bw[i][j]=de+(u[i-1][j+1]+u[i-1][j])*fe/2;
            bs[i][j]=dn+(v[i][j]+v[i][j-1])*fn/2;
            b[i][j]=be[i][j]+bn[i][j]+bw[i][j]+bs[i][j];

        }
   }
   for(j=1;j<n-2;j++) // Calculating Velocity (u,v)terms
   {
        for(i=1;i<m-1;i++)
        {
            un[i][j]=(ae[i][j]*un[i+1][j]+aw[i][j]*un[i-1][j]+an[i][j]*un[i][j+1]+as[i][j]*un[i][j-1]-dely*(p[i+1][j]-p[i][j]))/a[i][j];

        }
   }
   for(j=1;j<n-1;j++)
   {
        for(i=1;i<m-2;i++)
        {
            vn[i][j]=(be[i][j]*vn[i+1][j]+bw[i][j]*vn[i-1][j]+bn[i][j]*vn[i][j+1]+bs[i][j]*vn[i][j-1]-delx*(p[i][j+1]-p[i][j]))/b[i][j];

        }
   }
//calculation for stream function (Resolving Pressure)
         for(j=1;j<n-2;j++)
   {
        for(i=1;i<m-2;i++)
        {
            ff=un[i][j]*dely-un[i-1][j]*dely+vn[i][j]*delx-vn[i][j-1]*delx;
            ce[i][j]=dely*dely/a[i][j];
            cw[i][j]=dely*dely/a[i-1][j];
            cn[i][j]=delx*delx/b[i][j];
            cs[i][j]=delx*delx/b[i][j-1];
            c[i][j]=ce[i][j]+cw[i][j]+cs[i][j]+cn[i][j];
            pn[i][j]=(ce[i][j]*pn[i+1][j]+cw[i][j]*pn[i-1][j]+cn[i][j]*pn[i][j+1]+cs[i][j]*pn[i][j-1]-ff)/c[i][j];

        }
   }


    for(j=1;j<n-2;j++)  //Pressure Correction
   {
        for(i=1;i<m-2;i++)
        {
            p[i][j]=p[i][j]+pn[i][j];

        }
   }
    for(j=1;j<n-2;j++)     //Velocity correction u
   {
        for(i=1;i<m-1;i++)
        {
            u[i][j]=un[i][j]-delx*(p[i+1][j]-p[i-1][j])/b[i][j];

        }
   }
 for(j=1;j<n-1;j++)     //Velocity correction v
   {
        for(i=1;i<m-2;i++)
        {

            v[i][j]=vn[i][j]-dely*(p[i][j+1]-p[i][j-1])/a[i][j];


        }
   }
    for(j=1;j<n-1;j++)
	{
		for(i=1;i<m-1;i++)
		{
        sum1+=pow((u[i][j]-un[i][j]),2);//computing sum for omega
        sum2+=pow((v[i][j]-vn[i][j]),2); //computing sum for psi

       }
   }

    e1=sqrt(sum1/((m-2)*(n-2))); // calculating error for vorticity
    e2=sqrt(sum2/((m-2)*(n-2)));// calculating error for stream function
    k++;  //incrementing the no of iterations

}while(e1>pow(10,-6)||e2>pow(10,-6));  //comparing if the error is within the given limit for the maximum one
printf("%d\t%lf\t%lf\n",k,e1,e2);


return 0;
}
