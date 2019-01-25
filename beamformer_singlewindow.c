#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <xtensa/tie/xt_hifi2.h>
#include <xtensa/config/defs.h>
#include <xtensa/sim.h>

#define Length 1024
#define J 40 //Filter tap length
#define M 8 //Number of input channels
#define MJ 320
#define win_length 512
#define fft_length 1024
#define win_by2_plus1 257
#define win_inc 307
#define Lgs 37
#define LG 289
#define PI 3.1415926535
#define EPS 0.00001

#define abss(x) ((x>0)?(x):(-x))

void spline(double y[],double yy[]);
double bessi0( double x );
double bessi1( double x );
void FFT(int dir, long m, double xr[], double xi[]);

double cof_spline[M*M*LG];

int main( int argc, char **argv )
{
//Array Init
	int i,j,k;
	int m;
	int win_count = 1;
	double Noise[win_by2_plus1], Lambda[win_by2_plus1];

    double sigout[win_length][M],prev_sigout[win_length][M] = {0};
	double SIGOUTR[fft_length*8]={0},SIGOUTI[fft_length*8]={0};
	double Ysqr[win_by2_plus1],Yphase[win_by2_plus1];
	double gammanew[win_by2_plus1],xi[win_by2_plus1],nu[win_by2_plus1];
	double ctr[fft_length],cti[fft_length];

    //initialize WW
	double WW[MJ] = {0};
	//initalize YY
	double YY[Length];

	// MAKING SURE THAT SIMULATOR IS FAST FUNCTIONAL
	xt_iss_client_command("profile", "disable");
	xt_iss_switch_mode(XT_ISS_FUNCTIONAL);

	FILE *f_in,*f_out;
	if((f_in = fopen("noise_in.txt", "r")) == NULL)
	{
		printf("Error reading input file");
		return 0;
	}
	for(i=0;i<win_by2_plus1;i++)
		fscanf(f_in,"%lf",&Noise[i]);
    for(i=0;i<win_by2_plus1;i++)
		fscanf(f_in,"%lf",&Lambda[i]);
    fclose(f_in);

    if((f_in = fopen("window_in.txt", "r")) == NULL)
	{
		printf("Error reading input file");
		return 0;
	}
    for(j=0;j<M;j++)
        for(i=0;i<win_length;i++)
            fscanf(f_in,"%lf",&sigout[i][j]);

	if((f_in = fopen("window_prev.txt", "r")) == NULL)
	{
		printf("Error reading input file");
		return 0;
	}
	for(j=0;j<M;j++)
		for(i=0;i<win_inc;i++)
			fscanf(f_in,"%lf",&prev_sigout[i][j]);
	for(i=0;i<win_length-win_inc;i++)
		fscanf(f_in,"%lf",&YY[i]);
	fclose(f_in);

	// ENABLING CYCLE ACCURATE SIMULATION FROM HERE
	xt_iss_switch_mode(XT_ISS_CYCLE_ACCURATE);
	xt_iss_client_command("profile", "enable");

	double c = 340;
	double di = 0.05;
	double fs = 16000;

    double DOI;

	//Hamming window
	double wnd[win_length];
	for(j=0;j<win_length/2;j++)
	{
	    wnd[j] = 0.54 - 0.46*cos(2*PI*((double)j)/(win_length-1.0));
	    wnd[win_length-j-1] = wnd[j];
	}
	int cnt = 0,win_start;
	double Yr[win_length],Yi[win_length];

//DONE

	double alpha = 0.99, gamma1p5 = 0.886226925452758;
	double Gamma[win_by2_plus1],G[win_by2_plus1];
	for(i=0;i<win_by2_plus1;i++)
		Gamma[i] = G[i] = 1;

//Beamforming Init

	double delta = 1e-3;
	double maxdelay = 14;
	double u = .1;

	//form F
	double F[MJ] = {0};
	for(i=0;i<M;i++)
		F[i] = 0.125;

	//form P
	double P[MJ][MJ] = {0};
	for(i=0;i<MJ;i++)
		P[i][i] = 1;
	for(i=0;i<MJ;i+=8)
		for(j=i;j<i+8;j++)
			for(k=i;k<i+8;k++)
				P[j][k] -= 0.125;

//Start of Iterations

//MMSE

    double avg[8],sigma[8];
    for(i=0;i<8;i++)
    {
        for(j=0;j<win_length;j++)
        {
            Yr[j] = sigout[j][i]*wnd[j];
            Yi[j] = 0;
        }

        FFT(1,9,Yr,Yi);
        double Specr[win_length], Speci[win_length];
        for(j=0;j<win_by2_plus1;j++)
        {
            double mul;
            Ysqr[j] = Yr[j]*Yr[j] + Yi[j]*Yi[j];
            Yphase[j] = atan2(Yi[j],Yr[j]);
            gammanew[j] = Ysqr[j] / Lambda[j];
            if(gammanew[j]>1)   mul = gammanew[j] - 1;
            else mul = 0;
            xi[j] = alpha * (G[j]*G[j]) * Gamma[j] + (1-alpha)*mul;
            Gamma[j] = gammanew[j];
            nu[j] = Gamma[j] * xi[j] / (1+xi[j]);
            if(nu[j]/2 > 708)
                G[j] = xi[j] / (1+xi[j]);
            else
            {
                G[j] = gamma1p5*sqrt(nu[j])/Gamma[j]*exp(-nu[j]/2.0);
                G[j] *= ((1+nu[j])*bessi0(nu[j]/2) + nu[j]*bessi1(nu[j]/2));
            }

            Specr[j] = G[j] * sqrt(Ysqr[j]) * cos(Yphase[j]);
            Speci[j] = G[j] * sqrt(Ysqr[j]) * sin(Yphase[j]);
        }
        for(j=1;j<win_by2_plus1-1;j++)
        {
            Specr[win_length-j] = Specr[j];
            Speci[win_length-j] = -Speci[j];
        }
        FFT(-1,9,Specr,Speci);
        for(j=0;j<win_length;j++)
            sigout[j][i] = prev_sigout[j][i] + Specr[j];
        for(j=0;j<win_inc;j++)
            prev_sigout[j][i] = sigout[win_length-win_inc+j][i];
    }
//DONE

//DOA Estimation

    //Finding CCF
    for(i=0;i<M;i++)
    {
        avg[i] = sigma[i] = 0;
        for(j=0;j<win_length;j++)
            avg[i] += sigout[j][i];
        avg[i] /= win_length;
        for(j=0;j<win_length;j++)
            sigma[i] += ((sigout[j][i]-avg[i])*(sigout[j][i]-avg[i]));
        sigma[i] /= (win_length-1);
        sigma[i] = sqrt(sigma[i]);
    }
//DONE

    for(i=0;i<M;i++)
    {
        for(j=0;j<win_length;j++)
        {
            SIGOUTR[i*fft_length+j] = sigout[j][i] - avg[i];
            SIGOUTI[i*fft_length+j] = 0;
        }
        for(j=win_length;j<fft_length;j++)
        {
            SIGOUTR[i*fft_length+j] = 0;
            SIGOUTI[i*fft_length+j] = 0;
        }
        FFT(1,10,SIGOUTR+i*fft_length,SIGOUTI+i*fft_length);
    }
//DONE
    for(i=0;i<M;i++)
    {
        for(j=i+1;j<M;j++)
        {
            for(k=0;k<fft_length;k++)
            {
                ctr[k] = SIGOUTR[j*fft_length+k]*SIGOUTR[i*fft_length+k]+SIGOUTI[j*fft_length+k]*SIGOUTI[i*fft_length+k];
                cti[k] = SIGOUTI[j*fft_length+k]*SIGOUTR[i*fft_length+k]-SIGOUTR[j*fft_length+k]*SIGOUTI[i*fft_length+k];
            }
            FFT(-1,10,ctr,cti);
            for(k=0;k<fft_length;k++)
            {
                ctr[k] /= win_length;
                cti[k] = ctr[k]; //cti acts as a temp from hereon
            }
            for(k=win_length;k<fft_length;k++)
                ctr[k-win_length] = cti[k];
            for(k=0;k<win_length;k++)
                ctr[k+win_length] = cti[k]; //ctr now holds cov(i,j)
//DONE
            spline(ctr+win_length-18,cof_spline+i*M*LG+j*LG);
        }
    }
    //Evaluating MCC
    double lg,R[M][M],lag;
    double aa,bb,cc,dd,xx;
    double det,mini_det,mini_lag;
    for(lg=-2.5,m=0;lg<=2.5;lg+=0.125,m++)
    {
        //building R
        for(i=0;i<M;i++)
        {
            R[i][i] = 1;
            for(j=i+1;j<M;j++)
            {
                lag = lg*(j-i);
                k=0;
                double l;
                for(l=-18.0;l<=18.0;k++,l+=0.125)
                    if(l==lag) break;
                R[i][j] = cof_spline[i*M*LG + j*LG + k];
                R[i][j] /= (sigma[i]*sigma[j]);
                R[j][i] = R[i][j];
            }
        }

        //determinant(R)
        double L[M][M] = {0};
        for(i=0;i<M;i++)
            L[i][i] = 1;
        for(j=0;j<M;j++)
        {
            for(i=j+1;i<M;i++)
            {
                L[i][j] = R[i][j]/R[j][j];
                for(k=0;k<M;k++)
                    R[i][k] -= R[j][k]*L[i][j];
            }
        }
        det = 1;
        for(i=0;i<M;i++)
            det *= R[i][i];
        //find minimum det lag
        if(m==0)
        {
            mini_lag = lg;
            mini_det = abss(det);
        }
        else
        {
            if(abss(det) < mini_det)
            {
                mini_det = abss(det);
                mini_lag = lg;
            }
        }
    }
    //find DOI
    double coss = (c*mini_lag)/(di*16000);
    if(coss > 1)
        DOI = 7.5 * PI / 180.0;
    else if(coss < -1)
        DOI = (172.5) * PI / 180.0;
    else
        DOI = acos(coss);
    //printf("%lf %lf\n",DOI*180/PI,mini_lag);

//DONE

//Frost Beamforming

    double sample_delay[8];
    double out[win_length] = {0}; //Final output
    if(DOI <= PI/2.0)
        for(i=0;i<M;i++)
            sample_delay[i] = (di*cos(DOI)*fs/c*(7-i));
    else
        for(i=0;i<M;i++)
            sample_delay[i] = abss(di*cos(DOI)*fs/c*(i));
    // building XX
    double XX[M][J] = {0};
    int N,n;

    double h[J],temp;
    for(n=0;n<win_length;n++)
    {
        for(i=0;i<M;i++)
        {
            if(sample_delay[i] == 0)
                XX[i][0] = sigout[n][i];
            else
            {
            //fractional delay
                for(N=0;N<sample_delay[i];N++)
                    ;
                N *= 2;
                for(j=0;j<N+1;j++)
                    h[j] = 1;
                for(k=0;k<N+1;k++)
                    for(j=0;j<N+1;j++)
                        if(j!=k)
                        {
                            h[j] *= (sample_delay[i]-k)/(j-k);
                        }
                XX[i][0] = 0;
                if(N < j)
                    for(j=N;j>=0;j--)
                    {
                        temp = ((n-N+j)>=0?sigout[n-N+j][i]:prev_sigout[win_inc+n-N+j][i]);
                        XX[i][0] += h[N-j]*temp;
                    }
                else
                    for(j=0;j<N+1;j++)
                    {
                        XX[i][0] += h[N-j]*sigout[n-N+j][i];
                    }
            }
        }
//DONE
        // x = XX(:)
        // out[n] = x' * WW
        out[n] = 0;
        for(i=0;i<M;i++)
            for(j=0;j<J;j++)
            {
                out[n] += XX[i][j] * WW[i*J+j];
            }

        //WW = F+P*(WW-u*y(n)*x -u*delta*WW)
        //1. WW = WW-u*y(n)*x -u*delta*WW
        for(i=0;i<M;i++)
            for(j=0;j<J;j++)
            {
                WW[i*J+j] *= (1-u*delta);
                WW[i*J+j] -= u*out[n]*XX[i][j];
            }
        //2. WW = F + P * WW
        double tempW[MJ];
        for(i=0;i<MJ;i++)
            tempW[i] = WW[i];
        for(i=0;i<MJ;i++)
        {
            WW[i] = F[i];
            for(j=0;j<MJ;j++)
            {
                WW[i] += P[i][j] * tempW[j];
            }
        }
        //XX =[zeros(M,1) XX(:,2:end)];
        for(j=1;j<J;j++)
            for(i=0;i<M;i++)
                XX[i][j] = XX[i][j-1];
        for(i=0;i<M;i++)
            XX[i][0] = 0;
    }

//Output
    for(j=0;j<win_length;j++)
        YY[j] += out[j];

//Write output
    // MAKING SURE THAT SIMULATOR IS FAST FUNCTIONAL
	xt_iss_client_command("profile", "disable");
	xt_iss_switch_mode(XT_ISS_FUNCTIONAL);

	if ((f_out = fopen("window_out.txt","w")) == NULL)
	{
		printf("Error! Couldn't write file");
		return 0;
	}
	for (i=0;i<win_length;i++)
		fprintf(f_out,"%lf\n", YY[i]);
	fclose(f_out);

	if ((f_out = fopen("window_prev.txt","w")) == NULL)
	{
		printf("Error! Couldn't write file");
		return 0;
	}
    for(j=0;j<M;j++)
        for(i=0;i<win_inc;i++)
            fprintf(f_in,"%lf",prev_sigout[i][j]);
    for(j=0;j<M;j++)
        for(i=0;i<win_length;i++)
            fprintf(f_in,"%lf",sigout[i][j]);
    for(i=0;i<win_length-win_inc;i++)
        fprintf(f_in,"%lf",YY[i+win_inc]);
	fclose(f_out);

	return 0;
}

void spline(double y[],double yy[])
{
    int i;
    double a[Lgs],b[Lgs],c[Lgs],D[Lgs],m[Lgs],d0,zn;
    double div;
    for(i=0;i<Lgs;i++)
    {
        if(i==0 || i==(Lgs-1)) b[i] = 1;
        else b[i] = 4;
        if(i==(Lgs-1)) a[i] = -2;
        else a[i] = 1;
        if(!i) c[i] = -2;
        else c[i] = 1;
        if(i==0 || i==(Lgs-1))
            D[i] = 0;
        else
            D[i] = 6*(y[i+1]-2*y[i]+y[i-1]);
    }
    d0 = 1;
    zn = 1;
    for(i=1;i<Lgs;i++)
    {
        if(i == Lgs-2)
        {
            b[i] -= a[i]/b[i-1] * c[i-1];
            a[i+1]  -= zn/c[i-1] * c[i-1];
            D[i] -= a[i]/b[i-1] * D[i-1];
        }
        else
        {
            b[i] -= a[i]/b[i-1] * c[i-1];
            D[i] -= a[i]/b[i-1] * D[i-1];
        }
    }
    m[Lgs-1] = D[Lgs-1]/b[Lgs-1];
    for(i=Lgs-2;i>0;i--)
        m[i] = (D[i] - c[i]*m[i+1])/b[i];
    m[0] = (D[0]-c[0]*m[1]-d0*m[2])/b[0];
    int C = 0;
    double L,t,s0=0,s1=0,s2=0,s3=0;
    double lag;
    for(lag = -18;lag<=18;lag+=0.125,C++)
    {
        if(!(C%8))
        {
            i = C/8;
            yy[C] = y[i];
            L = C/8 - 18;
            if(L==18)
                continue;
            s1 = (y[i+1]-y[i]) - (2*m[i]+m[i+1])/6;
            s2 = m[i]/2;
            s3 = (m[i+1]-m[i])/6;

            continue;
        }
        s0 = y[i];
        t = lag - L;

        yy[C] = s0 + t*(s1 + t*(s2+s3*t));
    }
    return;
}


double bessi0( double x )
{
   double ax,ans;
   double y;
   if ((ax=fabs(x)) < 3.75)
   {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   }
   else
   {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

double bessi1( double x )
{
   double ax,ans;
   double y;
   if ((ax=fabs(x)) < 3.75)
   {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   }
   else
   {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}

void FFT(int dir, long m, double xr[], double xi[])
{
   long i, i1, i2,j, k, l, l1, l2, n;
   double txr, txi, t1r, t1i, ur, ui, cr, ci,temp,z;
   n = 1;
   for(i = 0; i < m; i++)
      n <<= 1;
   i2 = n >> 1;
   j = 0;
   for (i = 0; i < n-1 ; i++)
   {
      if (i < j)
      {
          temp = xr[i];
          xr[i] = xr[j];
          xr[j] = temp;
          temp = xi[i];
          xi[i] = xi[j];
          xi[j] = temp;
      }
      k = i2;
      while (k <= j)
	  {
         j -= k;
         k >>= 1;
      }
      j += k;
   }
   cr = -1.0;
   ci = 0.0;
   l2 = 1;
   for (l = 0; l < m; l++)
   {
      l1 = l2;
      l2 <<= 1;
      ur = 1.0;
      ui = 0.0;
      for (j = 0; j < l1; j++)
	  {
         for (i = j; i < n; i += l2)
		 {
            i1 = i + l1;
            t1r = ur*xr[i1] - ui*xi[i1];
            t1i = ur*xi[i1] + ui*xr[i1];
            xr[i1] = xr[i] - t1r;
            xi[i1] = xi[i] - t1i;
            xr[i] += t1r;
            xi[i] += t1i;
         }
         z = ur*cr - ui*ci;
         ui = ur*ci + ui*cr;
         ur = z;
      }
      ci = (sqrt((1.0 - cr) / 2.0));
      if (dir == 1)
         ci *= -1;
      cr = (sqrt((1.0 + cr) / 2.0));
   }
   if (dir == -1)
   {
      for (i = 0; i < n; i++)
      {
          xr[i] /= n;
          xi[i] /= n;
      }
   }
   return;
}
