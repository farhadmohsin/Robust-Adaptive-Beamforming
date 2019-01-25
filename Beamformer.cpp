/*
This code takes as input the audio from 8 channel.
The audio clips are taken as a structure, that consists of three elements.
int samples - the number of samples in each audio channel
int Fs - the sampling frequency
double sig[] - the audio samples

The maximum sample number for this code is 80000 which can be changed by changing the #define Length pre-processor line
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Length 80006
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

// Audio Structure
struct Audio
{
	int samples, Fs;
	double sig[Length];
}sigArray[8];

void spline(double y[],double yy[]);
double bessi0( double x );
double bessi1( double x );
void FFT(int dir, long m, double xr[], double xi[]);

double cof_spline[M*M*LG];

int main( int argc, char **argv )
{
	int i,j,k;
	int m;
	int win_count = 1;
	FILE *f_in,*f_out;
	if((f_in = fopen("sigArray.txt", "r")) == NULL)
	{
		printf("Error reading input file");
		return 0;
	}
	for(i=0;i<8;i++)
	{
		fscanf(f_in,"%d%d",&sigArray[i].Fs,&sigArray[i].samples);
		for(j=0;j<sigArray[i].samples;j++)
			fscanf(f_in,"%lf",&sigArray[i].sig[j]);
	}
	fclose(f_in);

//Array Init

	double c = 340;
	double di = 0.05;
	double fs = 16000;

    double DOI;

//Noise spectrum estimation
	double Noise[win_by2_plus1] = {0}, Lambda[win_by2_plus1] = {0};
	double Yabs;
	//Hamming window
	double wnd[win_length];
	for(j=0;j<win_length/2;j++)
	{
	    wnd[j] = 0.54 - 0.46*cos(2*PI*((double)j)/(win_length-1.0));
	    wnd[win_length-j-1] = wnd[j];
	}
	//Noise Estimation from 0.25 second
	/*The following portion estimates Noise for the MMSE algorithm
	from the silent period at the start
	Which we have taken to consist of 4000 samples*/
	int cnt = 0,win_start;
	double Yr[win_length],Yi[win_length];
    for(win_start = 0;win_start< fs/4-win_length;win_start += win_inc,cnt++)
	{
		for(j=0;j<win_length;j++)
		{
			Yr[j] = sigArray[0].sig[win_start+j]*wnd[j];
			Yi[j] = 0;
		}
		FFT(1,9,Yr,Yi);
		for(j=0;j<win_by2_plus1;j++)
		{
			Yabs = Yr[j]*Yr[j] + Yi[j]*Yi[j];
			Noise[j] += sqrt(Yabs);
			Lambda[j] += Yabs;
		}
	}
	for(j=0;j<win_by2_plus1;j++)
	{
		Noise[j] /= cnt;
		Lambda[j] /= cnt;
	}
//DONE

	double alpha = 0.99, gamma1p5 = 0.886226925452758;
	double Gamma[win_by2_plus1],G[win_by2_plus1];
	for(i=0;i<win_by2_plus1;i++)
		Gamma[i] = G[i] = 1;

//Beamforming Init

	double delta = 1e-3;
	double maxdelay = 14;
	double u = .1;

    /*The F and P matrix formed here are to be used in the Frost Beamforming portion
    They are constant for all the windows
    */
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
	//initialize WW
	/*WW will be the vector that contains the adaptive filter waits
	J taps for M channels each gives a total of M*J taps
	*/
	double WW[MJ] = {0};
	//initalize YY
	double YY[Length];
//Start of Iterations

	double sigout[win_length][M],prev_sigout[win_length][M];
	double SIGOUTR[fft_length*8]={0},SIGOUTI[fft_length*8]={0};
	double Yphase;
	double gammanew,xi,nu;
	double ctr[fft_length],cti[fft_length];

    int cnt_60=0,cnt_90=0,cnt_110=0,cnt_150=0;
    int window_cnt = 1;

    /*Each win_start starts processing for one window. We ignore the windows that
    constitues the first 4000 samples as they were previously considered as
    silent zone and hence do not give us speech data
    We assign 0 in the output data for those 4000 samples*/
    /*
    Each window contains of 512 samples, the win_start is incremented by 307 samples
    These numbers, 512 and 307 are chosen in order to find the best balance between
    efficiency and accuracy
    */
	for(win_start=4000;win_start<=(sigArray[0].samples-win_length);win_start+=win_inc)
	{

//MMSE Part

		double avg[8],sigma[8];
		for(i=0;i<8;i++)
		{
			for(j=0;j<win_length;j++)
			{
				sigout[j][i] = sigArray[i].sig[win_start+j];
				Yr[j] = sigout[j][i]*wnd[j];
				Yi[j] = 0;
			}

			FFT(1,9,Yr,Yi);
			double Specr[win_length], Speci[win_length];
			for(j=0;j<win_by2_plus1;j++)
			{
			    double mul;
				Yabs = Yr[j]*Yr[j] + Yi[j]*Yi[j];
				Yphase = atan2(Yi[j],Yr[j]);
				gammanew = Yabs / Lambda[j];
				if(gammanew>1)   mul = gammanew - 1;
				else mul = 0;
				xi = alpha * (G[j]*G[j]) * Gamma[j] + (1-alpha)*mul;
				Gamma[j] = gammanew;
				nu = Gamma[j] * xi / (1+xi);
				if(nu/2 > 708)
					G[j] = xi / (1+xi);
				else
				{
					G[j] = gamma1p5*sqrt(nu)/Gamma[j]*exp(-nu/2.0);
					G[j] *= ((1+nu)*bessi0(nu/2) + nu*bessi1(nu/2));
				}

				Specr[j] = G[j] * sqrt(Yabs) * cos(Yphase);
				Speci[j] = G[j] * sqrt(Yabs) * sin(Yphase);
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

//DOA Estimation by MCCC Part

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

                spline(ctr+win_length-18,cof_spline+i*M*LG+j*LG);
			}
		}
		//Evaluating MCCC matrix
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
					for(double l=-18.0;l<=18.0;k++,l+=0.125)
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
        printf("%lf %lf\n",DOI*180/PI,mini_lag);

//Frost Beamforming Part

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
/*The Following part contains the code for the matrix multiplications
in the CLMS algorithm for Frost Beamforming
*/
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

		if(win_start == 0)
		{
			for(i=win_start,j=0;i<win_start+win_length;i++,j++)
				YY[i] = out[j];
		}
		else
		{
			for(i=win_start,j=0;i<win_start+win_length;i++,j++)
				YY[i] += out[j];
		}
		/*the array out[] contains the processes audio speech
		for this window
		*/
	}
//Write output
	if ((f_out = fopen("Output.txt","w")) == NULL)
	{
		printf("Error! Couldn't write file");
		return 0;
	}
	fprintf(f_out,"%d %d\n",16000,win_start);
	for (i=0;i<win_start;i++)
		fprintf(f_out,"%lf\n", YY[i]);
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
    for(double lag = -18;lag<=18;lag+=0.125,C++)
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
