#include "rmt.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>



void RMT::swapD(double* l, int idx1, int idx2)
{
   double temp = l[idx1];
   l[idx1] = l[idx2];
   l[idx2] = temp;
}



void RMT::quickSortD(double* l, int size)
{
   if (size <= 1)
   {
      return;
   }
   int pivIdx = (int)size/1.618;
   double pivot = l[pivIdx];
   swapD(l,pivIdx,size-1);
   int leftPlace = 0;
   int i;
   for (i = 0;i<(size-1);++i)
   {
      if(l[i]<pivot)
      {
         swapD(l,i,leftPlace);
         ++leftPlace;
      }
   }
   swapD(l,size-1,leftPlace);
   quickSortD(l,leftPlace);
   quickSortD(&l[leftPlace+1],size-leftPlace-1);
}



double* RMT::unfolding(float* e, int size, int m)
{
   int count = 1;
   int i, j = 0;
   for(i = 0;i<(size-m);i+=m)
   {
      count++;
   }
   double* oX = (double*)malloc(sizeof(double)*count);
   double* oY = (double*)malloc(sizeof(double)*count);
   for(i = 0;i<(size-m);i+=m)
   {
      oX[j] = e[i];
      oY[j] = (i+1.0)/(double)size;
      j++;
   }
   oX[count-1] = e[size-1];
   oY[count-1] = 1;
   gsl_interp_accel* acc = gsl_interp_accel_alloc();
   gsl_spline* spline = gsl_spline_alloc(gsl_interp_akima,count);
   gsl_spline_init(spline,oX,oY,count);
   double* yy = (double*)malloc(sizeof(double)*size);
   for (i = 0;i<(size-2);++i)
   {
      yy[i+1] = gsl_spline_eval(spline,e[i+1],acc);
   }
   yy[0] = 0.0;
   yy[size-1] = 1.0;
   for (i = 0;i<(size-1);++i)
   {
      yy[i] = (yy[i+1]-yy[i])*size;
   }
   quickSortD(yy,size-1);
   gsl_spline_free(spline);
   gsl_interp_accel_free(acc);
   free(oX);
   free(oY);
   return yy;
}



float* RMT::degenerate(float* eigens, int size, int* newSize)
{
   int i, j = 0;
   int count = 1;
   int* flags;
   float* remDups;
   for (i = 0;i<size;++i)
   {
      if (fabs(eigens[i])<0.000001)
      {
         eigens[i] = 0.0;
      }
   }
   flags = (int*)malloc(sizeof(int)*size);
   memset(flags,0,size*sizeof(int));
   float temp = eigens[0];
   flags[0] = 1;
   for(i = 1;i<size;++i)
   {
      if(fabs(eigens[i]-temp)>0.000001)
      {
         count++;
         flags[i] = 1;
         temp = eigens[i];
      }
   }
   remDups = (float*)malloc(sizeof(float)*count);
   for(i = 0;i<size;++i)
   {
      if(flags[i]==1)
      {
         remDups[j] = eigens[i];
         j++;
      }
   }
   free(flags);
   *newSize = count;
   return remDups;
}



double RMT::getNNSDChiSquare(float* eigens, int size)
{
   double chiTest = 0;
   double avg_chiTest {0};
   int i = 0;
   int m;
   float * newE;
   int newSize;
   for (m = minUnfoldingPace;m<maxUnfoldingPace;++m)
   {
      newE = degenerate(eigens,size,&newSize);
      if (newSize<minEigenVectorSize)
      {
         continue;
      }
      if ((newSize/m)<5)
      {
         continue;
      }
      chiTest = getNNSDPaceChiSquare(newE,newSize,nnsdHistogramBin,m);
      avg_chiTest += chiTest;
      i++;
      free(newE);
   }
   return avg_chiTest/i;
}



double RMT::getNNSDPaceChiSquare(float* eigens, int size, double bin, int pace)
{
   double * edif;
   double obj;
   double expect;
   double chi = 0;
   int i, j, count;
   edif = unfolding(eigens,size,pace);
   --size;
   int n = (int)(3.0/bin)+1;
   for (i = 0;i<n;++i)
   {
      count = 0;
      for (j = 0;j<size;j++)
      {
         if (edif[j]>(i*bin)&&edif[j]<((i+1)*bin))
         {
            count++;
         }
      }
      expect = (exp(-1*i*bin)-exp(-1*(i+1)*bin))*size;
      obj = (double)count;
      chi += (obj-expect)*(obj-expect)/expect;
   }
   free(edif);
   return chi;
}
