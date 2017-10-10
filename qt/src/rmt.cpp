#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

#include "rmt.h"
#include "correlationmatrix.h"



using namespace std;






EAbstractAnalytic::ArgumentType RMT::getArgumentData(int argument)
{
}






QVariant RMT::getArgumentData(int argument, EAbstractAnalytic::Role role)
{
}






void RMT::setArgument(int argument, QVariant value)
{
}






void RMT::setArgument(int argument, QFile *file)
{
}






void RMT::setArgument(int argument, EAbstractData *data)
{
}






bool RMT::initialize()
{
}






void RMT::runSerial()
{
}






float RMT::determineThreshold()
{
   generateGeneThresholds();
   float chi {0.0};
   float threshold {_initialThreshold};
   QList<float> previousChi;
   QList<float> previousThresholds;
   while ( ( chi = determineChi(threshold) ) < 200.0 )
   {
      previousChi.push_back(chi);
      previousThresholds.push_back(threshold);
      threshold -= _thresholdStep;
      if ( threshold < _thresholdMinimum )
      {
         ;//ERROR!
      }
   }
   int i = previousChi.size()-1;
   while ( i > 0 && previousChi[i] > 100.0 )
   {
      --i;
   }
   if ( (i+1) < previousChi.size() )
   {
      ++i;
   }
   return previousThresholds[i];
}






float RMT::determineChi(float threshold)
{
   int size;
   unique_ptr<double> pruneMatrix {generatePruneMatrix(threshold,&size)};
   if ( size > 0 )
   {
      unique_ptr<float> eigens {generateMatrixEigens(pruneMatrix.get(),size)};
      float chi = getNNSDChiSquare(eigens.get(),size);
      if ( !isnan(chi) && !isinf(chi) )
      {
         return chi;
      }
   }
   return 0.0;
}






void RMT::generateGeneThresholds()
{
   _geneThresholds.reset(new float[_input->getGeneSize()]);
   for (int i = 0; i < _input->getGeneSize() ;++i)
   {
      _geneThresholds.get()[i] = -1.0;
   }
   int x {1};
   int y {0};
   CorrelationMatrix::Pair pair(_input);
   while ( x < _input->getGeneSize() )
   {
      pair.read(x,y);
      if ( !isnan(pair.at(0,0)) )
      {
         if ( pair.at(0,0) > _geneThresholds.get()[x] )
         {
            _geneThresholds.get()[x] = pair.at(0,0);
         }
         if ( pair.at(0,0) > _geneThresholds.get()[y] )
         {
            _geneThresholds.get()[y] = pair.at(0,0);
         }
      }
      CorrelationMatrix::increment(x,y);
   }
}






double* RMT::generatePruneMatrix(float threshold, int* size)
{
   QList<int> genes;
   for (int i = 0; i < _input->getGeneSize() ;++i)
   {
      if ( _geneThresholds.get()[i] >= threshold )
      {
         genes.push_back(i);
      }
   }
   unique_ptr<double> pruneMatrix {new double[genes.size()*genes.size()]};
   CorrelationMatrix::Pair pair(_input);
   for (int i = 0; i < genes.size() ;++i)
   {
      for (int j = 0; j < genes.size() ;++j)
      {
         int g1 {genes[i]};
         int g2 {genes[j]};
         if ( g1 == g2 )
         {
            pruneMatrix.get()[i*genes.size() + j] = 1.0;
         }
         else
         {
            if ( g2 > g1 )
            {
               swap(g1,g2);
            }
            pair.read(g1,g2);
            pruneMatrix.get()[i*genes.size() + j] = pair.at(0,0);
         }
      }
   }
   *size = genes.size();
   return pruneMatrix.release();
}






float* RMT::generateMatrixEigens(double* pruneMatrix, int size)
{
   unique_ptr<float> eigens {new float[size]};
   gsl_matrix_view m = gsl_matrix_view_array(pruneMatrix,size,size);
   gsl_vector* eval = gsl_vector_alloc(size);
   gsl_matrix* evec = gsl_matrix_alloc(size,size);
   gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
   gsl_eigen_symmv(&m.matrix,eval,evec,w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
   for (int i = 0; i < size ;i++)
   {
      eigens.get()[i] = gsl_vector_get(eval,i);
   }
   gsl_vector_free (eval);
   gsl_matrix_free (evec);
   return eigens.release();
}






// FOR ALL FUNCTIONS BELOW:
// These were all pulled from KINC version 1






void RMT::swapF(float* l, int idx1, int idx2) {
  float temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}






void RMT::quickSortF(float* l, int size){
  if (size <= 1) {
    return;
  }
  int pivIdx = (int) size / 1.618; //golden ratio
  float pivot = l[pivIdx];
  swapF(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for (i = 0; i < size - 1; i++) {
    if(l[i] < pivot){
      swapF(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapF(l, size - 1, leftPlace);
  quickSortF(l,leftPlace);
  quickSortF(&l[leftPlace + 1], size - leftPlace - 1);
  return;
}






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
   quickSortF(eigens,size);
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
