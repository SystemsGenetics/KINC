#include "rmt.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>



void RMT::input(Ace::Data* in)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyInputs>(!_in,f,__LINE__);
   Ace::assert<InvalidDataType>(in->type()==std::string("cmx"),f,__LINE__);
   _in = dynamic_cast<CMatrix*>(in);
}



void RMT::output(Ace::Data*)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<TooManyOutputs>(false,f,__LINE__);
}



void RMT::execute_cl(Ace::GetOpts& ops, Ace::Terminal& tm)
{
   tm << "OpenCL accelerated not implemented, routing to serial processing...\n";
   execute_pn(ops,tm);
}



void RMT::execute_pn(Ace::GetOpts& ops, Ace::Terminal& tm)
{
   static const char* f = __PRETTY_FUNCTION__;
   Ace::assert<NoDataInput>(_in,f,__LINE__);
   // Initialize all values and prepare vectors that will hold the history of chi and threshold
   // values.
   bool pass {true};
   float chi {0.0};
   float thresh {_chiStep};
   std::vector<float> pChi;
   std::vector<float> pThresh;
   // Keep getting new matrices and get its chi value from the current threshold, keep doing this,
   // reducing the threshold by an increment each time until the chi value is 200 or greater. Stop
   // this if the threshold is reduced below 50.
   while ( chi < 200.0 )
   {
      int size;
      std::unique_ptr<double> pMatrix {prune_matrix(thresh,size)};
      std::unique_ptr<float> eigens {matrix_eigens(pMatrix,size)};
      chi = getNNSDChiSquare(eigens.get(),size);
      pChi.push_back(chi);
      pThresh.push_back(thresh);
      thresh -= _chiStep;
      if ( thresh < _chiMinimum )
      {
         tm << "Error: Could not find acceptable threshold above " << _chiMinimum << ".\n";
         pass = false;
         break;
      }
   }
   if ( pass )
   {
      // Go back in the history and find the last threshold that had a chi value above 100. This is
      // done because there can be valleys or bumps while iterating chi values.
      int i = pChi.size()-1;
      while ( i > 0 && pChi[i] > 100.0 )
      {
         i--;
      }
      ++i;
      tm << "Found threshold! " << pThresh[i] << "\n";
   }
}



/// @brief Get eigenvalues of a matrix. The matrix must be square and symmetric.
///
/// @param matrix The matrix to get eigenvalues from.
/// @param size The size, n, of the n by n matrix.
/// @return New list of eigenvalues.
std::unique_ptr<float> RMT::matrix_eigens(const std::unique_ptr<double>& matrix, int size)
{
   // Use GSL to find the eigenvalues of the given matrix and return a new list of those
   // eigenvalues.
   std::unique_ptr<float> eigens {new float[size]};
   gsl_matrix_view m = gsl_matrix_view_array(matrix.get(),size,size);
   gsl_vector *eval = gsl_vector_alloc (size);
   gsl_matrix *evec = gsl_matrix_alloc (size,size);
   gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
   gsl_eigen_symmv(&m.matrix, eval, evec, w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
   for (int i = 0; i < size ;i++)
   {
      eigens.get()[i] = gsl_vector_get(eval,i);
   }
   gsl_vector_free (eval);
   gsl_matrix_free (evec);
   return eigens;
}



/// @brief Generate reduced gene matrix through threshold pruning.
///
/// Build a correlation matrix from the input matrix, but cut out any gene that does not have any
/// matches with another gene equalling or above the given threshold.
///
/// @param threshold The threshold a gene must have at least one correlation above or equal to.
/// @param size The number of genes in the new matrix.
/// @return A new correlation matrix only containing genes that meet the threshold.
std::unique_ptr<double> RMT::prune_matrix(float threshold, int& size)
{
   // Build list of all genes that have at least one match with another gene that meets the
   // threshold.
   std::vector<int> genes;
   for (int i = 0; i < _in->gene_size() ;++i)
   {
      if ( gene_has_matches(i,threshold) )
      {
         genes.push_back(i);
      }
   }
   // Create a row ordered n by n matrix where n is the number of genes found in previous search.
   // Populate the matrix with the correlations between the genes.
   std::unique_ptr<double> pMatrix {new double[genes.size()*genes.size()]};
   for (int i = 0; i < genes.size() ;++i)
   {
      for (int j = 0; j < genes.size() ;++j)
      {
         int g1 {genes[i]};
         int g2 {genes[j]};
         if ( g1 == g2 )
         {
            pMatrix.get()[g1*genes.size()+g2] = 1.0;
         }
         else
         {
            pMatrix.get()[g1*genes.size()+g2] = _in->find(g1,g2).corrs().find(0)[0];
         }
      }
   }
   // Return the size of the matrix n and the matrix itself.
   size = genes.size();
   return pMatrix;
}



/// @brief Does gene have given correlation threshold?
///
/// Determine if a gene in the input correlation matrix has at least one correlation that meets a
/// given threshold.
///
/// @param gene The indexed gene number to investigate.
/// @param threshold The threshold value the gene must meet.
/// @return Does the gene meet the given threshold?
bool RMT::gene_has_matches(int gene, float threshold)
{
   // Go through all possible gene combinations until you find a match is meets the threshold. If
   // a match is found, stop searching and return true. If no match is found return false.
   bool ret {false};
   int i {0};
   while ( i < _in->gene_size() )
   {
      if ( gene != i && _in->find(gene,i).corrs().find(0)[0] >= threshold )
      {
         ret = true;
         break;
      }
      ++i;
   }
   return ret;
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
