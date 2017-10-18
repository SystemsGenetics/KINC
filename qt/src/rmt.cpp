#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <ace/core/metadata.h>

#include "rmt.h"
#include "correlationmatrix.h"
#include "datafactory.h"
#include <iostream>



using namespace std;






EAbstractAnalytic::ArgumentType RMT::getArgumentData(int argument)
{
   // use type declaration
   using Type = EAbstractAnalytic::ArgumentType;

   // figure out which argument is being queried and return its type
   switch (argument)
   {
   case InputData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   default: return Type::Bool;
   }
}






QVariant RMT::getArgumentData(int argument, EAbstractAnalytic::Role role)
{
   // use role declaration
   using Role = EAbstractAnalytic::Role;

   // figure out which role is being queried
   switch (role)
   {
   case Role::CommandLineName:
      // figure out which argument is being queried and return command line name
      switch (argument)
      {
      case InputData: return QString("input");
      case OutputFile: return QString("output");
      default: return QVariant();
      }
   case Role::Title:
      // figure out which argument is being queried and return title
      switch (argument)
      {
      case InputData: return tr("Input:");
      case OutputFile: return tr("Output:");
      default: return QVariant();
      }
   case Role::WhatsThis:
      // figure out which argument is being queried and return "What's This?" text
      switch (argument)
      {
      case InputData: return tr("Input correlation matrix that will be used to output correlation"
                                " above a threshold determined by Random Matrix Theory.");
      case OutputFile: return tr("Output text file that will hold listing of all gene correlations"
                                 " above a certain threshold.");
      default: return QVariant();
      }
   case Role::DataType:
      // see if this is input data and return type else return nothing
      switch (argument)
      {
      case InputData: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case Role::FileFilters:
      // see if this is output file and return filter else return nothing
      switch (argument)
      {
      case OutputFile: return tr("Text File %1").arg("(*.txt)");
      default: return QVariant();
      }
   default:
      return QVariant();
   }
}






void RMT::setArgument(int argument, QFile *file)
{
   // set output argument if this is correct argument
   if ( argument == OutputFile )
   {
      _output = file;
   }
}






void RMT::setArgument(int argument, EAbstractData *data)
{
   // set input argument if this is correct argument
   if ( argument == InputData )
   {
      _input = dynamic_cast<CorrelationMatrix*>(data);
   }
}






bool RMT::initialize()
{
   // make sure input and output were set properly
   if ( !_input || !_output )
   {
      // report argument error and fail
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Argument Error"));
      e.setDetails(QObject::tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // nothing to pre-allocate
   return false;
}






void RMT::runSerial()
{
   // determine threshold using RMT and make sure interruption is not requested
   float threshold {determineThreshold()};
   if ( isInterruptionRequested() )
   {
      return;
   }

   // initialize text stream for output and get gene names
   QTextStream stream(_output);
   const EMetadata::List* names {_input->getGeneNames().toArray()};

   // initialize pair iterator, xy gene indexes, and last percent complete
   CorrelationMatrix::Pair pair(_input);
   int x {1};
   int y {0};
   int lastPercent {66};

   // iterate through all gene pairs until end is reached
   while ( x < _input->getGeneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read gene pair and check if it meets threshold
      pair.read(x,y);
      if ( pair.at(0,0) >= threshold )
      {
         // write correlation to output
         stream << *(names->at(x)->toString()) << "\t" << *(names->at(y)->toString()) << "\tco\t"
                << pair.at(0,0) << "\n";
      }

      // increment to next gene pair
      CorrelationMatrix::increment(x,y);

      // determine new percentage complete and check if new
      int newPercent {66 + 34*x/_input->getGeneSize()};
      if ( newPercent != lastPercent )
      {
         // update to new percentage and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}






float RMT::determineThreshold()
{
   // generate list of maximum threshold for each gene
   generateGeneThresholds();

   // initialize chi, threshold, and history of both
   float chi {0.0};
   float threshold {_initialThreshold};
   QList<float> previousChi;
   QList<float> previousThresholds;

   // initialize last percent, steps, and total steps
   int lastPercent {33};
   int steps {0};
   int totalSteps = (_initialThreshold - _thresholdMinimum)/_thresholdStep;

   // continue while chi is less than 200
   while ( ( chi = determineChi(threshold) ) < 200.0 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // record previous chi and threshold
      previousChi.push_back(chi);
      previousThresholds.push_back(threshold);

      // decrement threshold by step and make sure minimum is not reached
      threshold -= _thresholdStep;
      if ( threshold < _thresholdMinimum )
      {
         // report no threshold could be found and fail
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("RMT Threshold Error"));
         e.setDetails(QObject::tr("Could not find non-random threshold above minimum."));
         throw e;
      }

      // determine new percent complete and check if it is new
      int newPercent {33 + 33*steps/totalSteps};
      if ( newPercent != lastPercent )
      {
         // update new percentage and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }

      // increment steps taken
      ++steps;
   }

   // go back into history of chi values and find first one that was above 100
   int i = previousChi.size()-1;
   while ( i > 0 && previousChi[i] > 100.0 )
   {
      --i;
   }
   if ( (i+1) < previousChi.size() )
   {
      ++i;
   }

   // return threshold where chi was first above 100
   return previousThresholds[i];
}






float RMT::determineChi(float threshold)
{
   // initialize size and generate pruned matrix based off threshold
   int size {0};
   unique_ptr<double> pruneMatrix {generatePruneMatrix(threshold,&size)};

   // check and make sure matrix is not empty
   if ( size > 0 )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // generate eigen vector of pruned matrix and make sure interruption is not requested
      unique_ptr<float> eigens {generateMatrixEigens(pruneMatrix.get(),size)};
      if ( isInterruptionRequested() )
      {
         return 0.0;
      }

      // generate chi from eigen vector and return it if it is a real number
      float chi = getNNSDChiSquare(eigens.get(),size);
      if ( !isnan(chi) && !isinf(chi) )
      {
         return chi;
      }
   }

   // if no real chi was found return zero
   return 0.0;
}






void RMT::generateGeneThresholds()
{
   // allocate memory for thresholds matrix and initialize all to the minimum
   _geneThresholds.reset(new float[_input->getGeneSize()]);
   for (int i = 0; i < _input->getGeneSize() ;++i)
   {
      _geneThresholds.get()[i] = -1.0;
   }

   // initialize percent complete and steps
   int lastPercent {0};
   qint64 steps {0};
   qint64 totalSteps {_input->getGeneSize()*(_input->getGeneSize() - 1)/2};

   // xy gene indexes and iterator
   int x {1};
   int y {0};
   CorrelationMatrix::Pair pair(_input);

   // iterate through all gene pairs
   while ( x < _input->getGeneSize() )
   {
      // make sure interruption is not requested
      if ( isInterruptionRequested() )
      {
         return;
      }

      // read in gene pair and check if it is a real number
      pair.read(x,y);
      if ( !isnan(pair.at(0,0)) )
      {
         // if value is greater than current max of gene x set it to new max
         if ( pair.at(0,0) > _geneThresholds.get()[x] )
         {
            _geneThresholds.get()[x] = pair.at(0,0);
         }

         // if value is greater than current max of gene y set it to new max
         if ( pair.at(0,0) > _geneThresholds.get()[y] )
         {
            _geneThresholds.get()[y] = pair.at(0,0);
         }
      }

      // increment to next gene pair, steps, and compute new percent complete
      CorrelationMatrix::increment(x,y);
      ++steps;
      qint64 newPercent {33*steps/totalSteps};

      // check to see if percent complete changed
      if ( newPercent != lastPercent )
      {cout << newPercent << "\n";
         // update percent complete and emit progressed signal
         lastPercent = newPercent;
         emit progressed(lastPercent);
      }
   }
}






double* RMT::generatePruneMatrix(float threshold, int* size)
{
   // generate vector of gene indexes that have max threshold above given threshold
   QVector<int> genes;
   for (int i = 0; i < _input->getGeneSize() ;++i)
   {
      if ( _geneThresholds.get()[i] >= threshold )
      {
         genes.push_back(i);
      }
   }

   // allocate new pruned matrix with gene vector size and initialize iterator
   unique_ptr<double> pruneMatrix {new double[genes.size()*genes.size()]};
   CorrelationMatrix::Pair pair(_input);

   // populate the pruned matrix with gene correlation values
   for (int i = 0; i < genes.size() ;++i)
   {
      for (int j = 0; j < genes.size() ;++j)
      {
         // make sure interruption is not requested
         if ( isInterruptionRequested() )
         {
            return nullptr;
         }

         // get indexes for both genes
         int g1 {genes[i]};
         int g2 {genes[j]};

         // if genes have same index set ij matrix value to 1
         if ( g1 == g2 )
         {
            pruneMatrix.get()[i*genes.size() + j] = 1.0;
         }

         // else get correlation value
         else
         {
            // make sure first gene index is bigger
            if ( g2 > g1 )
            {
               swap(g1,g2);
            }

            // read gene pair and set correlation value to pruned matrix
            pair.read(g1,g2);
            pruneMatrix.get()[i*genes.size() + j] = pair.at(0,0);
         }
      }
   }

   // set size to size of pruned matrix and return pointer
   *size = genes.size();
   return pruneMatrix.release();
}






float* RMT::generateMatrixEigens(double* pruneMatrix, int size)
{
   // allocate new array to hold eigen values
   unique_ptr<float> eigens {new float[size]};

   // use GSL library for compute eigen values for pruned matrix
   gsl_matrix_view m = gsl_matrix_view_array(pruneMatrix,size,size);
   gsl_vector* eval = gsl_vector_alloc(size);
   gsl_matrix* evec = gsl_matrix_alloc(size,size);
   gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
   gsl_eigen_symmv(&m.matrix,eval,evec,w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

   // grab eigen values from GSL and set to eigen array
   for (int i = 0; i < size ;i++)
   {
      eigens.get()[i] = gsl_vector_get(eval,i);
   }

   // cleanup rest of GSL resoruces and return eigen array
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
