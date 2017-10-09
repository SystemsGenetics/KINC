#include "rmt.h"



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
   unique_ptr<float> pruneMatrix {generatePruneMatrix(threshold,&size)};
   if ( size > 0 )
   {
      unique_ptr<float> eigens {generateMatrixEigens(pruneMatrix.get(),&size)};
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
}






float *RMT::generatePruneMatrix(float threshold, int *size)
{
}






float *RMT::generateMatrixEigens(float *pruneMatrix, int *size)
{
}






float RMT::getNNSDChiSquare(float *eigens, int size)
{
}
