#include "similarity.h"
#include "similarity_input.h"
#include "similarity_resultblock.h"
#include "similarity_serial.h"
#include "similarity_workblock.h"
#include "similarity_opencl.h"



using namespace std;






int Similarity::size() const
{
   const qint64 totalPairs {(qint64) _input->getGeneSize() * (_input->getGeneSize() - 1) / 2};
   const qint64 WORK_BLOCK_SIZE { 32 * 1024 };

   return (totalPairs + WORK_BLOCK_SIZE - 1) / WORK_BLOCK_SIZE;
}






std::unique_ptr<EAbstractAnalytic::Block> Similarity::makeWork(int index) const
{
   const qint64 totalPairs {(qint64) _input->getGeneSize() * (_input->getGeneSize() - 1) / 2};
   const qint64 WORK_BLOCK_SIZE { 32 * 1024 };

   qint64 start {index * WORK_BLOCK_SIZE};
   qint64 size {min(totalPairs - start, WORK_BLOCK_SIZE)};

   return unique_ptr<EAbstractAnalytic::Block>(new WorkBlock(index, start, size));
}






std::unique_ptr<EAbstractAnalytic::Block> Similarity::makeWork() const
{
   return unique_ptr<EAbstractAnalytic::Block>(new WorkBlock);
}






std::unique_ptr<EAbstractAnalytic::Block> Similarity::makeResult() const
{
   return unique_ptr<EAbstractAnalytic::Block>(new ResultBlock);
}






void Similarity::process(const EAbstractAnalytic::Block* result)
{
   const ResultBlock* resultBlock {result->cast<ResultBlock>()};

   // iterate through all pairs in result block
   Pairwise::Index index {resultBlock->start()};

   for ( auto& pair : resultBlock->pairs() )
   {
      // save clusters whose correlations are within thresholds
      if ( pair.K > 1 )
      {
         CCMatrix::Pair ccmPair(_ccm);

         for ( qint8 k = 0; k < pair.K; ++k )
         {
            float corr = pair.correlations[k];

            if ( !isnan(corr) && _minCorrelation <= abs(corr) && abs(corr) <= _maxCorrelation )
            {
               ccmPair.addCluster();

               for ( int i = 0; i < _input->getSampleSize(); ++i )
               {
                  ccmPair.at(ccmPair.clusterSize() - 1, i) = (pair.labels[i] >= 0)
                     ? (k == pair.labels[i])
                     : -pair.labels[i];
               }
            }
         }

         if ( ccmPair.clusterSize() > 0 )
         {
            ccmPair.write(index);
         }
      }

      // save correlations that are within thresholds
      if ( pair.K > 0 )
      {
         CorrelationMatrix::Pair cmxPair(_cmx);

         for ( qint8 k = 0; k < pair.K; ++k )
         {
            float corr = pair.correlations[k];

            if ( !isnan(corr) && _minCorrelation <= abs(corr) && abs(corr) <= _maxCorrelation )
            {
               cmxPair.addCluster();
               cmxPair.at(cmxPair.clusterSize() - 1, 0) = corr;
            }
         }

         if ( cmxPair.clusterSize() > 0 )
         {
            cmxPair.write(index);
         }
      }

      ++index;
   }
}






EAbstractAnalytic::Input* Similarity::makeInput()
{
   return new Input(this);
}






EAbstractAnalytic::Serial* Similarity::makeSerial()
{
   return new Serial(this);
}






EAbstractAnalytic::OpenCL* Similarity::makeOpenCL()
{
   return new OpenCL(this);
}






void Similarity::initialize()
{
   // make sure input and output are valid
   if ( !_input || !_ccm || !_cmx )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input and/or output arguments."));
      throw e;
   }

   // make sure cluster range is valid
   if ( _maxClusters < _minClusters )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Minimum clusters must be less than or equal to maximum clusters."));
      throw e;
   }

   // make sure kernel size is a power of 2
   int pow2 {1};
   while ( pow2 < _kernelSize )
   {
      pow2 *= 2;
   }

   if ( _kernelSize != pow2 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Kernel size must be a power of two."));
      throw e;
   }

   // initialize cluster matrix
   _ccm->initialize(_input->getGeneNames(), _maxClusters, _input->getSampleNames());

   // initialize correlation matrix
   EMetaArray correlations;
   correlations.append(_corrModel->getName());

   _cmx->initialize(_input->getGeneNames(), _maxClusters, correlations);
}
