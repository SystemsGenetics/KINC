#include "powerlaw.h"
#include "powerlaw_input.h"
#include "correlationmatrix.h"



using namespace std;
using RawPair = CorrelationMatrix::RawPair;






/*!
 * Return the total number of blocks this analytic must process as steps
 * or blocks of work.
 */
int PowerLaw::size() const
{
   EDEBUG_FUNC(this);

   return 1;
}






/*!
 * Process the given index with a possible block of results if this analytic
 * produces work blocks. This analytic implementation has no work blocks.
 *
 * @param result
 */
void PowerLaw::process(const EAbstractAnalytic::Block*)
{
   EDEBUG_FUNC(this);

   // initialize log text stream
   QTextStream stream(_logfile);

   // load raw correlation data, row-wise maximums
   std::vector<RawPair> pairs {_input->dumpRawData()};
   std::vector<float> maximums {computeMaximums(pairs)};

   // continue until network is sufficiently scale-free
   float threshold {_thresholdStart};

   while ( true )
   {
      qInfo("\n");
      qInfo("threshold: %8.3f", threshold);

      // compute adjacency matrix based on threshold
      size_t size;
      std::vector<bool> adjacencyMatrix {computeAdjacencyMatrix(pairs, maximums, threshold, &size)};

      qInfo("adjacency matrix: %lu", size);

      // make sure that adjacency matrix is not empty
      float correlation {0};

      if ( size > 0 )
      {
         // compute degree distribution of matrix
         std::vector<int> histogram {computeDegreeDistribution(adjacencyMatrix, size)};

         // compute correlation of degree distribution
         correlation = computeCorrelation(histogram);

         qInfo("correlation: %8.3f", correlation);
      }

      // output to log file
      stream << threshold << "\t" << size << "\t" << correlation << "\n";

      // TODO: break if network is sufficently scale-free

      // decrement threshold and fail if minimum threshold is reached
      threshold -= _thresholdStep;
      if ( threshold < _thresholdStop )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("Power-law Threshold Error"));
         e.setDetails(tr("Could not find scale-free network above stopping threshold."));
         throw e;
      }
   }

   // write final threshold
   stream << threshold << "\n";
}






/*!
 * Make a new input object and return its pointer.
 */
EAbstractAnalytic::Input* PowerLaw::makeInput()
{
   EDEBUG_FUNC(this);

   return new Input(this);
}






/*!
 * Initialize this analytic. This implementation checks to make sure the input
 * data object and output log file have been set, and that various integer
 * arguments are valid.
 */
void PowerLaw::initialize()
{
   EDEBUG_FUNC(this);

   // make sure input and output were set properly
   if ( !_input || !_logfile )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Did not get valid input or logfile arguments."));
      throw e;
   }

   // make sure threshold arguments are valid
   if ( _thresholdStart <= _thresholdStop )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Starting threshold must be greater than stopping threshold."));
      throw e;
   }
}






/*!
 * Compute the row-wise maximums of a correlation matrix.
 *
 * @param pairs
 */
std::vector<float> PowerLaw::computeMaximums(const std::vector<RawPair>& pairs)
{
   EDEBUG_FUNC(this,&pairs);

   // initialize elements to minimum value
   std::vector<float> maximums(_input->geneSize(), 0);

   // compute maximum correlation of each row
   for ( auto& pair : pairs )
   {
      int i = pair.index.getX();

      for ( size_t k = 0; k < pair.correlations.size(); ++k )
      {
         float correlation = fabs(pair.correlations[k]);

         if ( maximums[i] < correlation )
         {
            maximums[i] = correlation;
         }
      }
   }

   // return row-wise maximums
   return maximums;
}






/*!
 * Compute the adjacency matrix of a correlation matrix with a given threshold.
 * This function uses the pre-computed row-wise maximums for faster computation.
 * Additionally, all zero-columns removed. The number of rows in the adjacency
 * matrix is returned as a pointer argument.
 *
 * @param pairs
 * @param maximums
 * @param threshold
 * @param size
 */
std::vector<bool> PowerLaw::computeAdjacencyMatrix(const std::vector<RawPair>& pairs, const std::vector<float>& maximums, float threshold, size_t* size)
{
   EDEBUG_FUNC(this,&pairs,&maximums,threshold,size);

   // generate vector of row indices that have a correlation above threshold
   std::vector<int> indices(_input->geneSize(), -1);
   size_t pruneSize = 0;

   for ( size_t i = 0; i < maximums.size(); ++i )
   {
      if ( maximums[i] >= threshold )
      {
         indices[i] = pruneSize;
         pruneSize++;
      }
   }

   // extract adjacency matrix from correlation matrix
   std::vector<bool> adjacencyMatrix(pruneSize * pruneSize);

   // initialize diagonal
   for ( size_t i = 0; i < pruneSize; ++i )
   {
      adjacencyMatrix[i * pruneSize + i] = 1;
   }

   // iterate through all pairs
   for ( auto& pair : pairs )
   {
      // get indices into pruned matrix
      int i = indices[pair.index.getX()];
      int j = indices[pair.index.getY()];

      // skip pair if it was pruned
      if ( i == -1 || j == -1 )
      {
         continue;
      }

      // select correlation from pair
      float correlation = pair.correlations[0];

      // save correlation if it is above threshold
      if ( fabs(correlation) >= threshold )
      {
         adjacencyMatrix[i * pruneSize + j] = 1;
         adjacencyMatrix[j * pruneSize + i] = 1;
      }
   }

   // save size of adjacency matrix
   *size = pruneSize;

   // return adjacency matrix
   return adjacencyMatrix;
}






/*!
 * Compute the degree distribution of an adjacency matrix.
 *
 * @param matrix
 * @param size
 */
std::vector<int> PowerLaw::computeDegreeDistribution(const std::vector<bool>& matrix, size_t size)
{
   EDEBUG_FUNC(this,&matrix,size);

   // compute degree of each node
   std::vector<int> degrees(size);

   for ( size_t i = 0; i < size; i++ )
   {
      for ( size_t j = 0; j < size; j++ )
      {
         degrees[i] += matrix[i * size + j];
      }
   }

   // compute max degree
   int max {0};

   for ( size_t i = 0; i < degrees.size(); i++ )
   {
      if ( max < degrees[i] )
      {
         max = degrees[i];
      }
   }

   // compute histogram of degrees
   std::vector<int> histogram(max);

   for ( size_t i = 0; i < degrees.size(); i++ )
   {
      if ( degrees[i] > 0 )
      {
         histogram[degrees[i] - 1]++;
      }
   }

   return histogram;
}






/*!
 * Compare a degree distribution to a power-law distribution. The goodness-of-fit
 * is measured by the Pearson correlation of the log-transformed histogram.
 *
 * @param histogram
 */
float PowerLaw::computeCorrelation(const std::vector<int>& histogram)
{
   EDEBUG_FUNC(this,&histogram);

   // compute log-log transform of histogram data
   const int n = histogram.size();
   std::vector<float> x(n);
   std::vector<float> y(n);

   for ( int i = 0; i < n; i++ )
   {
      x[i] = log(i + 1);
      y[i] = log(histogram[i] + 1);
   }

   // visualize log-log histogram
   qInfo("histogram:");

   for ( int i = 0; i < 10; i++ )
   {
      float sum {0};
      for ( int j = i * n / 10; j < (i + 1) * n / 10; j++ )
      {
         sum += y[j];
      }

      int len {(int)(sum / log((float) _input->geneSize()))};
      QString bin(len, '#');

      qInfo(" | %s", bin.toStdString().c_str());
   }

   // compute Pearson correlation of x, y
   float sumx = 0;
   float sumy = 0;
   float sumx2 = 0;
   float sumy2 = 0;
   float sumxy = 0;

   for ( int i = 0; i < n; ++i )
   {
      sumx += x[i];
      sumy += y[i];
      sumx2 += x[i] * x[i];
      sumy2 += y[i] * y[i];
      sumxy += x[i] * y[i];
   }

   return (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
}
