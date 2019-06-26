#ifndef RMT_H
#define RMT_H
#include <ace/core/core.h>
#include "correlationmatrix.h"



/*!
 * This class implements the RMT analytic. This analytic takes a correlation
 * matrix and attempts to find a threshold which, when applied to the correlation
 * matrix, produces a non-random network. This analytic uses Random Matrix Theory
 * (RMT), which involves computing the eigenvalues of a thresholded correlation
 * matrix, computing the nearest-neighbor spacing distribution (NNSD) of the eigenvalues,
 * and comparing the distribution to a Poisson distribution using a chi-squared
 * test. This process is repeated at each threshold step from the starting threshold;
 * as the threshold decreases, the NNSD transitions from a Poisson distribution
 * to a Gaussian orthogonal ensemble (GOE) distribution, which causes the
 * chi-squared value to increase sharply. The final threshold is chosen as the
 * lowest threshold whose chi-squared value was below the critical value.
 */
class RMT : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalyticBlock* result) override final;
   virtual EAbstractAnalyticInput* makeInput() override final;
   virtual void initialize();
private:
   /*!
    * Defines the reduction methods this analytic supports.
    */
   enum class ReductionMethod
   {
      /*!
       * Select the first cluster
       */
      First
      /*!
       * Select the cluster with the highest absolute correlation
       */
      ,MaximumCorrelation
      /*!
       * Select the cluster with the largest sample size
       */
      ,MaximumSize
      /*!
       * Select a random cluster
       */
      ,Random
   };
private:
   std::vector<float> computeMaximums(const std::vector<CorrelationMatrix::RawPair>& pairs);
   std::vector<float> computePruneMatrix(const std::vector<CorrelationMatrix::RawPair>& pairs, const std::vector<float>& maximums, float threshold, size_t* size);
   std::vector<float> computeEigenvalues(std::vector<float>* pruneMatrix, size_t size);
   std::vector<float> computeUnique(const std::vector<float>& values);
   float computeChiSquare(const std::vector<float>& eigens);
   float computeChiSquareHelper(const std::vector<float>& values);
   std::vector<float> computeSpline(const std::vector<float>& values, int pace);
   std::vector<float> computeSpacings(const std::vector<float>& values);
   /*!
    * Pointer to the input correlation matrix.
    */
   CorrelationMatrix* _input {nullptr};
   /*!
    * Pointer to the output log file.
    */
   QFile* _logfile {nullptr};
   /*!
    * The reduction method to use. Pairwise reduction is used to select pairwise
    * correlations when there are multiple correlations per pair. By default, the
    * first cluster is selected from each pair.
    */
   ReductionMethod _reductionMethod {ReductionMethod::First};
   /*!
    * The starting threshold.
    */
   float _thresholdStart {0.99f};
   /*!
    * The threshold decrement.
    */
   float _thresholdStep {0.001f};
   /*!
    * The stopping threshold. The analytic will fail if it cannot find a
    * proper threshold before reaching the stopping threshold.
    */
   float _thresholdStop {0.5f};
   /*!
    * The critical value for the chi-squared test, which is dependent on the
    * degrees of freedom and the alpha-value of the test. This particular
    * value is based on df = 60 and alpha = 0.001. Note that since the degrees
    * of freedom corresponds to the number of histogram bins, this value
    * must be re-calculated if the number of histogram bins is changed.
    */
   float _chiSquareThreshold1 {99.607f};
   /*!
    * The final chi-squared threshold. Once the chi-squared test goes below the
    * chi-squared critical value, it must go above this value in order for the
    * analytic to find a proper threshold.
    */
   float _chiSquareThreshold2 {200};
   /**
    * The number of threads to use during eigenvalue computation.
    */
   int _numThreads {1};
   /*!
    * The minimum difference required between an eigenvalue and the previous
    * eigenvalue in ascending order for the eigenvalue to be considered unique.
    */
   float _uniqueEpsilon {1e-6f};
   /*!
    * The minimum number of unique eigenvalues which must exist in a pruned matrix
    * for the analytic to compute the NNSD of the eigenvalues. If the number of
    * unique eigenvalues is less, the chi-squared test for that threshold is skipped.
    */
   int _minUniqueEigenvalues {50};
   /*!
    * Whether to perform spline interpolation on each set of eigenvalues before
    * computing the spacings. If this option is enabled then the chi-squared value
    * for each set of eigenvalues will be the average of multiple tests in which
    * the spline pace is varied (according to the minimum and maximum spline pace);
    * otherwise, only one test is performed for each set of eigenvalues.
    */
   bool _splineInterpolation {true};
   /*!
    * The minimum pace of the spline interpolation.
    */
   int _minSplinePace {10};
   /*!
    * The maximum pace of the spline interpolation.
    */
   int _maxSplinePace {40};
   /*!
    * The number of histogram bins in the NNSD of eigenvalues. This value
    * corresponds to the degrees of freedom in the chi-squared test, therefore
    * it affects the setting of the chi-squared critical value.
    */
   int _histogramBinSize {60};
};



#endif
