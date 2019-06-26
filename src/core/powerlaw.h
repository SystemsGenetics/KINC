#ifndef POWERLAW_H
#define POWERLAW_H
#include <ace/core/core.h>
#include "correlationmatrix.h"



/*!
 * This class implements the Power-law thresholding analytic. This analytic takes
 * a correlation matrix and attempts to find a threshold which, when applied to
 * the correlation matrix, produces a scale-free network. Each thresholded network
 * is evaluated by comparing the degree distribution of the network to a power-law
 * distribution. This process is repeated at each threshold step from the starting
 * threshold.
 */
class PowerLaw : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalyticBlock* result) override final;
   virtual EAbstractAnalyticInput* makeInput() override final;
   virtual void initialize();
private:
   std::vector<float> computeMaximums(const std::vector<CorrelationMatrix::RawPair>& pairs);
   std::vector<bool> computeAdjacencyMatrix(const std::vector<CorrelationMatrix::RawPair>& pairs, const std::vector<float>& maximums, float threshold, size_t* size);
   std::vector<int> computeDegreeDistribution(const std::vector<bool>& matrix, size_t size);
   float computeCorrelation(const std::vector<int>& histogram);
   /*!
    * Pointer to the input correlation matrix.
    */
   CorrelationMatrix* _input {nullptr};
   /*!
    * Pointer to the output log file.
    */
   QFile* _logfile {nullptr};
   /*!
    * The starting threshold.
    */
   float _thresholdStart {0.99f};
   /*!
    * The threshold decrement.
    */
   float _thresholdStep {0.01f};
   /*!
    * The stopping threshold. The analytic will fail if it cannot find a
    * proper threshold before reaching the stopping threshold.
    */
   float _thresholdStop {0.5f};
};



#endif
