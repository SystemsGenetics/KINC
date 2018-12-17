#ifndef POWERLAW_H
#define POWERLAW_H
#include <ace/core/core.h>
#include "correlationmatrix.h"



/*!
 * This class implements the Power-law thresholding analytic. This analytic takes
 * a correlation matrix and attempts to find a threshold which, when applied to
 * the correlation matrix, produces a scale-free network. Each thresholded network
 * is evaluted by comparing the degree distribution of the network to a power-law
 * distribution. This process is repeated at each threshold step from the starting
 * threshold.
 */
class PowerLaw : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   QVector<float> computeMaximums(const QVector<CorrelationMatrix::RawPair>& pairs);
   QVector<bool> computeAdjacencyMatrix(const QVector<CorrelationMatrix::RawPair>& pairs, const QVector<float>& maximums, float threshold, int* size);
   QVector<int> computeDegreeDistribution(const QVector<bool>& matrix, int size);
   float computeCorrelation(const QVector<int>& histogram);
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
   float _thresholdStart {0.99};
   /*!
    * The threshold decrement.
    */
   float _thresholdStep {0.01};
   /*!
    * The stopping threshold. The analytic will fail if it cannot find a
    * proper threshold before reaching the stopping threshold.
    */
   float _thresholdStop {0.5};
};



#endif
