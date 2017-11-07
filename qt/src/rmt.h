#ifndef RMT_H
#define RMT_H
#include <ace/core/AceCore.h>



class CorrelationMatrix;



class RMT : public EAbstractAnalytic
{
public:
   enum Arguments
   {
      InputData = 0
      ,OutputFile
      ,LogFile
      ,FilterSize
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   virtual quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
private:
   float determineThreshold();
   float determineChi(float threshold, int* size);
   void generateGeneThresholds();
   QVector<double> generatePruneMatrix(float threshold, int* size);
   QVector<double> generateMatrixEigens(QVector<double>* pruneMatrix, int size);
   float getChiSquare(QVector<double>* eigens);
   float getPaceChiSquare(const QVector<double>& eigens, int pace);
   QVector<double> unfold(const QVector<double>& eigens, int pace);
   void degenerate(QVector<double>* eigens);
   constexpr static float _nnsdHistogramBin {0.05};
   constexpr static int _minUnfoldingPace {10};
   constexpr static int _maxUnfoldingPace {41};
   constexpr static int _minEigenVectorSize {100};
   constexpr static float _chiSquareBinSize {0.05};
   constexpr static double _minimumEigenValue {0.000001};
   CorrelationMatrix* _input {nullptr};
   QFile* _output {nullptr};
   QFile* _logfile {nullptr};
   int _filterSize {10};
   float _initialThreshold {0.99607};
   float _thresholdStep {0.001};
   float _thresholdMinimum {0.5};
   QVector<float> _geneThresholds;
};



#endif
