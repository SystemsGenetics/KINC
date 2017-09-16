#ifndef SPEARMAN_H
#define SPEARMAN_H
#include <ace/core/AceCore.h>
#include <ace/core/AceOpenCL.h>



class ExpressionMatrix;
class CorrelationMatrix;



class Spearman : public EAbstractAnalytic
{
   Q_OBJECT
public:
   ~Spearman();
   enum Arguments
   {
      InputData = 0
      ,OutputData
      ,Minimum
      ,BlockSize
      ,KernelSize
      ,PreAllocate
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int /*argument*/, QFile* /*file*/) override final {}
   virtual void setArgument(int argument, EAbstractData* data) override final;
   quint32 getCapabilities() const override final
      { return Capabilities::Serial|Capabilities::OpenCL; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
   virtual int getBlockSize() override;
   virtual bool runBlock(int block) override final;
   virtual void finish() override final {}
private:
   struct Block
   {
      enum
      {
         Start
         ,Load
         ,Execute
         ,Read
         ,Done
      };
      Block(EOpenCLDevice& device, int size, int kernelSize)
      {
         references = device.makeBuffer<cl_int>(2*kernelSize).release();
         answers = device.makeBuffer<cl_float>(kernelSize).release();
         workBuffer = device.makeBuffer<cl_float>(2*size*kernelSize).release();
         rankBuffer = device.makeBuffer<cl_int>(2*size*kernelSize).release();
         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }
      ~Block()
      {
         delete references;
         delete answers;
         delete workBuffer;
         delete rankBuffer;
      }
      int state {Start};
      int x;
      int y;
      EOpenCLEvent event;
      EOpenCLBuffer<cl_int>* references;
      EOpenCLBuffer<cl_float>* answers;
      EOpenCLBuffer<cl_float>* workBuffer;
      EOpenCLBuffer<cl_int>* rankBuffer;
   };
   void initializeKernel();
   void initializeBlockExpressions();
   void initializeKernelArguments();
   void runStartBlock(Block& block);
   void runLoadBlock(Block& block);
   void runExecuteBlock(Block& block);
   void runReadBlock(Block& block);
   ExpressionMatrix* _input {nullptr};
   CorrelationMatrix* _output {nullptr};
   int _minimum {30};
   int _blockSize {4};
   int _kernelSize {4096};
   bool _allocate {false};
   Block** _blocks {nullptr};
   EOpenCLProgram* _program {nullptr};
   EOpenCLKernel* _kernel {nullptr};
   EOpenCLBuffer<cl_float>* _expressions {nullptr};
   int _x {1};
   int _y {0};
   qint64 _totalPairs;
   qint64 _pairsComplete {0};
   int _lastPercent {0};
};



#endif
