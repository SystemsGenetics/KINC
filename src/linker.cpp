#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"
#include "spearman.h"



class KINCFactory : public Ace::Factory
{
   Ace::Analytic* build_analytic(const std::string&) override final { return nullptr; }
   Ace::Data* build_data(const std::string& type) override final
   {
      Ace::Data* ret {nullptr};
      if (type==std::string("emx"))
      {
         ret = new EMatrix;
      }
      return ret;
   }
};
