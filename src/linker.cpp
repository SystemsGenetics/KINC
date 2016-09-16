#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"
#include "spearman.h"
#include "rmt.h"



class KINCFactory : public Ace::Factory
{
   Ace::Analytic* build_analytic(const std::string& type) override final {
      Ace::Analytic* ret {nullptr};
      if (type==std::string("spearman"))
      {
         ret = new Spearman;
      }
      else if (type==std::string("rmt"))
      {
         ret = new RMT;
      }
      return ret;
   }
   Ace::Data* build_data(const std::string& type) override final
   {
      Ace::Data* ret {nullptr};
      if (type==std::string("emx"))
      {
         ret = new EMatrix;
      }
      else if (type==std::string("cmx"))
      {
         ret = new CMatrix;
      }
      return ret;
   }
};
