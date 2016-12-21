#include <ace.h>
#include "pearson.h"
#include "ematrix.h"
#include "cmatrix.h"
#include "spearman.h"
#include "rmt.h"




// This is where the different plugin/plugin utilities are invoked through the
// ACE terminal interface

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
      else if (type==std::string("pearson"))
      {
        ret = new Pearson;
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
