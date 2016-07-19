#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"
#include "spearman.h"



ACE_BEGIN_DATA
ACE_DATA_PLUGIN(emx,EMatrix)
ACE_DATA_PLUGIN(cmx,CMatrix)
ACE_END_DATA



ACE_BEGIN_ANALYTIC
ACE_ANALYTIC_PLUGIN(spearman,Spearman)
ACE_END_ANALYTIC
