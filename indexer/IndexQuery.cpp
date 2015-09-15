#include "IndexQuery.h"

/**
 *
 */
IndexQuery::IndexQuery(char * indexdir, EMatrix * ematrix) {
  this->indexdir = indexdir;
  this->ematrix = ematrix;
  this->th = 0.0;
}

/**
 *
 */
IndexQuery::~IndexQuery() {

}
