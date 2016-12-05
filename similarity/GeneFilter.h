#ifndef _GENEFILTER_
#define _GENEFILTER_

// This structure is used for containing the set of genes for filtering.
struct geneFilter {
  char * filename;
  int * indicies;
  int num_genes;
};

#endif
