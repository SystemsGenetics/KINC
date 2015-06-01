#ifndef _PWCLUSTERWRITER_
#define _PWCLUSTERWRITER_

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>

#include "PairWiseClusterList.h"

using namespace std;

/**
 *
 */
class PairWiseClusterWriter {
  private:
    // An array of file pointers.
    //FILE ** fps;
    ofstream ** fps;
    // Specifies the correlation method: pc, mi, sc
    char * method;
    // The prefix for the filename.
    char * fileprefix;
    // A unique id to differentiate between parallel executions.
    int id;

    // Opens and creates file pointers for all of the
    void openOutFiles();
    void closeOutFiles();
  public:
    // Constructor.
    PairWiseClusterWriter(char * method, char * fileprefix, int id);
    // Destructor.
    ~PairWiseClusterWriter();
    // Writes a PairWiseCluster to the proper file.
    void writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2);
};

#endif
