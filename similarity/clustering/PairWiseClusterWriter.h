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
    // The number of samples that will be written.
    int num_samples;
    // An array of file pointers.
    //FILE ** fps;
    fstream ** fps;
    // Specifies the correlation method: pc, mi, sc
    char ** method;
    // Indicates the number of methods.
    int num_methods;
    // The prefix for the filename.
    char * fileprefix;
    // A unique id to differentiate between parallel executions.
    int job_index;

    int recovery_x;
    int recovery_y;

    // Opens and creates file pointers for all of the
    void openOutFiles();
    void closeOutFiles();
    void findLastPositions();

  public:
    // Constructor.
    PairWiseClusterWriter(char ** method, int num_methods, char * fileprefix,
        int id, int num_samples);
    // Destructor.
    ~PairWiseClusterWriter();
    // Getters.
    int getRecoveryX() {return recovery_x;}
    int getRecoveryY() {return recovery_y;}
    // Writes a PairWiseCluster to the proper file.
    void writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2, float *min_sim);

};

#endif
