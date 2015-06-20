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
    char * method;
    // The prefix for the filename.
    char * fileprefix;
    // A unique id to differentiate between parallel executions.
    int job_index;

    // Holds the last x and y values that exist in each file before
    // the job terminated without completion.
    int * last_x;
    int * last_y;
    // The position in the file where the last_x and last_y are found.
    int * last_seek;

    // Opens and creates file pointers for all of the
    void openOutFiles();
    void closeOutFiles();
    void findLastPositions();

  public:
    // Constructor.
    PairWiseClusterWriter(char * method, char * fileprefix, int id, int num_samples);
    // Destructor.
    ~PairWiseClusterWriter();
    // Writes a PairWiseCluster to the proper file.
    void writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2);
};

#endif
