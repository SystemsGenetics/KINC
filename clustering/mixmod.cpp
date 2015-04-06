#include "mixmod.h"


/**
 * Constructor for MixMod class.
 *
 * @param EMatrix *ematrix;
 */
MixMod::MixMod(EMatrix * ematrix) {
  // The data we are using is qualitative, so set the data type.
  dataType = XEM::QuantitativeData;

  // Create the Gaussian Data object and set the dataDescription object.
  gdata = new XEM::GaussianData(ematrix->getNumGenes(), ematrix->getNumSamples(), ematrix->getMatrix());
  dataDescription = new XEM::DataDescription(gdata);

  // Set the number of clusters to be tested to 9.
  nbCluster.push_back(2);
  nbCluster.push_back(3);
  nbCluster.push_back(4);
  nbCluster.push_back(5);
  nbCluster.push_back(6);
  nbCluster.push_back(7);
  nbCluster.push_back(8);
  nbCluster.push_back(9);

  // Create the ClusteringInput object.
  cInput = new XEM::ClusteringInput(nbCluster, *dataDescription);
  cOutput = NULL;

  // Finalize input: run a series of sanity checks on it
  cInput->finalize();
}

/**
 * Destructor for the MixMod class.
 */
MixMod::~MixMod() {
  delete gdata;
  delete dataDescription;
  delete cInput;
  if (cOutput) {
    delete cOutput;
  }
}

/**
 * Executes mixture model clustering.
 */
void MixMod::run() {
   // Create XEM::ClusteringMain
   XEM::ClusteringMain cMain(cInput);

   // Run XEM::ClusteringMain
   cMain.run();

   // Create a new XEM::ClusteringOutput object
   cOutput = cMain.getOutput();

   if (cOutput->atLeastOneEstimationNoError()) {
     vector<XEM::ClusteringModelOutput*> cMOutput = cOutput->getClusteringModelOutput();
     for (int i=0; i<3; i++) {
      cout << "MODEL " << i << " " << cMOutput[i]->getNbCluster() << " " << cMOutput[i]->getLikelihood() << " " << cMOutput[i]->getCriterionOutput(0).getValue() << endl;
     }
   }
}
