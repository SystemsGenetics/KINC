#include "mixmod.h"


/**
 * Constructor for MixMod class.
 *
 * @param EMatrix *ematrix;
 */
MixMod::MixMod(double *a, double *b, int n) {
  size = n;

  // The data we are using is qualitative, so set the data type.
  dataType = XEM::QuantitativeData;

  // Create the Gaussian Data object and set the dataDescription object.
  data = (double **) malloc(sizeof(double *) * n);
  for (int i = 0; i < n ; i++) {
    data[i] = (double *) malloc(sizeof(double) * 2);
    data[i][0] = a[i];
    data[i][1] = b[i];
  }

  gdata = new XEM::GaussianData(n, 2, data);
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
  // Free the data array.
  for (int i = 0; i < size ; i++) {
    free(data[i]);
  }
  free(data);
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

   //
   cOutput->sort (XEM::BIC);

   if (cOutput->atLeastOneEstimationNoError()) {
   // get the best XEM::ClusteringModelOutput
     XEM::ClusteringModelOutput* cMOutput = cOutput->getClusteringModelOutput().front();
     XEM::ParameterDescription* paramDescription = cMOutput->getParameterDescription();

     cout << "-----------------------------------------------------------------------" << endl;
     cout << "Best model is " << endl;
     cout << " - nbCluster : " << paramDescription->getNbCluster() << endl;
     cout << "-----------------------------------------------------------------------" << endl << endl;

     cout << "-----------------------------------------------------------------------" << endl;
     cout << "Parameters display" << endl;

     XEM::Parameter* param = paramDescription->getParameter();
     // print out parameters
     param->edit();
     // print out criterion values
     for (int64_t iCriterion = 0; iCriterion < cInput->getCriterionName().size(); iCriterion++)
      cMOutput->getCriterionOutput (cInput->getCriterionName (iCriterion)).editTypeAndValue (std::cout);
    }
    cout << "-----------------------------------------------------------------------" << endl;
}
