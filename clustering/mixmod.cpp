#include "mixmod.h"


/**
 * Constructor for MixMod class.
 *
 * @param EMatrix *ematrix;
 */
MixModClusters::MixModClusters(PairWiseSet *pwset, int min_obs) {
  this->pwset = pwset;
  this->min_obs = min_obs;
  this->labels = NULL;
  this->pwcl = NULL;
  this->gdata = NULL;
  this->dataDescription = NULL;
  this->cInput = NULL;
  this->cOutput = NULL;
  this->data = NULL;

  // Create the PairWiseClusterList object
  this->pwcl = new PairWiseClusterList(pwset);

  // Make sure we have the correct number of observations before preparing
  // for the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // The data we are using is qualitative, so set the data type.
  this->dataType = XEM::QuantitativeData;

  // Create the Gaussian Data object and set the dataDescription object.
  this->data = (double **) malloc(sizeof(double *) * pwset->n_clean);
  for (int i = 0; i < pwset->n_clean ; i++) {
    this->data[i] = (double *) malloc(sizeof(double) * 2);
    this->data[i][0] = pwset->x_clean[i];
    this->data[i][1] = pwset->y_clean[i];
  }

  this->gdata = new XEM::GaussianData(pwset->n_clean, 2, data);
  this->dataDescription = new XEM::DataDescription(gdata);

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
  this->cInput = new XEM::ClusteringInput(nbCluster, *dataDescription);
  this->cOutput = NULL;

  // Finalize input: run a series of sanity checks on it
  cInput->finalize();
}

/**
 * Destructor for the MixMod class.
 */
MixModClusters::~MixModClusters() {
  if (gdata) {
    delete gdata;
  }
  if (dataDescription) {
    delete dataDescription;
  }
  if (cInput) {
    delete cInput;
  }
  if (cOutput) {
    delete cOutput;
  }
  // Free the data array.
  if (data) {
    for (int i = 0; i < pwset->n_clean ; i++) {
      free(data[i]);
    }
    free(data);
  }
  if (labels) {
    delete labels;
  }
  delete pwcl;
}

/**
 * Executes mixture model clustering.
 */
void MixModClusters::run() {
  // Make sure we have the correct number of observations before performing
  // the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // Create XEM::ClusteringMain
  XEM::ClusteringMain cMain(cInput);

  // Run XEM::ClusteringMain
  cMain.run();

  // Create a new XEM::ClusteringOutput object
  cOutput = cMain.getOutput();

  //
  cOutput->sort(XEM::BIC);

  if (cOutput->atLeastOneEstimationNoError()) {
    // Get the best XEM::ClusteringModelOutput
    XEM::ClusteringModelOutput* cMOutput = cOutput->getClusteringModelOutput().front();
//    XEM::ParameterDescription* paramDescription = cMOutput->getParameterDescription();
//
//    cout << "-----------------------------------------------------------------------" << endl;
//    cout << "Best model is " << endl;
//    cout << " - nbCluster : " << paramDescription->getNbCluster() << endl;
//    cout << "-----------------------------------------------------------------------" << endl << endl;
//
//    cout << "-----------------------------------------------------------------------" << endl;
//    cout << "Parameters display" << endl;
//
//    XEM::Parameter* param = paramDescription->getParameter();
//    // print out parameters
//    param->edit();
//    // print out criterion values
//    for (int64_t iCriterion = 0; iCriterion < cInput->getCriterionName().size(); iCriterion++) {
//      cMOutput->getCriterionOutput (cInput->getCriterionName (iCriterion)).editTypeAndValue (std::cout);
//    }

    // Get the classification (cluster) labels for the samples.
    XEM::LabelDescription * ldescription = cMOutput->getLabelDescription();
    XEM::Label * label = ldescription->getLabel();
    this->labels = label->getTabLabel();

    // Create the cluster objects.
    int cluster_num = 1;
    int cluster_samples[this->pwset->n_clean];
    bool done = false;
    while (!done) {
      done = true;
      for (int i = 0; i < this->pwset->n_clean; i++) {
        if (this->labels[i] == cluster_num) {
          done = false;
          cluster_samples[i] = 1;
        }
        else {
          cluster_samples[i] = 0;
        }
      }
      // if we found samples with the current cluster_num then create a
      // cluster and add it to the list
      if (!done) {
        PairWiseCluster * cluster = new PairWiseCluster(this->pwset);
        cluster->setClusterSamples(cluster_samples, true);
        cluster->performSimilarity("sc", this->min_obs);
        //cluster->printCluster();
        this->pwcl->addCluster(cluster);
      }
      cluster_num++;
    }

    // Second,
    //cout << endl;
  }
  //cout << "-----------------------------------------------------------------------" << endl;
}
