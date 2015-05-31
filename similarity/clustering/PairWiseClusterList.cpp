#include "PairWiseClusterList.h"

/**
 *
 */
PairWiseClusterList::PairWiseClusterList(PairWiseSet * pwset) {
  this->pwset = pwset;
  this->num_clusters = 0;
  this->head = NULL;
}
/**
 *
 */
PairWiseClusterList::~PairWiseClusterList() {
  PairWiseCluster * next;

  PairWiseCluster * curr = head;
  while (curr) {
    next = curr->neighbor;
    delete curr;
    curr = next;
  }
}

/**
 * Adds a cluster to the cluster list.
 *
 * Upon deletion of the list the memory assocaited with the cluster will
 * be removed as well.
 */
void PairWiseClusterList::addCluster(PairWiseCluster * pwc) {

  this->num_clusters++;
  pwc->index = this->num_clusters;

  // Check the list to see if it is empty. If so, then make this item the
  // new head.
  if (head == NULL) {
    head = pwc;
    return;
  }

  // Traverse the list to the end and then add the new item
  PairWiseCluster * curr = head;
  while (curr->neighbor != NULL) {
    curr = curr->neighbor;
  }
  curr->neighbor = pwc;
}
