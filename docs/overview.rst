Software Overview
=================


In KINC, GMM clusters undergo correlation analysis independent of one another. This is different from traditional approaches where all samples, which may not meet the assumptions of the correlation test, are analyzed together.  Methods include:

- *Pearson*
- *Spearman*

.. note::

  KINC does not implement Mutual Information, another common "association" method for GCN construction as literature has shown no significant improvement over methods such as Pearson's or Spearman.

**Filtering**

- *Correlation Power Analysis* removes non-significant edges. Correlations from GMM clusters that had insufficient samples to "trust" a given correlation are removed.

**Thresholding**

- *Power-law* thresholding iterates through decreasing correlation scores, and checks if the resulting network follows a power-law distribution based on node degree. A network that follows a power-law distribution is known as a `scale-free` network. A threshold can be identified at a given correlation score where the resulting network does not appear scale-free.
- *Random Matrix Theory (RMT)* also iterates through decreasing correlation scores, the resulting similarity matrix is examined for properties of random matricies. If the resulting matrix ceases to look non-random a threshold is identified.

.. note::

  A user of KINC may decide to apply one or more of the above thresholding approaches.

**Module Discovery**




How does KINC work?
-------------------






Pairwise Clustering
```````````````````
- *Gaussian Mixture Models (GMMs)* is performed prior to correlation analysis, modes of gene expression are identified. These modes potentially represent condition-specific gene relationships and become `edges` in teh network if significant.
