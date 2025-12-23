Citation
========

If you use KINC-Python in your research, please cite the following publications:

Primary Citations
-----------------

**Original KINC Publication:**

   Ficklin, S.P., Dunwoodie, L.J., Poehlman, W.L. et al. (2017).
   "Discovering Condition-Specific Gene Co-Expression Patterns Using Gaussian Mixture Models: A Cancer Case Study."
   *Scientific Reports* 7, 8617.
   https://doi.org/10.1038/s41598-017-09094-4

**GPU Implementation:**

   Shealy, B.T., Burns, J.J.R., Smith, M.C., Feltus, F.A., & Ficklin, S.P. (2019).
   "GPU Implementation of Pairwise Gaussian Mixture Models for Multi-Modal Gene Co-Expression Networks."
   *IEEE Access* 7, 160845-160857.
   https://doi.org/10.1109/ACCESS.2019.2951284

**Addressing Noise in GCNs:**

   Burns, J.J.R., Shealy, B.T., Greer, M.S., Hadish, J.A., McGowan, M.T., Biggs, T., Smith, M.C., Feltus, F.A., & Ficklin, S.P. (2022).
   "Addressing Noise in Co-Expression Network Construction."
   *Briefings in Bioinformatics* 23(1), bbab495.
   https://doi.org/10.1093/bib/bbab495

BibTeX Entries
--------------

Original KINC:

.. code-block:: bibtex

   @article{ficklin2017discovering,
     title={Discovering Condition-Specific Gene Co-Expression Patterns Using Gaussian Mixture Models: A Cancer Case Study},
     author={Ficklin, Stephen P and Dunwoodie, Liam J and Poehlman, William L and Watson, Carolyn and Roche, Katherine E and Feltus, F Alex},
     journal={Scientific Reports},
     volume={7},
     number={1},
     pages={8617},
     year={2017},
     publisher={Nature Publishing Group},
     doi={10.1038/s41598-017-09094-4}
   }

GPU Implementation:

.. code-block:: bibtex

   @article{shealy2019gpu,
     title={GPU Implementation of Pairwise Gaussian Mixture Models for Multi-Modal Gene Co-Expression Networks},
     author={Shealy, Ben T and Burns, John JR and Smith, Melissa C and Feltus, F Alex and Ficklin, Stephen P},
     journal={IEEE Access},
     volume={7},
     pages={160845--160857},
     year={2019},
     publisher={IEEE},
     doi={10.1109/ACCESS.2019.2951284}
   }

Addressing Noise:

.. code-block:: bibtex

   @article{burns2022addressing,
     title={Addressing Noise in Co-Expression Network Construction},
     author={Burns, John JR and Shealy, Ben T and Greer, Melissa S and Hadish, John A and McGowan, Matthew T and Biggs, Tyler and Smith, Melissa C and Feltus, F Alex and Ficklin, Stephen P},
     journal={Briefings in Bioinformatics},
     volume={23},
     number={1},
     pages={bbab495},
     year={2022},
     publisher={Oxford University Press},
     doi={10.1093/bib/bbab495}
   }

Related Publications
--------------------

These publications use or extend KINC methodology:

**KINC.R Package:**

   Poehlman, W.L., Rynge, M., Baxter, D.J., Balamurugan, D., Hadish, J.A., Petersen, T.N., & Feltus, F.A. (2019).
   "OSG-KINC: High-throughput gene co-expression network construction using the Open Science Grid."
   *Bioinformatics and Biology Insights* 13, 1177932219828776.

**GEMmaker:**

   Hadish, J.A., Biggs, T.D., Shealy, B.T., Bender, M.R., McKnight, C.J., Wytko, C., Meng, A., Kinser, S.R., Feltus, F.A., & Ficklin, S.P. (2022).
   "GEMmaker: process big RNA-seq data on heterogeneous computational infrastructure."
   *BMC Bioinformatics* 23(1), 156.

Software Citation
-----------------

If you want to cite the KINC-Python software specifically:

   KINC-Python Development Team (2024).
   "KINC-Python: Python Implementation of Knowledge Independent Network Construction."
   Version 3.4.2.
   https://github.com/SystemsGenetics/KINC-Python

Dependencies
------------

KINC-Python builds on several excellent open-source libraries that should also be acknowledged:

- **NumPy**: Harris et al. (2020). *Nature* 585, 357-362.
- **SciPy**: Virtanen et al. (2020). *Nature Methods* 17, 261-272.
- **scikit-learn**: Pedregosa et al. (2011). *Journal of Machine Learning Research* 12, 2825-2830.
- **pandas**: McKinney (2010). *Proceedings of the 9th Python in Science Conference*.
- **HDF5/h5py**: The HDF Group (1997-2024).

How to Cite in Different Contexts
----------------------------------

**For General KINC Use:**

"We constructed gene co-expression networks using KINC (Ficklin et al., 2017), a tool that uses Gaussian Mixture Models to identify condition-specific correlations."

**For GMM-based Networks:**

"We employed the GMM-based approach implemented in KINC (Ficklin et al., 2017; Shealy et al., 2019) to identify condition-specific gene co-expression relationships."

**For Python Implementation:**

"We used KINC-Python, a Python implementation of the KINC software (Ficklin et al., 2017), which leverages modern scientific Python libraries for gene co-expression network construction."

**For Outlier Removal:**

"Outliers were identified and removed using Tukey's method as described in Burns et al. (2022)."

Acknowledgments
---------------

KINC-Python development was supported by:

- National Science Foundation (NSF)
- Washington State University
- Clemson University
- The SystemsGenetics research group

The original C++ KINC was developed with support from NSF awards:

- NSF IOS-1127112
- NSF ABI-1062472
- NSF DBI-1659300

Contact
-------

For questions about citing KINC-Python or collaborations:

- **Stephen Ficklin**: stephen.ficklin@wsu.edu
- **GitHub Issues**: https://github.com/SystemsGenetics/KINC-Python/issues

License
-------

KINC-Python is released under the MIT License, making it free for academic and commercial use.

See the `LICENSE file <https://github.com/SystemsGenetics/KINC-Python/blob/master/LICENSE>`_ for details.
