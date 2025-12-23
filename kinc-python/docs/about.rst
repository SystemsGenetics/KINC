About KINC-Python
=================

KINC-Python is a complete Python rewrite of the Knowledge Independent Network Construction (KINC) software for constructing gene co-expression networks from gene expression data.

Project Overview
----------------

KINC-Python was developed to:

1. **Modernize KINC**: Leverage modern Python scientific libraries for easier maintenance and development
2. **Improve Accessibility**: Make KINC more accessible to researchers familiar with Python
3. **Enhance Interoperability**: Use standard HDF5 formats for better integration with Python data science tools
4. **Maintain Compatibility**: Provide statistically equivalent results to the original C++ version

Key Differences from C++ KINC
------------------------------

**Advantages of Python Version:**

- **Easier to Install**: Simple pip installation, no complex compilation
- **More Readable**: Python code is easier to understand and modify
- **Better Data Formats**: HDF5 instead of custom binary formats
- **Simpler Dependencies**: Uses well-maintained Python packages
- **Easier to Extend**: Add new correlation methods or clustering algorithms easily

**Trade-offs:**

- **Performance**: ~2-5x slower than C++ in serial mode (but parallelization helps)
- **Numerical Differences**: Uses sklearn GMM which is statistically equivalent but not bit-exact

What is Gene Co-expression Network Construction?
-------------------------------------------------

Gene co-expression networks (GCNs) are constructed by calculating correlation (similarity) between expression profiles of genes across multiple samples. These networks help identify:

- **Co-regulated genes**: Genes that work together in biological processes
- **Functional modules**: Groups of genes involved in specific pathways
- **Gene function prediction**: Infer function of unknown genes from known neighbors

Traditional vs. GMM-based Construction
---------------------------------------

**Traditional Approach:**

- Computes a single correlation value per gene pair
- Assumes all samples share the same relationship
- Simpler but may miss condition-specific relationships

**GMM-based Approach (KINC):**

- Uses Gaussian Mixture Models to identify clusters of samples
- Computes separate correlations for each cluster
- Captures **condition-specific** relationships
- More complex but reveals hidden patterns

Example: A gene pair might have:

- **Cluster 1**: High positive correlation in drought conditions
- **Cluster 2**: No correlation in normal conditions
- **Cluster 3**: Negative correlation in heat stress

Development Team
----------------

KINC-Python was developed by the SystemsGenetics group at Washington State University.

**Original KINC (C++):**

- Ben Shealy
- Josh Burns
- Stephen Ficklin
- Frank Feltus

**Python Implementation:**

- Stephen Ficklin
- Contributors welcome!

License
-------

KINC-Python is distributed under the same license as the original KINC.

See the `LICENSE file <https://github.com/SystemsGenetics/KINC-Python/blob/master/LICENSE>`_ for details.

Links
-----

- **GitHub Repository**: https://github.com/SystemsGenetics/KINC-Python
- **Original C++ KINC**: https://github.com/SystemsGenetics/KINC
- **Bug Reports**: https://github.com/SystemsGenetics/KINC-Python/issues
- **SystemsGenetics Lab**: https://ficklin.github.io/

Related Projects
----------------

- **KINC.R**: R package for post-processing KINC networks
- **GEMmaker**: Workflow for creating gene expression matrices
- **GEMprospector**: Interactive visualization of GEMs
