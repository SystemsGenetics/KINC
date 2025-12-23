Examples
========

This page provides real-world examples of using KINC-Python for gene co-expression network construction.

Example 1: Basic Workflow
--------------------------

Simple workflow for small datasets without clustering.

**Dataset**: 100 genes, 50 samples

**Step 1: Import Data**

.. code-block:: bash

   kinc import-emx --input data.txt --output data.emx.h5

**Step 2: Compute Correlations**

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --clusmethod none \
     --mincorr 0.85

**Step 3: Extract Network**

.. code-block:: bash

   kinc extract \
     --input data.cmx.h5 \
     --output network.txt \
     --mincorr 0.85

Example 2: GMM-based Condition-Specific Network
------------------------------------------------

Full workflow with GMM clustering to find condition-specific correlations.

**Dataset**: Rice RNA-seq data, 20K genes, 100 samples

**Step 1: Import**

.. code-block:: bash

   kinc import-emx \
     --input rice_rnaseq.txt \
     --output rice.emx.h5 \
     --nan "NA"

**Step 2: GMM-based Similarity**

.. code-block:: bash

   kinc similarity \
     --input rice.emx.h5 \
     --cmx rice.cmx.h5 \
     --ccm rice.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --minclus 1 \
     --maxclus 5 \
     --crit ICL \
     --minsamp 30 \
     --preout true \
     --postout true \
     --mincorr 0.85

**Step 3: Extract Network with Sample Masks**

.. code-block:: bash

   kinc extract \
     --input rice.cmx.h5 \
     --ccm rice.ccm.h5 \
     --output rice_network.txt \
     --mincorr 0.85 \
     --format text

**Output**:

.. code-block:: text

   Source   Target   Cluster   Correlation   SampleMask
   LOC_Os01g01010  LOC_Os01g01020  0  0.923  111111111100000000...
   LOC_Os01g01010  LOC_Os01g01020  1  0.887  000000000011111111...

Example 3: Python API Workflow
-------------------------------

Complete workflow using Python API:

.. code-block:: python

   from kinc.analytics import import_emx, Similarity, extract_network
   from kinc.data import ExpressionMatrix, CorrelationMatrix

   # 1. Import expression matrix
   print("Importing expression matrix...")
   import_emx("data.txt", "data.emx.h5")

   # 2. Inspect the data
   with ExpressionMatrix("data.emx.h5", mode='r') as emx:
       print(f"Genes: {emx.n_genes}")
       print(f"Samples: {emx.n_samples}")
       print(f"First 5 genes: {emx.gene_names[:5]}")

   # 3. Compute similarity with GMM
   print("Computing pairwise similarity...")
   sim = Similarity(
       input_file="data.emx.h5",
       cmx_file="data.cmx.h5",
       ccm_file="data.ccm.h5",
       correlation_method='spearman',
       clustering_method='gmm',
       min_clusters=1,
       max_clusters=5,
       criterion='ICL',
       min_samples=30,
       min_correlation=0.85
   )
   sim.run()

   # 4. Inspect results
   with CorrelationMatrix("data.cmx.h5", mode='r') as cmx:
       print(f"Total pairs saved: {cmx.pair_count}")

   # 5. Extract network
   print("Extracting network...")
   extract_network(
       cmx_file="data.cmx.h5",
       output_file="network.txt",
       ccm_file="data.ccm.h5",
       min_correlation=0.85
   )

   print("Done!")

Example 4: Expression Thresholding
-----------------------------------

Filter samples by expression level before computing correlations:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --minexpr 1.0 \
     --maxexpr 15.0 \
     --minsamp 25 \
     --mincorr 0.80

This excludes samples where either gene has expression < 1.0 or > 15.0.

Example 5: Custom Analysis
---------------------------

Access data directly for custom analysis:

.. code-block:: python

   import numpy as np
   import pandas as pd
   from kinc.data import ExpressionMatrix, CorrelationMatrix, CCMatrix

   # Read expression data
   emx = ExpressionMatrix("data.emx.h5", mode='r')
   expressions = emx.read_expressions()

   # Calculate custom statistics
   mean_expr = np.mean(expressions, axis=1)
   std_expr = np.std(expressions, axis=1)

   # Create DataFrame
   stats_df = pd.DataFrame({
       'gene': emx.gene_names,
       'mean': mean_expr,
       'std': std_expr
   })

   print(stats_df.head())

   # Read correlation results
   cmx = CorrelationMatrix("data.cmx.h5", mode='r')

   # Extract high-confidence edges
   edges = []
   for gene_x, gene_y, correlations in cmx.iter_pairs():
       for cluster, corr in enumerate(correlations):
           if not np.isnan(corr) and abs(corr) > 0.90:
               edges.append({
                   'source': emx.gene_names[gene_x],
                   'target': emx.gene_names[gene_y],
                   'correlation': corr,
                   'cluster': cluster
               })

   edges_df = pd.DataFrame(edges)
   print(f"High-confidence edges: {len(edges_df)}")

Example 6: Batch Processing
----------------------------

Process multiple datasets:

.. code-block:: python

   import glob
   from kinc.analytics import import_emx, Similarity

   # Find all GEM files
   gem_files = glob.glob("datasets/*.txt")

   for gem_file in gem_files:
       print(f"Processing {gem_file}...")

       # Generate output filenames
       base = gem_file.replace(".txt", "")
       emx_file = f"{base}.emx.h5"
       cmx_file = f"{base}.cmx.h5"
       ccm_file = f"{base}.ccm.h5"

       # Import
       import_emx(gem_file, emx_file)

       # Similarity
       sim = Similarity(
           input_file=emx_file,
           cmx_file=cmx_file,
           ccm_file=ccm_file,
           correlation_method='spearman',
           clustering_method='gmm',
           min_correlation=0.85
       )
       sim.run()

       print(f"Completed {gem_file}")

Example 7: Network Analysis
----------------------------

Analyze extracted network using NetworkX:

.. code-block:: python

   import pandas as pd
   import networkx as nx

   # Load network
   network = pd.read_csv("network.txt", sep='\t')

   # Create graph
   G = nx.Graph()

   # Add edges (using first cluster only)
   for _, row in network[network['Cluster'] == 0].iterrows():
       G.add_edge(
           row['Source'],
           row['Target'],
           weight=row['Correlation']
       )

   # Network statistics
   print(f"Nodes: {G.number_of_nodes()}")
   print(f"Edges: {G.number_of_edges()}")
   print(f"Density: {nx.density(G):.4f}")

   # Find highly connected genes
   degree = dict(G.degree())
   top_genes = sorted(degree.items(), key=lambda x: x[1], reverse=True)[:10]

   print("\nTop 10 hub genes:")
   for gene, deg in top_genes:
       print(f"  {gene}: {deg} connections")

   # Community detection
   from networkx.algorithms import community
   communities = community.greedy_modularity_communities(G)

   print(f"\nCommunities detected: {len(communities)}")
   for i, comm in enumerate(communities[:5]):
       print(f"  Community {i+1}: {len(comm)} genes")

Example 8: Visualization
-------------------------

Visualize network using matplotlib:

.. code-block:: python

   import pandas as pd
   import networkx as nx
   import matplotlib.pyplot as plt

   # Load network
   network = pd.read_csv("network.txt", sep='\t')
   network = network[network['Cluster'] == 0]  # First cluster only

   # Create graph
   G = nx.Graph()
   for _, row in network.iterrows():
       if abs(row['Correlation']) > 0.90:  # High correlations only
           G.add_edge(row['Source'], row['Target'],
                     weight=abs(row['Correlation']))

   # Draw network
   plt.figure(figsize=(12, 12))
   pos = nx.spring_layout(G, k=0.5, iterations=50)

   # Node size based on degree
   node_size = [300 * G.degree(node) for node in G.nodes()]

   # Edge width based on correlation
   edges = G.edges()
   weights = [G[u][v]['weight'] for u, v in edges]
   edge_width = [3 * w for w in weights]

   nx.draw_networkx(
       G, pos,
       node_size=node_size,
       width=edge_width,
       with_labels=True,
       font_size=8,
       node_color='lightblue',
       edge_color='gray',
       alpha=0.7
   )

   plt.title("Gene Co-expression Network (r > 0.90)")
   plt.axis('off')
   plt.tight_layout()
   plt.savefig("network.png", dpi=300, bbox_inches='tight')
   plt.show()

Example 9: Comparison Analysis
-------------------------------

Compare traditional vs GMM-based networks:

.. code-block:: python

   from kinc.analytics import Similarity, extract_network

   # Traditional (no clustering)
   print("Computing traditional network...")
   sim_trad = Similarity(
       input_file="data.emx.h5",
       cmx_file="traditional.cmx.h5",
       ccm_file="traditional.ccm.h5",
       correlation_method='spearman',
       clustering_method='none',
       min_correlation=0.85
   )
   sim_trad.run()

   # GMM-based
   print("Computing GMM-based network...")
   sim_gmm = Similarity(
       input_file="data.emx.h5",
       cmx_file="gmm.cmx.h5",
       ccm_file="gmm.ccm.h5",
       correlation_method='spearman',
       clustering_method='gmm',
       min_clusters=1,
       max_clusters=5,
       min_correlation=0.85
   )
   sim_gmm.run()

   # Extract both
   extract_network("traditional.cmx.h5", "traditional_network.txt", min_correlation=0.85)
   extract_network("gmm.cmx.h5", "gmm_network.txt", min_correlation=0.85)

   # Compare
   import pandas as pd

   trad = pd.read_csv("traditional_network.txt", sep='\t')
   gmm = pd.read_csv("gmm_network.txt", sep='\t')

   print(f"Traditional edges: {len(trad)}")
   print(f"GMM edges: {len(gmm)}")
   print(f"GMM found {len(gmm) - len(trad)} additional edges")

Complete Workflow Script
-------------------------

Save this as ``kinc_workflow.sh``:

.. code-block:: bash

   #!/bin/bash
   #
   # Complete KINC workflow for gene co-expression network construction
   #

   # Configuration
   INPUT_GEM="rice_expression.txt"
   PREFIX="rice"
   MIN_CORR=0.85

   # Create output directory
   mkdir -p results

   echo "=== KINC Workflow Started ==="
   date

   # Step 1: Import
   echo "Step 1: Importing expression matrix..."
   kinc import-emx \
     --input "$INPUT_GEM" \
     --output "results/${PREFIX}.emx.h5"

   # Step 2: Similarity with GMM
   echo "Step 2: Computing pairwise similarity..."
   kinc similarity \
     --input "results/${PREFIX}.emx.h5" \
     --cmx "results/${PREFIX}.cmx.h5" \
     --ccm "results/${PREFIX}.ccm.h5" \
     --corrmethod spearman \
     --clusmethod gmm \
     --minclus 1 \
     --maxclus 5 \
     --crit ICL \
     --minsamp 30 \
     --preout true \
     --postout true \
     --mincorr $MIN_CORR

   # Step 3: Extract network
   echo "Step 3: Extracting network..."
   kinc extract \
     --input "results/${PREFIX}.cmx.h5" \
     --ccm "results/${PREFIX}.ccm.h5" \
     --output "results/${PREFIX}_network.txt" \
     --mincorr $MIN_CORR \
     --format text

   echo "=== Workflow Complete ==="
   date
   echo "Output files in results/ directory"

Make it executable and run:

.. code-block:: bash

   chmod +x kinc_workflow.sh
   ./kinc_workflow.sh

See Also
--------

- :doc:`quickstart` for a quick introduction
- :doc:`usage` for detailed parameter documentation
- :doc:`api` for Python API reference
