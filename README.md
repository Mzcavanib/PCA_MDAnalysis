# PCA_MDanalysis
Principal Component Analysis using MDAnalysis and Numpy to generate a .cvs of the PC1 and PC2 for Cytoscape analysis.

# Residue Contact Network with PCA Contributions

This repository contains a Python script that builds a **residue–residue contact network** from a molecular dynamics trajectory and enriches it with **principal component analysis (PCA) contributions**.

Features
- Loads protein trajectory data using **MDAnalysis**.
- Constructs a **contact network** of residues based on Cα atom distances.
- Counts how often each contact occurs across the trajectory (edge weights).
- Integrates **PCA contributions** from an external CSV file.
- Exports the annotated network in **GraphML format** for visualization in tools like **Gephi** or **Cytoscape**.

For visualization of the residue contribution to PCA there is another program called 'network_pca.py'.

Requirements
Make sure you have the following Python packages installed:

- [MDAnalysis](https://www.mdanalysis.org/)
- [NetworkX](https://networkx.org/)
- [Pandas](https://pandas.pydata.org/)

Then:

- Execute 'pca_csv.py' to generate the PCA contributions on '.csv' format.
- Execute 'network_pca.py' to generate the residue contact network with PCA contributions.
- Execute 'pca_res_analysis.py' to visualize the PCA contributions on PC1 and PC2.



