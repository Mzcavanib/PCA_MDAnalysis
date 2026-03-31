import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import pca
import pandas as pd

# 1. Cargar universo y seleccionar Cα
u = mda.Universe("final.gro", "final.xtc")
ca = u.select_atoms("protein and name CA")

# 2. Ejecutar PCA
pca_analysis = pca.PCA(u, select="protein and name CA").run()

# 3. Obtener componentes principales
try:
    p_components = pca_analysis.results.p_components
except AttributeError:
    p_components = pca_analysis.p_components

# 4. Reorganizar en vectores (x,y,z) por átomo
n_atoms = len(ca)
pc1 = p_components[:, 0].reshape(n_atoms, 3)
pc2 = p_components[:, 1].reshape(n_atoms, 3)

# 5. Magnitud de contribución por átomo
pc1_magnitudes = np.linalg.norm(pc1, axis=1)
pc2_magnitudes = np.linalg.norm(pc2, axis=1)

# 6. Crear tabla con residuo y contribuciones
resids = ca.resids
resnames = ca.resnames
df = pd.DataFrame({
    "resid": resids,
    "resname": resnames,
    "PC1_contrib": pc1_magnitudes,
    "PC2_contrib": pc2_magnitudes
})

# 7. Guardar como CSV
df.to_csv("residue_contributions.csv", index=False)
print("Archivo 'residue_contributions.csv' generado.")

