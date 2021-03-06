# Clustering into regions of synchronous variations

Joint work with  [Christophe Giraud](http://www.cmap.polytechnique.fr/~giraud/).

---

This R package implements a clustering method described in detail in [this article](http://www.cmap.polytechnique.fr/~giraud/SynchronousPop.pdf).

The goal is to divide a map of observational sites into regions of similar intra-variability over the years, but 
with distinct dynamics compared to other regions.
The problem has both spatial (the sites) and temporal (measurements every year) aspects, and is difficult because 
there are no a priori indications about the regions.

---

Two (heuristic) methods are available in the package, either direct graph clustering or parameters estimation 
in a model-based formulation of the problem.

Main function to cluster populations data: `findSyncVarRegions()`.
