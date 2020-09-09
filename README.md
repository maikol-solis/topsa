# TopSA (Topological Sensitivity Analysis)

To install the latest version use the instruction:

```
remotes::install_github("maikol-solis/topsa@develop")
```

This package estimates geometrical sensitivity indices reconstructing the embedding manifold of the data. The reconstruction is done via a Vietoris Rips with a fixed radius. With the homology of order 2 we estimate symmetric reflections of the those manifold to determine a sensitivity index on of the input model.  

A detail explanation of the method is covered [here](https://www.dropbox.com/s/0kcrjhhkl7899n1/article-symmetric-reflection.pdf?dl=0)
