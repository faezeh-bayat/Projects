
Todo:
- Model improvements:
    + [4] Add improvement that involves combining all zero-value bins.
    + [6] Speedup from the binning improvement #3 described in theory.tex. (Or #4, which is easier)
- Quantitative evaluation:
    + [2] Measure the quality of fit of a transformation. Transform two replicates, R1 and R2 to produce t(R1) and t(R2). Plot: x = t(R1). y = t(R2) - t(R1). There should not be a trend (sloping up or down) to the points. 
    + [6] Can we figure out what is the theoretically maximum possible differential expression evaluation value?
- Qualitative evaluation (visualization etc):
    + [3] New plot (similar to the previous one): Each point is a gene. Horizontal axis = ChIP-seq signal in cell type X. Vertical axis = ChIP-seq signal in cell type Y. Color = difference in gene expression between X and Y (orange if X>>Y, teal if Y>>X, black if X~=Y). Two panels: one for fold-enrichment, one for variance-stabilized signal. 
    + [3-Brian] UCSC tracks: The UCSC genome browser is an online tool for visualizing genomics data. Please convert a H3K4me3 ChIP-seq track (from your favorite cell type) into BigBed format and upload it to the UCSC genome browser. Do this for both fold enrichment and variance-stabilized data. Neda is doing this for the SSM model so she can help you; I can also help when we meet. I will help you find a region that illustrates our method well and we will take a screenshot to include in the paper. 

- Presentation (writing the paper, publishing code etc): 
    + [2-Max] Revise the manuscript. 

Finished:
- [1] Unify the code so that you pick a set of parameters and then the code computes both the differential expression and gene expression analyses, and produces all the plots.
- [1] Organize the project according to the principles in [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)
- Model improvements:
    + [1] Try the 3 different input/output methods listed in the theory.tex document. 
    + [4] Think about whether we can come up with a method that treats both replicates as the same.
    + [4] Think about different modes: (1) The user has done just one experiment (no replicates) and needs to use a 
    mean-variance curve from something else or (2) the user has done multiple replicates and wants to use 
    their own mean-variance curve.
    + [4] What should the model do if there are more than two replicates (3, 4, etc)?
    + [4] Try a version where we learn the mean-variance relationship on one data set and use it to transform another data set. 
    + [3] Parallelize the code on Compute Canada. 
    + 
- Quantitative evaluation:
- Qualitative evaluation (visualization etc):
    + [2] New plot: This plot is meant to be a visual version of the diff-expr evaluation. Each point is a gene. Horizontal axis = difference in ChIP-seq signal between cell type X and Y. (Don't take the absolute value; there should be both positive and negative values). Vertical axis = difference in gene expression value between X and Y. Two panels: one for fold-enrichment, one for variance-stabilized signal. Again, use the same ChIP-seq track as in your evaluation and pick your favorite cell type pair. 




Meaning of priority numbers:
- 1: Key parts of our core contribution. Examples: Implementing the main method, doing the main evaluation against our main competitors. 
- 2: Important supporting contributions. Examples: Comparing main variants of our method, delving into why we achieve better performance. 
- 3: Steps that aren’t part of our main contribution, but that are still necessary. Examples: Making sure that our objective function is being optimized properly; documenting our code.
- 4: The next step for improving our contribution. Not necessary if we have a publishable unit, but this should be our next step otherwise. 
- 5: Ideas for improving our contribution. If what we’re doing is working well then we can save these for a follow-up paper. If it’s not working, we should look here for ideas. 
- 6: Not worth the time now, but recording the idea in case it comes in handy later. 
