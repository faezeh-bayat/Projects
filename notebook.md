Research project : Prediction of transcription factor binding sites by introducing new signal thresholds
______
Proposed method:Using mean-variance stabilization method for introducing new signals for each reads
_________
Evaluation methods:

Motif occurence 
Differentially expression of the genes in different cell types
Correlation of gene expression data with tag read signals
______

2018-05-05-Motif occurence

We had a comparison of our method with fold enrichment method and spp.

Result:
Our method showed better results in comparison with the other two methods.
_______
2018-05-22-Differentially expression of the genes in different cell types

We used the the MAnorm method for our comparison.
We have devided our comparison into four groups:

1. MAnorm on raw read
2. MAnorm on variance-stabilization signal
3. MAnorm on fold enrichment signal
4. Difference of variance-stabilization signals
________
2018-06-21-Correlation of gene expression data with transcription factors

For this part of the analysis, we are going to identify the correlation between the gene expression data and signals that we have from different methods. For this purpose, we are going to compute the signals in the window size of 10000 from the TSS,transcriptional start site,(8kb toward upstream and 2kb toward the downstream). Then we calculate the summation of each base's signal in the window. Finally we are going to plot the window's signal versus the gene expression. 
