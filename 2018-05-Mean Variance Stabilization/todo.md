
Todo:
- Model improvements:
    + [4] Think about whether we can come up with a method that treats both replicates as the same.
    + What value should we transform? Should we transform the rep1 value, the rep2 value, or mean(rep1, rep2)? 
    + Think about different modes: (1) The user has done just one experiment (no replicates) and needs to use a 
    mean-variance curve from something else or (2) the user has done multiple replicates and wants to use 
    their own mean-variance curve.
    + What should the model do if there are more than two replicates (3, 4, etc)?

- Quantitative evaluation (gene expression, enhanacer):

- Qualitative evaluation (visualization etc):

- Comparison methods

- Presentation (writing the paper, publishing code etc): 

Meaning of priority numbers:
- 1: Key parts of our core contribution. Examples: Implementing the main method, doing the main evaluation against our main competitors. 
- 2: Important supporting contributions. Examples: Comparing main variants of our method, delving into why we achieve better performance. 
- 3: Steps that aren’t part of our main contribution, but that are still necessary. Examples: Making sure that our objective function is being optimized properly; documenting our code.
- 4: The next step for improving our contribution. Not necessary if we have a publishable unit, but this should be our next step otherwise. 
- 5: Ideas for improving our contribution. If what we’re doing is working well then we can save these for a follow-up paper. If it’s not working, we should look here for ideas. 
- 6: Not worth the time now, but recording the idea in case it comes in handy later. 
