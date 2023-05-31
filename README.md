# Bacterial lifestyle shapes pangenomes
Authors: Anna Dewar, Chunhui Hao, Laurence Belcher, Melanie Ghoul, Stuart West
Affiliation: Department of Biology, University of Oxford, United Kingdom

## Overview
This repository contains code and data for the research article 'Bacterial lifestyle shapes pangenomes'. For any queries please contact Anna Dewar at anna.dewar@biology.ox.ac.uk.

## Supplementary Material 1

S1 contains all supplementary methods and results. The document was compiled from an Rmarkdown file. The original .Rmd file, including code for all models and results in the S1 document, is available in this repository as: 'Code_S1.Rmd'.

### Data & tree
The .Rmd file requires the data file 'pangenome_lifestyles.csv', which contains all data used in our analyses, and the phylogeny 'panX_tree.nex'. For ease of use we recommend saving all files, including the data, tree and .Rmd file, into a folder alongside an Rstudio project. 

### Compiling 'Code_S2.Rmd' as a pdf
If you wish to locally compile the file 'Code_S1.pdf' into a PDF identical to the S1 document, please also download the file 'figure_order_header.tex' and include this within the same folder as the .Rmd file.

## Phylogeny
To build our phylogeny, we used a recently published maximum likelihood tree generated with 16S ribosomal protein data as the basis for our phylogeny (Hug et al. (2016), 'A new view of the tree of life', Nat. Microbiol.). We used the R package ‘ape’ to identify all branches that matched either a species or a genus in our dataset. In cases where we had multiple species within a single genus, we used the R package ‘phytools’ to add these species as additional branches in the tree. We used published phylogenies from the literature to add any within-genus clustering of species’ branches (details and references of these phylogenies are available in 'Supp_material_3.xlsx').

The code for how did this is available to download as 'tree_script.R'. It requires the files 'Original_tree.txt' and 'pangenome_species.csv'. Each line of the script edits the tree to produce the final tree, so the lines of code should be run in the order they are in the script.

