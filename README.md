# mito-filtering
When we were doing my first scRNA-seq analyses, it became clear that some cells had a very high amount of expression of mitochondrial genes, which almost always indicates that the cell has been ruptured/damaged during the isolation process. 
People usually throw out all cells with more than x% of reads from mtRNA, usually 5% or 10%, but that was extremely stringent for our tumors and didn't leave us with much to analyze.
We ended up looking at library complexity along with % mitochondria of individual cells, and generating a linear mixture model to determine if a cell was likely displaying characteristics of a healthy or a compromised cell.
Now, we're explanding that work into a tool called miQC, which has numerous advantages over the simple cutoff approach.

This analysis once lived in greenelab/sc-cancer-hgsc/analyses/differential_expression/smart_QC.Rmd.
