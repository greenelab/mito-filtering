## Scripts:
##	prepare_our_data.R: takes the raw matrix of gene by cell counts returned
##		by alevin, and prepares it as a SingleCellExperiment object with
##		basic QC run and percent mitochondrial reads calculated.
##	generate_model.R: runs a mixture model on the standardized data, generates
##		some explanatory plots, and outputs a keep/don't keep filter per cell.
##	generate_model_UMAP.R: 
## Parameters:
##	-f, --file <16030X2>: which tumor you want to run the given script on.
##	-m, --model <linear>: what kind of mixture model you want to use to filter. Current
##	             options are linear, spline, and polynomial.

#knit RMarkdown

#Rscript prepare_our_data.R -f 16030X2
#Rscript prepare_our_data.R -f 16030X3
#Rscript prepare_our_data.R -f 16030X4

Rscript generate_model.R -f 16030X2 -m linear
Rscript generate_model.R -f 16030X2 -m spline
Rscript generate_model.R -f 16030X2 -m polynomial
Rscript generate_model_UMAP.R -f 16030X2 -m linear
Rscript generate_model_UMAP.R -f 16030X2 -m spline
Rscript generate_model_UMAP.R -f 16030X2 -m polynomial
Rscript generate_cutoff_UMAP.R -f 16030X2


Rscript generate_model.R -f 16030X3 -m linear
Rscript generate_model.R -f 16030X3 -m spline
Rscript generate_model.R -f 16030X3 -m polynomial
Rscript generate_model_UMAP.R -f 16030X3 -m linear
Rscript generate_model_UMAP.R -f 16030X3 -m spline
Rscript generate_model_UMAP.R -f 16030X3 -m polynomial
Rscript generate_cutoff_UMAP.R -f 16030X3


Rscript generate_model.R -f 16030X4 -m linear
Rscript generate_model.R -f 16030X4 -m spline
Rscript generate_model.R -f 16030X4 -m polynomial
Rscript generate_model_UMAP.R -f 16030X4 -m linear
Rscript generate_model_UMAP.R -f 16030X4 -m spline
Rscript generate_model_UMAP.R -f 16030X4 -m polynomial
Rscript generate_cutoff_UMAP.R -f 16030X4


Rscript generate_model.R -f Bach -m linear
Rscript generate_model.R -f Bach -m spline
Rscript generate_model.R -f Bach -m polynomial
Rscript generate_model_UMAP.R -f Bach -m linear
Rscript generate_model_UMAP.R -f Bach -m spline
Rscript generate_model_UMAP.R -f Bach -m polynomial
Rscript generate_cutoff_UMAP.R -f Bach


Rscript generate_model.R -f Buettner -m linear
Rscript generate_model.R -f Buettner -m spline
Rscript generate_model.R -f Buettner -m polynomial
Rscript generate_model_UMAP.R -f Buettner -m linear
Rscript generate_model_UMAP.R -f Buettner -m spline
Rscript generate_model_UMAP.R -f Buettner -m polynomial
Rscript generate_cutoff_UMAP.R -f Buettner


Rscript generate_model.R -f Campbell -m linear
Rscript generate_model.R -f Campbell -m spline
Rscript generate_model.R -f Campbell -m polynomial
Rscript generate_model_UMAP.R -f Campbell -m linear
Rscript generate_model_UMAP.R -f Campbell -m spline
Rscript generate_model_UMAP.R -f Campbell -m polynomial
Rscript generate_cutoff_UMAP.R -f Campbell


Rscript generate_model.R -f Kolodziejczyk -m linear
Rscript generate_model.R -f Kolodziejczyk -m spline
Rscript generate_model.R -f Kolodziejczyk -m polynomial
Rscript generate_model_UMAP.R -f Kolodziejczyk -m linear
Rscript generate_model_UMAP.R -f Kolodziejczyk -m spline
Rscript generate_model_UMAP.R -f Kolodziejczyk -m polynomial
Rscript generate_cutoff_UMAP.R -f Kolodziejczyk


Rscript generate_model.R -f Lawlor -m linear
Rscript generate_model.R -f Lawlor -m spline
Rscript generate_model.R -f Lawlor -m polynomial
Rscript generate_model_UMAP.R -f Lawlor -m linear
Rscript generate_model_UMAP.R -f Lawlor -m spline
Rscript generate_model_UMAP.R -f Lawlor -m polynomial
Rscript generate_cutoff_UMAP.R -f Lawlor


Rscript generate_model.R -f Lun -m linear
Rscript generate_model.R -f Lun -m spline
Rscript generate_model.R -f Lun -m polynomial
Rscript generate_model_UMAP.R -f Lun -m linear
Rscript generate_model_UMAP.R -f Lun -m spline
Rscript generate_model_UMAP.R -f Lun -m polynomial
Rscript generate_cutoff_UMAP.R -f Lun


Rscript generate_model.R -f Macosko -m linear
Rscript generate_model.R -f Macosko -m spline
Rscript generate_model.R -f Macosko -m polynomial
Rscript generate_model_UMAP.R -f Macosko -m linear
Rscript generate_model_UMAP.R -f Macosko -m spline
Rscript generate_model_UMAP.R -f Macosko -m polynomial
Rscript generate_cutoff_UMAP.R -f Macosko


Rscript generate_model.R -f Messmer -m linear
Rscript generate_model.R -f Messmer -m spline
Rscript generate_model.R -f Messmer -m polynomial
Rscript generate_model_UMAP.R -f Messmer -m linear
Rscript generate_model_UMAP.R -f Messmer -m spline
Rscript generate_model_UMAP.R -f Messmer -m polynomial
Rscript generate_cutoff_UMAP.R -f Messmer


Rscript generate_model.R -f Nestorowa -m linear
Rscript generate_model.R -f Nestorowa -m spline
Rscript generate_model.R -f Nestorowa -m polynomial
Rscript generate_model_UMAP.R -f Nestorowa -m linear
Rscript generate_model_UMAP.R -f Nestorowa -m spline
Rscript generate_model_UMAP.R -f Nestorowa -m polynomial
Rscript generate_cutoff_UMAP.R -f Nestorowa


Rscript generate_model.R -f Richard -m linear
Rscript generate_model.R -f Richard -m spline
Rscript generate_model.R -f Richard -m polynomial
Rscript generate_model_UMAP.R -f Richard -m linear
Rscript generate_model_UMAP.R -f Richard -m spline
Rscript generate_model_UMAP.R -f Richard -m polynomial
Rscript generate_cutoff_UMAP.R -f Richard


Rscript generate_model.R -f Shekhar -m linear
Rscript generate_model.R -f Shekhar -m spline
Rscript generate_model.R -f Shekhar -m polynomial
Rscript generate_model_UMAP.R -f Shekhar -m linear
Rscript generate_model_UMAP.R -f Shekhar -m spline
Rscript generate_model_UMAP.R -f Shekhar -m polynomial
Rscript generate_cutoff_UMAP.R -f Shekhar

Rscript generate_model.R -f Zeisel -m linear
Rscript generate_model.R -f Zeisel -m spline
Rscript generate_model.R -f Zeisel -m polynomial
Rscript generate_model_UMAP.R -f Zeisel -m linear
Rscript generate_model_UMAP.R -f Zeisel -m spline
Rscript generate_model_UMAP.R -f Zeisel -m polynomial
Rscript generate_cutoff_UMAP.R -f Zeisel
