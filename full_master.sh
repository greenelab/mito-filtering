## Scripts:
##	prepare_our_data.R: takes the raw matrix of gene by cell counts returned
##	by alevin, and prepares it as a SingleCellExperiment object with
##	basic QC run and percent mitochondrial reads calculated.
##	generate_model.R: runs a mixture model on the standardized data, generates
##	some explanatory plots, and outputs a keep/don't keep recommendation per cell.
## Parameters:
##	-f, --file <sample id>: which tumor you want to run the given script on.
##	-m, --model: what kind of mixture model you want to use to filter. Current
##	             options are linear, spline, and polynomial.

Rscript prepare_our_data.R -f 16030X2
Rscript prepare_our_data.R -f 16030X3
Rscript prepare_our_data.R -f 16030X4

Rscript generate_model.R -f 16030X2 -m linear
Rscript generate_model.R -f 16030X2 -m spline
Rscript generate_model.R -f 16030X2 -m polynomial

Rscript generate_model.R -f 16030X3 -m linear
Rscript generate_model.R -f 16030X3 -m spline
Rscript generate_model.R -f 16030X3 -m polynomial

Rscript generate_model.R -f 16030X4 -m linear
Rscript generate_model.R -f 16030X4 -m spline
Rscript generate_model.R -f 16030X4 -m polynomial

