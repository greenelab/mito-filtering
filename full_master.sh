## Scripts:
##	prepare_our_data.R: takes the raw matrix of gene by cell counts returned
##		by alevin, and prepares it as a SingleCellExperiment object with
##		basic QC run and percent mitochondrial reads calculated.
##	generate_model.R: runs a mixture model on the standardized data, generates
##		some explanatory plots, and outputs a keep/don't keep filter per cell.
##	generate_model_UMAP.R: create UMAP of data before and after smart filtering.
##	generate_cutoff_UMAP.R: create UMAP of data before and after simple cutoff.
##	generate_two_way_table.R: counts number of cells excluded and included by 
##		each method; we're particularly interested in discordant decisions.
## Parameters:
##	-f, --file <16030X2>: which tumor you want to run the given script on.
##	-m, --model <linear>: what kind of mixture model you want to use to filter. Current
##		options are linear, spline, and polynomial.
##	-p, --percent <20>: what percent of mitochondrial reads to use as a simple cutoff
##		for comparison analyses.

#Rscript -e 'rmarkdown::render("prepare_public_data.Rmd")'
#Rscript prepare_our_data.R -f 16030X2
#Rscript prepare_our_data.R -f 16030X3
#Rscript prepare_our_data.R -f 16030X4

samples=('16030X2' '16030X3' '16030X4' 'Bach' 'Buettner' 'Campbell' 'Kolodziejczyk' 'Lawlor' 'Macosko' 'Messmer' 'Richard' 'Shekhar' 'Zeisel')
models=('linear' 'spline' 'polynomial')
cutoffs=(10 15 20 25)

parallel -j3 Rscript generate_model.R -m {1} -f {2} ::: ${models[*]} ::: ${samples[*]}
#parallel -j3 Rscript generate_model_UMAP.R -m {1} -f {2} ::: ${models[*]} ::: ${samples[*]}
#parallel -j3 Rscript generate_cutoff_UMAP.R -p {1} -f {2} ::: ${cutoffs[*]} ::: ${samples[*]}
#parallel -j3 Rscript generate_two_way_table.R -p {1} -m {2} -f {3} ::: ${cutoffs[*]} ::: ${models[*]} ::: ${samples[*]}
