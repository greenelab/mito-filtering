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

