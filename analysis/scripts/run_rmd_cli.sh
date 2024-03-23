RMD_FILE="/workspace/analysis/rmd/Seurat_v_Seurat.Rmd"  # "/workspace/analysis/test_cli.Rmd"

for seu2_read_fraction_after_downsampling in 0_16 0_08 0_04 0_02 0_01; do
    for seu2_read_downsample_seed in 100 101 102; do
        for data_input in default seu1; do
        	# Use Rscript to call rmarkdown::render() with parameters
        	Rscript -e "rmarkdown::render('${RMD_FILE}', output_file='/workspace/analysis/output/rmd_html/seu_read_input_${data_input}_frac${seu2_read_fraction_after_downsampling}_seed${seu2_read_downsample_seed}.html', params=list(seu2_read_fraction_after_downsampling='${seu2_read_fraction_after_downsampling}', seu2_read_downsample_seed='${seu2_read_downsample_seed}', data_input='${data_input}'))"
        done
    done
done