# How to run
1. Snakemake to align Fastqs with cellranger --> count matrices (you can skip this step if you downloaded the count matrices directly)
2. Run through `notebooks/Prepare_Data.ipynb` using jupyter notebook. (you can also skip this step if you aquired the file `exp2_lung_tcells.h5`)
3. Produce the figures in the paper with `notebooks/Produce_Figures_4IJK_S15.ipynb` using jupyter notebook.

# exp2_lung_tcells
You can find this file as an rds file here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200596.
You will need to convert it to h5 if you want to use it in scanpy. Check out https://github.com/mojaveazure/seurat-disk as a tool for that.
