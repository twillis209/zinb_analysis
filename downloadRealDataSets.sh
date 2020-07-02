#  Allen dataset is available via the scRNASeq R package

# Patel
wget -P real_data/patel https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57872/suppl/GSE57872_GBM_data_matrix.txt.gz

# Zheng 
wget -P real_data/zheng http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz
wget -P real_data/zheng http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_raw_gene_bc_matrices.tar.gz

# Fletcher
wget -P real_data/fletcher https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95601/suppl/GSE95601_oeHBCdiff_Cufflinks_eSet_counts_table.txt.gz 

# Zeisel
wget -P real_data/zeisel https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
wget -P real_data/zeisel https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt

# Kolodziejczyk

wget -P real_data/kolodziejczyk https://espresso.teichlab.sanger.ac.uk/static/counttable_es.csv
wget -P real_data/kolodziejczyk https://espresso.teichlab.sanger.ac.uk/static/counttable_es_norm.csv
