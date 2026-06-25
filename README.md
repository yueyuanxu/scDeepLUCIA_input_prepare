## Generate input data for scDeepLUCIA loop prediction

### 25bp-resolution epigenomic feature npy file
Before running, generate pseudo-bulk ATAC bigWig files by following [scHiCAR downstream analysis codes](https://github.com/yueyuanxu/scHiCAR/edit/dev/5_downstream_analysis/README.md?plain=1#L101)
```bash
python -c "import sys, gzip, numpy, pyBigWig, itertools, cooler, pysam; from scipy.signal import convolve2d; from pathlib import Path; print('Package check passed')"
python ./python_scripts/build_tiled_epigenomic_feature.py ./genomic_bin/mm10_genomic_bin_with_mark.bed.gz 5000 25 ./pseudobulk_ATAC_bigwig/Astrocyte_sig.pval.signal.bigwig epigenome_pval_npy/mm10/astrocyte/R1.5kb.npy
for i in chr{1..19} chrX
do
echo "Generating R1.5kb.${i}.npy..."
python ./python_scripts/build_sliced_multimark_epigenomic_feature.py ./genomic_bin/mm10_genomic_bin_with_mark.bed.gz $i epigenome_pval_npy/mm10/astrocyte R1 sliced_epi_array/mm10/astrocyte/R1.5kb.${i}.npy
done
```

### 5kb-resolution contact matrix npy file
```bash
for i in chr{1..19} chrX
do
echo "Generating convolved_window5.${i}.npy..."
python ./python_scripts/cool_to_npy_with_convolution.py ./pseudobulk_contact_cool/Astrocyte.5000.cool $i 5 sliced_con_array/mm10/astrocyte/convolved_window5.${i}.npy
done
```

### one-hot encoded genomic feature npy file
```bash
python3 ./python_scripts/build_genomic_feature.py ./genomic_bin/mm10_genomic_bin_with_mark.bed.gz 5000 mm10.genome.fa genome_onehot_npy/mm10/seq_onehot.5kb.npy
for i in chr{1..19} chrX
do
echo "Generating seq_onehot.5kb.${i}.npy..."
python3 ./python_scripts/build_sliced_genomic_feature.py ./anchor_list/mm10_anchor_list.txt.gz $i genome_onehot_npy/mm10/seq_onehot.5kb.npy sliced_seq_array/mm10/seq_onehot.5kb.${i}.npy
done
```
