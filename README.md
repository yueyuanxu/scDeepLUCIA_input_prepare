### 25bp-resolution epigenomic feature npy file (mm10, astrocyte, chr10)
```bash
python -c "import sys, gzip, numpy, pyBigWig, itertools; from pathlib import Path; print('Package check passed')"
python ./python_scripts/build_tiled_epigenomic_feature.py ./genomic_bin/mm10/genomic_bin_with_mark.bed.gz 5000 25 ./pseudobulk_ATAC_bigwig/Astrocyte_sig.pval.signal.bigwig epigenome_pval_npy/mm10/astrocyte/R1.5kb.npy
python ./python_scripts/build_sliced_multimark_epigenomic_feature.py ./genomic_bin/mm10/genomic_bin_with_mark.bed.gz chr10 epigenome_pval_npy/mm10/astrocyte R1 sliced_epi_array/mm10/astrocyte/R1.5kb.chr10.npy
```

### 5kb-resolution contact matrix npy file
```bash
python -c "import sys, numpy, cooler; from scipy.signal import convolve2d; from pathlib import Path; print('Package check passed')"
python ./python_scripts/cool_to_npy_with_convolution.py ./pseudobulk_contact_cool/Astrocyte.5000.cool chr10 5 sliced_con_array/mm10/astrocyte/convolved_window5.chr10.npy
```
### one-hot encoded genomic feature npy file
./sliced_seq_array/
