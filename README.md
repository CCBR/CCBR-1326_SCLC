# CCBR 1326

Integrate ATAC-seq and RNA-seq data in small cell lung cancer

[![doi](https://img.shields.io/badge/doi-10.1016%2Fj.xcrm.2026.103027-blue)](https://doi.org/10.1016%2Fj.xcrm.2026.103027)

## Citation

This work was published in the following paper:

> Kumar, R., Desai, P., El Meskini, R., Mohindroo, C., Benő, A.Á., Atkinson, D., Nirula, M., Bassel, L., Takahashi, N., Sovacool, K.L., Meinhardt, A., Arora, M., Parmar, K., Nausome, D., Nichols, S., Sciuto, L., Swift, S., Nemes, K., Sharma, A.K., Pate, N., Rajan, A., Levy, E., Kassin, M., Aldana, C.F., Koparde, V., Galloux, M., Shreshta, R., Ohler, Z.W., Pongor, L., Thomas, A., 2026. A stem-like chromatin program in small-cell lung cancer is associated with poor outcomes after chemoimmunotherapy. Cell Reports Medicine 103027. https://doi.org/10.1016/j.xcrm.2026.103027

```bibtex
@article{kumar_stem-like_2026,
  title = {A Stem-like Chromatin Program in Small-Cell Lung Cancer Is Associated with Poor Outcomes after Chemoimmunotherapy},
  author = {Kumar, Rajesh and Desai, Parth and El Meskini, Rajaa and Mohindroo, Chirayu and Ben{\H o}, Alexandra {\'A}gnes and Atkinson, Devon and Nirula, Michael and Bassel, Laura and Takahashi, Nobuyuki and Sovacool, Kelly L. and Meinhardt, Anna and Arora, Mohit and Parmar, Kanak and Nausome, Darryl and Nichols, Samantha and Sciuto, Linda and Swift, Shannon and Nemes, Kolos and Sharma, Ajit Kumar and Pate, Nathan and Rajan, Arun and Levy, Elliot and Kassin, Michael and Aldana, Christopher Febres and Koparde, Vishal and Galloux, Melissa and Shreshta, Roshan and Ohler, Zoe Weaver and Pongor, Lorinc and Thomas, Anish},
  year = 2026,
  month = sep,
  journal = {Cell Reports Medicine},
  pages = {103027},
  issn = {26663791},
  doi = {10.1016/j.xcrm.2026.103027},
  url = {https://linkinghub.elsevier.com/retrieve/pii/S2666379126004441},
  urldate = {2026-09-14},
  langid = {english}
}
```

## Usage

on Biowulf

### integrate ATAC & RNAseq DGE

```sh
module load nextflow
nextflow run main.nf -resume
```

### tobias TF footprinting

```sh
bash /data/CCBR_Pipeliner/Pipelines/CCBR_tobias/sovacool-dev-tobias/run_tobias.bash \
  -w=/data/SCLCgenomics/sovacoolkl/CCBR-1326/ccbr_tobias \
  -m=dryrun
```
