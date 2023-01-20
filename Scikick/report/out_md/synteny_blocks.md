
```r
library(pacman)
p_load(data.table, rtracklayer, ggplot2, ggthemes,
       plyranges, ggpubr, BRGenomics, reshape2, plotly,
       dplyr, gplots, Biostrings, scales)
```


```r
# nolint start
PROJECT <- "/home/daniele/Desktop/IV_course/I_semester/Kursinis_projektas/"
INPUTS <- paste0(PROJECT, "Inputs/")
# nolint end
```


```
[1] "Reading peak file: 1"
[1] "Total peak count: 137787"
[1] "Overlap count: 12211"
[1] "Percentage: 8.86222938303323"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chrX 52988454-52988755      * |        Tbx5      P70326   52988596        19  21.83139
  [2]     chrX 53056855-53057437      * |        Tbx5      P70326   53057322        14   7.74618
  [3]     chrX 53058442-53058874      * |        Tbx5      P70326   53058615        25  30.00492
  [4]     chrX 53099965-53100078      * |        Tbx5      P70326   53100006         6   4.81963
  [5]     chrX 53166272-53166387      * |        Tbx5      P70326   53166321         7   6.08279
  [6]     chrX 53167212-53167393      * |        Tbx5      P70326   53167264         8   7.25168
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]        11.79810  19.03095    137004
  [2]         4.61562   5.47864    137005
  [3]        14.76271  27.00569    137006
  [4]         3.97458   2.77291    137007
  [5]         4.68707   3.93676    137008
  [6]         5.30914   5.01078    137009
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         5         757
  [2]         6         757
  [3]         7         757
  [4]         8         757
  [5]         9         757
  [6]        10         757
  -------
  queryLength: 137787 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 2"
[1] "Total peak count: 163190"
[1] "Overlap count: 17005"
[1] "Percentage: 10.4203688951529"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chrX 53173628-53174120      * |        Tbx5      P70326   53173961        26   17.4045
  [2]     chrX 53179665-53180571      * |        Tbx5      P70326   53180019        87   97.0674
  [3]     chrX 53187933-53189405      * |        Tbx5      P70326   53188230        27   18.4172
  [4]     chrX 53197732-53198122      * |        Tbx5      P70326   53197942        48   42.5784
  [5]     chrX 53199809-53200682      * |        Tbx5      P70326   53200281        23   14.4628
  [6]     chrX 53201636-53202019      * |        Tbx5      P70326   53201805        22   13.5164
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         7.23993   15.2035    161180
  [2]        23.59682   93.4227    161181
  [3]         7.50808   16.1838    161182
  [4]        13.13914   39.7632    161183
  [5]         6.43550   12.3633    161184
  [6]         6.16735   11.4528    161185
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         5         303
  [2]         6         303
  [3]         7         303
  [4]         8         303
  [5]         9         303
  [6]        10         303
  -------
  queryLength: 163190 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 3"
[1] "Total peak count: 135532"
[1] "Overlap count: 14171"
[1] "Percentage: 10.4558333087389"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chrX 53161218-53161673      * |        Tbx5      P70326   53161426        42  28.27330
  [2]     chrX 53166806-53167921      * |        Tbx5      P70326   53167263        29  13.44138
  [3]     chrX 53179637-53180375      * |        Tbx5      P70326   53180011        46  36.73520
  [4]     chrX 53188042-53188437      * |        Tbx5      P70326   53188282        20  13.42739
  [5]     chrX 53197724-53198123      * |        Tbx5      P70326   53197925        27  15.26699
  [6]     chrX 53200132-53200531      * |        Tbx5      P70326   53200375        20   8.44185
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         8.47291  24.99133    133841
  [2]         4.95868  10.90113    133842
  [3]        11.05882  33.11527    133843
  [4]         6.56250  10.88777    133844
  [5]         6.02151  12.61735    133845
  [6]         4.07767   6.23147    133846
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         5         303
  [2]         6         303
  [3]         7         303
  [4]         8         303
  [5]         9         303
  [6]        19         757
  -------
  queryLength: 135532 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 4"
[1] "Total peak count: 117366"
[1] "Overlap count: 12303"
[1] "Percentage: 10.4825929144727"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chrX 53131364-53131773      * |        Tbx5      P70326   53131464        13   5.33706
  [2]     chrX 53142757-53143134      * |        Tbx5      P70326   53142939        12   5.38386
  [3]     chrX 53161180-53161623      * |        Tbx5      P70326   53161471        18   9.62079
  [4]     chrX 53167081-53167880      * |        Tbx5      P70326   53167550        27  15.76371
  [5]     chrX 53187918-53188154      * |        Tbx5      P70326   53187986        10   4.01532
  [6]     chrX 53200319-53200707      * |        Tbx5      P70326   53200470        15   5.65540
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         3.42969   3.35252    115961
  [2]         3.57886   3.37801    115962
  [3]         4.93827   7.24308    115963
  [4]         6.24442  12.95829    115964
  [5]         3.02826   2.16376    115965
  [6]         3.36700   3.62811    115966
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         4         303
  [2]         5         303
  [3]        12         757
  [4]        13         757
  [5]        14         757
  [6]        15         757
  -------
  queryLength: 117366 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 5"
[1] "Total peak count: 66176"
[1] "Overlap count: 6553"
[1] "Percentage: 9.90238152804642"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chrX 53015653-53015791      * |        Tbx5      P70326   53015719        10   4.88396
  [2]     chrX 53045772-53045878      * |        Tbx5      P70326   53045802         9   4.09712
  [3]     chrX 53057374-53057577      * |        Tbx5      P70326   53057534        20   5.98094
  [4]     chrX 53166985-53167468      * |        Tbx5      P70326   53167169        11   4.79155
  [5]     chrX 53179689-53180280      * |        Tbx5      P70326   53180040        34  26.53271
  [6]     chrX 53197609-53198215      * |        Tbx5      P70326   53197951        62  50.19518
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         3.54677   2.59276     65453
  [2]         3.18725   1.92314     65454
  [3]         3.05566   3.57697     65455
  [4]         3.36606   2.52437     65456
  [5]         9.69529  23.05856     65457
  [6]        12.25681  46.27299     65458
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         4         757
  [2]         5         757
  [3]         6         757
  [4]         7         757
  [5]         8         757
  [6]         9         757
  -------
  queryLength: 66176 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 6"
[1] "Total peak count: 13520"
[1] "Overlap count: 1493"
[1] "Percentage: 11.042899408284"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chr9 60371927-60372186      * |        Tbx5      P70326   60372062     15.69  11.05920
  [2]     chr9 60508391-60508734      * |        Tbx5      P70326   60508519     16.19  11.80336
  [3]     chrX 12773243-12773407      * |        Tbx5      P70326   12773312     12.14   6.68330
  [4]     chrX 52474414-52474650      * |        Tbx5      P70326   52474526     15.69   9.31879
  [5]     chrX 52506063-52506272      * |        Tbx5      P70326   52506163     13.66   9.03319
  [6]     chrX 53167424-53167592      * |        Tbx5      P70326   53167469     10.63   6.22160
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         6.63926   7.04380     12968
  [2]         6.65317   7.75516     12969
  [3]         4.38128   3.28892     13403
  [4]         5.56193   5.53531     13429
  [5]         5.83393   5.27094     13430
  [6]         4.62592   2.85409     13432
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1         755
  [2]         5         133
  [3]         6         133
  [4]         7         133
  [5]         7         818
  [6]         8         818
  -------
  queryLength: 13520 / subjectLength: 1403
[1] " "
[1] "Reading peak file: 7"
[1] "Total peak count: 13440"
[1] "Overlap count: 1255"
[1] "Percentage: 9.33779761904762"
GRanges object with 6 ranges and 8 metadata columns:
      seqnames            ranges strand |        name  uniprot_id abs_summit    pileup   p_value
         <Rle>         <IRanges>  <Rle> | <character> <character>  <numeric> <numeric> <numeric>
  [1]     chr9 60526131-60526435      * |        Tbx5      P70326   60526271        28  11.66788
  [2]     chr9 66381745-66382345      * |        Tbx5      P70326   66382182        39  20.46344
  [3]     chrX 11343967-11344297      * |        Tbx5      P70326   11344116        37  20.37194
  [4]     chrX 12880358-12880654      * |        Tbx5      P70326   12880476        36  19.48640
  [5]     chrX 12978490-12978872      * |        Tbx5      P70326   12978674        41  24.02334
  [6]     chrX 13072767-13072977      * |        Tbx5      P70326   13072866        20   7.14124
      fold_enrichment   q_value        ID
            <numeric> <numeric> <numeric>
  [1]         4.42152   8.35448     13069
  [2]         6.09864  16.73770     13122
  [3]         6.37086  16.65122     13384
  [4]         6.20321  15.79929     13385
  [5]         7.04148  20.16832     13386
  [6]         3.52074   4.13549     13387
  -------
  seqinfo: 36 sequences from an unspecified genome
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1         757
  [2]         2           7
  [3]         3         133
  [4]        14        1283
  [5]        15         182
  [6]        15         694
  -------
  queryLength: 13440 / subjectLength: 1403
[1] " "
```


<!-- Use HTML for horizontal lines instead of markdown (---) so
that rmarkdown render does not treat the markdown
as a YAML header -->

<hr></hr>

<details> <summary> Click to see page metadata </summary>

<div>



Computation Started: `2023-01-19 15:35:19`

Finished in `43.119 secs`

<hr></hr>


**Git Log** 



No git history available for this page

---


**Packages**


|package              |version   |date       |
|:--------------------|:---------|:----------|
|utils                |4.2.2     |2022-11-12 |
|bitops               |1.0-7     |2022-04-28 |
|matrixStats          |0.63.0    |2022-12-13 |
|bit64                |4.0.5     |2022-04-29 |
|RColorBrewer         |1.1-3     |2022-05-02 |
|httr                 |1.4.4     |2022-12-13 |
|GenomeInfoDb         |1.32.4    |2022-12-13 |
|tools                |4.2.2     |2022-11-12 |
|backports            |1.4.1     |2022-04-29 |
|utf8                 |1.2.2     |2022-04-28 |
|R6                   |2.5.1     |2022-04-28 |
|KernSmooth           |2.23-20   |NA         |
|lazyeval             |0.2.2     |2022-04-28 |
|DBI                  |1.1.3     |2022-10-07 |
|BiocGenerics         |0.42.0    |2022-04-28 |
|colorspace           |2.0-3     |2022-04-28 |
|withr                |2.5.0     |2022-04-28 |
|tidyselect           |1.2.0     |2022-12-13 |
|DESeq2               |1.36.0    |2022-04-30 |
|git2r                |0.30.1    |2022-05-27 |
|bit                  |4.0.5     |2022-12-13 |
|compiler             |4.2.2     |2022-11-12 |
|cli                  |3.4.1     |2022-10-07 |
|Biobase              |2.56.0    |2022-04-29 |
|pacman               |0.5.1     |2022-04-29 |
|datasets             |4.2.2     |2022-11-12 |
|plotly               |4.10.1    |2022-12-13 |
|DelayedArray         |0.22.0    |2022-04-30 |
|rtracklayer          |1.56.1    |2022-12-13 |
|base                 |4.2.2     |2022-11-12 |
|caTools              |1.18.2    |2022-04-28 |
|scales               |1.2.1     |2022-12-13 |
|genefilter           |1.78.0    |2022-04-30 |
|BRGenomics           |1.8.0     |2022-04-30 |
|digest               |0.6.31    |2022-12-13 |
|stringr              |1.5.0     |2022-12-13 |
|Rsamtools            |2.12.0    |2022-04-30 |
|XVector              |0.36.0    |2022-04-29 |
|htmltools            |0.5.4     |2022-12-13 |
|pkgconfig            |2.0.3     |2022-04-28 |
|MatrixGenerics       |1.8.1     |2022-10-07 |
|fastmap              |1.1.0     |2022-04-28 |
|grDevices            |4.2.2     |2022-11-12 |
|htmlwidgets          |1.5.4     |2022-04-28 |
|rlang                |1.0.6     |2022-10-07 |
|ggthemes             |4.2.4     |2022-04-30 |
|RSQLite              |2.2.19    |2022-12-13 |
|BiocIO               |1.6.0     |2022-04-30 |
|generics             |0.1.3     |2022-10-07 |
|jsonlite             |1.8.4     |2022-12-13 |
|gtools               |3.9.4     |2022-12-13 |
|BiocParallel         |1.30.4    |2022-12-13 |
|dplyr                |1.0.10    |2022-12-13 |
|car                  |3.1-1     |2022-12-13 |
|RCurl                |1.98-1.9  |2022-10-07 |
|magrittr             |2.0.3     |2022-04-28 |
|GenomeInfoDbData     |1.2.8     |2022-04-29 |
|Matrix               |1.5-3     |NA         |
|Rcpp                 |1.0.9     |2022-10-07 |
|munsell              |0.5.0     |2022-04-28 |
|S4Vectors            |0.34.0    |2022-04-28 |
|fansi                |1.0.3     |2022-04-28 |
|abind                |1.4-5     |2022-04-28 |
|lifecycle            |1.0.3     |2022-12-13 |
|stringi              |1.7.8     |2022-10-07 |
|yaml                 |2.3.6     |2022-12-13 |
|carData              |3.0-5     |2022-04-30 |
|SummarizedExperiment |1.26.1    |2022-12-13 |
|zlibbioc             |1.42.0    |2022-04-29 |
|gplots               |3.1.3     |2022-04-28 |
|plyr                 |1.8.8     |2022-12-13 |
|grid                 |4.2.2     |2022-11-12 |
|blob                 |1.2.3     |2022-04-29 |
|parallel             |4.2.2     |2022-11-12 |
|crayon               |1.5.2     |2022-10-07 |
|methods              |4.2.2     |2022-11-12 |
|lattice              |0.20-45   |NA         |
|Biostrings           |2.64.1    |2022-12-13 |
|splines              |4.2.2     |2022-11-12 |
|annotate             |1.74.0    |2022-04-30 |
|KEGGREST             |1.36.3    |2022-12-13 |
|locfit               |1.5-9.6   |2022-10-07 |
|knitr                |1.41      |2022-12-13 |
|pillar               |1.8.1     |2022-12-13 |
|ggpubr               |0.5.0     |2022-12-13 |
|GenomicRanges        |1.48.0    |2022-04-30 |
|rjson                |0.2.21    |2022-04-30 |
|ggsignif             |0.6.4     |2022-12-13 |
|geneplotter          |1.74.0    |2022-04-30 |
|reshape2             |1.4.4     |2022-04-28 |
|codetools            |0.2-18    |NA         |
|stats4               |4.2.2     |2022-11-12 |
|XML                  |3.99-0.13 |2022-12-13 |
|glue                 |1.6.2     |2022-04-28 |
|evaluate             |0.19      |2022-12-14 |
|data.table           |1.14.6    |2022-12-13 |
|vctrs                |0.5.1     |2022-12-13 |
|png                  |0.1-8     |2022-12-13 |
|graphics             |4.2.2     |2022-11-12 |
|gtable               |0.3.1     |2022-10-07 |
|purrr                |0.3.5     |2022-12-13 |
|tidyr                |1.2.1     |2022-12-13 |
|assertthat           |0.2.1     |2022-04-28 |
|cachem               |1.0.6     |2022-04-28 |
|ggplot2              |3.4.0     |2022-12-13 |
|xfun                 |0.35      |2022-12-13 |
|xtable               |1.8-4     |2022-04-28 |
|broom                |1.0.1     |2022-12-13 |
|restfulr             |0.0.15    |2022-12-13 |
|rstatix              |0.7.1     |2022-12-13 |
|viridisLite          |0.4.1     |2022-10-07 |
|survival             |3.4-0     |NA         |
|tibble               |3.1.8     |2022-12-13 |
|GenomicAlignments    |1.32.1    |2022-12-13 |
|stats                |4.2.2     |2022-11-12 |
|AnnotationDbi        |1.58.0    |2022-04-29 |
|plyranges            |1.16.0    |2022-04-30 |
|memoise              |2.0.1     |2022-04-29 |
|IRanges              |2.30.1    |2022-10-07 |

---

**System Information**


|         |systemInfo                                  |
|:--------|:-------------------------------------------|
|version  |R version 4.2.2 Patched (2022-11-10 r83330) |
|platform |x86_64-pc-linux-gnu (64-bit)                |
|locale   |lt_LT.utf8                                  |
|OS       |Linux Mint 20.3                             |
|UI       |X11                                         |

**Scikick Configuration**


```bash
cat scikick.yml
```

```
### Scikick Project Workflow Configuration File

# Directory where Scikick will store all standard notebook outputs
reportdir: report

# --- Content below here is best modified by using the Scikick CLI ---

# Notebook Execution Configuration (format summarized below)
# analysis:
#  first_notebook.Rmd:
#  second_notebook.Rmd: 
#  - first_notebook.Rmd 	# must execute before second_notebook.Rmd
#  - functions.R 	        # file is used by second_notebook.Rmd
#
# Each analysis item is executed to generate md and html files, E.g.:
# 1. <reportdir>/out_md/first_notebook.md
# 2. <reportdir>/out_html/first_notebook.html
analysis: !!omap
- synteny_blocks.Rmd:
version_info:
  snakemake: 7.0.0
  ruamel.yaml: 0.17.21
  scikick: 0.2.0
```

---

**Functions**





</div>
</details>

<!-- Project map inserted below after knitting -->

