# truQuant
truQuant annotation creates PRO-Seq driven nucleotide resolution annotations of the hg38 genome using [tsrFinder](https://github.com/P-TEFb/tsrFinder). The script works by mapping transcription start regions to genes within a search region from 1,000 bp upstream of GENCODE V31 annotations to the most upstream AUG. The new 5’ end of the gene is labeled as the max transcription start site in the max transcription start region. Pause regions are defined as ±75 bp from the 5’ end and all other transcription start regions are expanded and blacklisted. Gene bodies are defined as the end of the pause region to the 3’ end of the GENCODE V31 annotation. truQuant quantitation sums the 5’ reads in the pause region, blacklists 3’ reads inside blacklisted regions, and sums the 3’ reads in the gene bodies to output two files. 

truQuant annotation can be run using python3 from the command line using the following syntax

```
python3 truQuant_annotation.py –t [TSR file.bed] –s [Sequencing file.bed] –a [Annotation file ] [optional parameters]
```

Optional parameters:
```
-n [int] = TSR window size
-e [int] = Annotation extension
-b [int] = Blacklist extension amount
-p [int 0-1] = Minimum percentage for blacklisting TSRs
```

truQuant annotation will output four files: an annotation file, a pause region file, a gene bodies region file, and a blacklist region file. These files are used for quantitation part of truQuant.

truQuant quantitation can be run using python3 from the command line using the following syntax

```
python3 truQuant_quantitation.py –p [Pause region file.bed] –g [Gene bodies file.bed] –s [Sequencing file.bed] –b [Blacklist file.bed]
```
