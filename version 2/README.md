# truQuant version 2
truQuant annotation creates PRO-Seq driven nucleotide resolution annotations of the hg38 genome using [tsrFinder](https://github.com/P-TEFb/tsrFinder). 
The script works by first blacklisting the data for non-Pol II RNAs (version 2 includes LSU and SSU RNAs!), 
then maps transcription start regions to genes within a search region from 1,000 bp upstream of GENCODE V31 annotations to the most upstream AUG. 
The new 5’ end of the gene is labeled as the max transcription start site in the max transcription start region. 
Pause regions are defined as ±75 bp from the weighted center of the TSR (avgTSS from tsrFinder) and all other transcription start regions are expanded and blacklisted. 
Gene bodies are defined as the end of the pause region to the 3’ end of the GENCODE V31 annotation. Pause regions are quantified using the sum of 5’ ends in the pause region and 
gene bodies are quantified using 3’ end reads. Importantly, TSRs inside gene bodies are blacklisted.

truQuant can be run using python3 from the command line using the following syntax

```
python3 truQuant.py <Sequencing File For Annotation> <Sequencing Files (Optional)>
```

Optional parameters found in version 1 can be changed by modifying the source code, however we do not recommend changing the values. 

truQuant will make five files: a tsrFinder file, three regions files, and one master file. The three regions files provide genomic coordinates 
for the pause regions, gene bodies, and blacklisted regions. The master file will provide the coordinates of gene, 
statistics on the TSR generating the pause region, and the quantification of the provided datasets. 