# Code repository for eukrhythmic paper

Code and data repository for 2022 eukrhythmic manuscript.

## Snakemake workflows

All `Snakemake` workflows used to further process output from `eukrhythmic` are stored in the `snakemake-workflows` directory. These include:

1. `blast-snake`: conduct a `BLAST` search of the output proteins from the `eukrhythmic` assembly against the original proteins of the "designer" assembly. This is to safeguard against the scenario in which contigs don't cluster together but have suitable sequence-based match results.
2. `busco-snake`: assess the `BUSCO` completeness of each of the MMETSP transcriptomes that were included in the analysis to build the designer metatranscriptomes.
3. `clustering-snake`: use `mmseqs2` to cluster the contigs from the designer transcriptomes at some level of coverage and percentage sequence identity.
4. `eggnog-snake`: use `eggNOG-mapper` to assign functional annotations to protein outputs from assembly and original designer proteins.
5. `eukulele-snake`: use `EUKulele` to assign a taxonomic annotation to contigs from the output assemblies as well as determine how many of the original designer contigs can be assigned an annotation via the same method.
6. `process-four-way-contigs-snake`: process contigs that just are in common between _n_ assemblers (e.g. four, as the directory is named), as though that is the only eukrhythmic output that we're using (cut out anything that isn't agreed upon between some >1 number of assemblers). 
7. `salmon-snake`: calculate coverage of raw reads by assemblies using the `Salmon` software.
8. `sourmash-snake`: calculate similarity within and between communities using the `sourmash` software.
9. `transdecoder-snake`: identify proteins using `TransDecoder`.