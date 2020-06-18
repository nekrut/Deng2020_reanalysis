# Analysis of [Deng 2020](https://science.sciencemag.org/content/early/2020/06/05/science.abb9263) data

A recent publication by Deng et al. [(2020)](https://science.sciencemag.org/content/early/2020/06/05/science.abb9263) analyzed transmission patters of SARS-CoV2 in Northern California. This is one of very rare SARS-CoV2 publications in that it makes primary data (sequencing reads) available (BioProject [PRJNA629889](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA629889)).     

We are a part of the GalaxyProject [effort](https://covid19.galaxyproject.org) aimed at developing open and freely accessible practices for the computational analyses of SARS-CoV2 data. As a part of this effort we already automatically downloaded and processed these data. Naturally, we decided to look if the results of our pipeline agree with the results reported by Deng et al. Here is what we discovered. 

## Data and associated metadata discrepancy: 36 patients, 44 SRA datasets, and 66 samples.

The paper states that the authors analyzed 36 infected patents. The corresponding short read archive data corresponding to [PRJNA629889](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA629889) lists 43 samples with the following metadata (an example representing one sample):

|    | study_accession   | experiment_accession   | experiment_title                                | experiment_desc                                 |   organism_taxid  | organism_name                                   | library_strategy   | library_source   | library_selection   | sample_accession   |   sample_title | instrument            |   total_spots |   total_size | run_accession   |   run_total_spots |   run_total_bases | run_alias     | sra_url_alt1                                                              | sra_url_alt2                                                              | sra_url                                                                   |   experiment_alias |   Titer (Ct value) | Sequencing Platform   |   Distinguishing Info | strain         | isolate        | collected_by                           | collection_date   | geo_loc_name   | host         | host_disease   | isolation_source   | lat_lon              | BioSampleModel   | sra_url_alt                                                                             | ena_fastq_url                                                                   | ena_fastq_ftp                                                                      |
|---:|:------------------|:-----------------------|:------------------------------------------------|:------------------------------------------------|------------------:|:------------------------------------------------|:-------------------|:-----------------|:--------------------|:-------------------|---------------:|:----------------------|--------------:|-------------:|:----------------|------------------:|------------------:|:--------------|:--------------------------------------------------------------------------|:--------------------------------------------------------------------------|:--------------------------------------------------------------------------|-------------------:|-------------------:|:----------------------|----------------------:|:---------------|:---------------|:---------------------------------------|:------------------|:---------------|:-------------|:---------------|:-------------------|:---------------------|:-----------------|:----------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|
|  0 | SRP265005         | SRX8409213             | Severe acute respiratory syndrome coronavirus 2 | Severe acute respiratory syndrome coronavirus 2 |           2697049 | Severe acute respiratory syndrome coronavirus 2 | AMPLICON           | METAGENOMIC      | PCR                 | SRS6721570         |            nan | Illumina HiSeq 1500   |         13478 |      1016143 | SRR11859166     |             13478 |           1886920 | UC13.fastq.gz | https://storage.googleapis.com/sra-pub-src-9/SRR11859166/UC13.fastq.gz.1  | https://sra-pub-src-9.s3.amazonaws.com/SRR11859166/UC13.fastq.gz.1        | https://sra-download.ncbi.nlm.nih.gov/traces/sra74/SRR/011581/SRR11859166 |                nan |               27.5 | HiSeq                 |                    13 | Not applicable | Not applicable | University of Califronia San Francisco | 2020-03-02        | USA:California | Homo sapiens | Acue infection | clinical sample    | 36.7783 N 119.4179 W | Pathogen.cl      | nan                                                                                     | nan                                                                             | nan                                                                                |

One can see that many useful columns are either empty (e.g., `sample_title`) or do not contain useful information (e.g., `geo_loc_name`). This is unfortunate as their supplemental Table S1 contains Sample Lab IDs and it is impossible to figure out which is which. In addition, Table S1 lists 66 samples. So there is a discrepancy: 36 patients, 44 SRA datasets, and 66 samples.

Downloading data resulted in 25 single end and 18 paired end datasets. We performed initial QC using [`fastp`](https://github.com/OpenGene/fastp), to remove low quality reads as well as to detect and trim Illumina adapter sequences:

![](fastp_qc.png)




