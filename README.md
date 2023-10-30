# Mitogenome-extraction

## <b>Mitogenome-extraction</b> is a Nextflow workflow that assembles and annotates mitogenomes from Illumina short reads.

This pipeline is specialized to retrieve the mitogenome from low-quality assemblies, which may result in the mitogenome being either distributed among multiple contigs and/or assembled with low coverage. If the mitogenome retrieval was successful, all identified genes will be annotated and the results streamlined. 

### Input for this pipeline are:
- the reads obtained through whole-genome-sequencing (WGS)
- the corresponding assembled contigs
- a reference mitogenome from a related species in FASTA format


### Mitogenome-extraction will attempt to:
1) Identify all contigs containing sequences belonging to the mitochondrion (Method: Blastn).

2) Recover the mitogenome as a single contig. If the mitogenome is split onto multiple contigs, attempt reassembly (Method: NOVOPlasty).
   
3) Ensure the strand direction of the assembly is the same as in the reference genome.

4) Annotate the mitogenome (Method: MITOS). Save each gene as a separate file with the header `>identifier@gene`. Adjust gene order to list cox1 first.

5) Obtain barcoding sequence ( --barcode <reference_sequence>.fasta ).
-----
## Installation

### Installing Nextflow and Anaconda

[Nextflow](https://www.nextflow.io/) and [Anaconda](https://www.anaconda.com/) are required to run this pipeline. For the
installation, simply follow the [Nextflow installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation)
and [Anaconda installation instructions](https://docs.anaconda.com/anaconda/install/index.html) in their respective documentation.

### Installing better-fasta-grep

better-fasta-grep is required to run this pipeline. Available sources are [pip](https://pypi.org/project/better-fasta-grep/) and [gitlab](https://gitlab.com/fethalen/bfg).

### Downloading this Pipeline

To download this pipeline, go to the folder where you wish to store this
pipeline in and then run the following:

```
$ git clone https://github.com/ThiloSchulze/mitogenome-extraction
```

## Parameter list

```
usage: nextflow run main.nf --reads '<reads>.fasta' --contigs <contigs>.fasta --mitogenome <reference_Mito>.fasta 
       [--help] [--species_id '<species_id>'] [--output '<dir>'] [--min_blast_wordsize <size>] [--max_blast_wordsize <size>]
       [--nucleotide_size <size>] [--min_size <size>] [--max_size <size>] [--kmer_size <size>] [--read_length <size>]
       [--insert_size <size>] [--genetic_code '<genetic code>'] [--mitos_reference '<dir>']

required arguments:

  --reads '<reads>.fasta'           Illumina short reads from your species
  --contigs <contigs>.fasta         Assembled fasta contigs to be searched to find the mitogenome
  --mitogenome <relatedMito>.fasta  Reference mitogenome from a closely related species
  
optional arguments:

   //General
  --help                            Print help menu and exit
  --species_id '<species_id>'       Identifier used as header for annotated genes, e.g. >species_id@cox1
  --output '<dir>'                  Specify output directory (default: 'mitogenome-extraction')
  
  //Blastn
  --min_blast_wordsize <size>       Minimum word size word size used for blastn search (default: 11)
  --max_blast_wordsize <size>       Maximum word size word size used for blastn search (default: 25)
  --nucleotide_size <size>          Estimated nucleotide count of target mitochondrion (default: 16500)
  
  //NOVOPlasty (USAGE: nextflow run main.nf)
   --min_size <size>                 Estimated minimum nucleotide count of target mitochondrion (default: 12000)
   --max_size <size>                 Estimated maximum nucleotide count of target mitochondrion (default: 22000)
   --kmer_size <size>                K-mer size used for reassembly (default: 33)
   --read_length <size>              Read length of Illumina short reads (default: 151)
   --insert_size <size>              Insert size of paired-end reads (default: 300)
   
  //PRICE      (USAGE: nextflow run price.nf)
   To be added soon.
  
  // Genetic code (default: Invertebrates)
   --genetic_code '<genetic code>'  Select the genetic code appropriate for the target mitochondrion (default: '05')
                                    (NCBI genetic code table numeration: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
   --mitos_reference '<dir>'        Select a reference for gene annotation (default: 'refseq63m')
 ```


### Required arguments

#### `--reads`

* Paired-end reads in FASTQ format are required.
* Forward and reverse reads should be located in the same folder.
* The following file extensions are supported: `.fastq.gz`, `.fq.gz`, `.fastq` and `.fq`.

#### `--contigs`

* Contigs need to be ptovided as a single file in FASTA format.

#### `--mitogenome`

* A mitochondria genome of a closely-related species in FASTA format is required.
* Default is set for the annotation of mitochondria within invertebrates. To switch to a different organism group, the `--genetic_code` and `--mitos_reference` flags need to be modified.

### Example usage & output files

As an example, we have a designated directory for our sample `4Alsu_Alitta_succinea`, which contains the forward and reverse reads `NG-29255_4AlSu_lib572143_7903_2_1.fastq.gz` & `NG-29255_4AlSu_lib572143_7903_2_2.fastq.gz` as well as the contigs file `final_contigs.fasta`. The reference mitogenome is the file `platynereis_dumerilii_complete_mitogenome.fna`.

In a directory structure like this:

```
mitogenome-extraction
├──platynereis_dumerilii_complete_mitogenome.fna
└──4Alsu_Alitta_succinea
   ├── final_contigs.fasta
   └── raw_reads
       ├── NG-29255_4AlSu_lib572143_7903_2_1.fastq.gz
       └── NG-29255_4AlSu_lib572143_7903_2_2.fastq.gz
```

The pipeline could be started with the following command:

```
nextflow run main.nf --mitogenome platynereis_dumerilii_complete_mitogenome.fna --contigs 4Alsu_Alitta_succinea/final_contigs.fasta --reads '4Alsu_Alitta_succinea/raw_reads/NG-29255_4AlSu_lib572143_7903_2_{1,2}.fastq.gz' --species_id '4Alsu_Alitta_succinea' --output '4Alsu_out'
```

The following files are created:

```
4Alsu_out
├── mitogenome_extraction
│   └── single_contig_mitogenome.fa
│   └── stats.txt
│   └── unique_mito_seqid.txt
├── strand_control
│   └── single_contig_mitogenome.fa
└── MITOS_annotation
    ├── individual_genes_nuc
    │   └── cox1.fna
    │   └── cob.fna
    │   └── ...
    └── individual_genes_prot
    │   └── cox1.faa
    │   └── cob.faa
    │   └── ...
    ├── mitos_output
    │   └──adjusted_result.geneorder
    │   └──[MITOS files]
    └── mitos_output.txt
    
```


##### `/mitogenome_extraction/`

 - `stats.txt` may hold relevant additional information if the pipeline fails. It contains general information on the nucleotide count match for each of the 14 blastn searches. The nucleotide count of the recovered contig is listed as 'single_contig_mitogenome.fa'. If a reassembly with NOVOPlasty / PRICE was attempted, the input file is listed as 'split_mitogenome.fa'.

 - `single_contig_mitogenome.fa` contains the retrieved mitochondrial genome (if the recovery was successful).

 - `unique_mito_seqid.txt` lists all original contigs used for recovery of the mitogenome.

##### `/strand_control/`
 - `single_contig_mitogenome.fa` contains the retrieved mitochondrial genome. Strand direction is now matching the reference mitogenome.

##### `/MITOS_annotation/`

 - `individual_genes_nuc/` contains each annotated gene as a file in fasta format. Fasta headers are specified with the `--species_id` flag, e.g.: `>species_id@cox1`. If the `--species_id` is not set, fasta header ids will be taken from the forward read file. All genes listed in a single file are provided as `result.fas`

 - `individual_genes_prot/` is the same as `individual_genes_nuc/`, but contains protein coding gene files.

 - `adjusted_result.geneorder` will be created if cox1 was identified. It contains a shifted gene order with cox1 listed first.

 - `mitos_output.txt` contains MITOS comments on the annotation results. Genes that are missing, duplicated, or interrupted by stop-codons will be listed here.


### Options if reassembly is incomplete / fails

#### Mitogenome extraction method
14 Blastn searches are performed ranging from word size 11 to 25 (Blastn). For each blastn, all contigs with matches to the reference mitogenome are concatenated into a unique file. Contigs containing mitochondrial genomes usually have high coverage. To account for this, all contigs below a certain coverage count are filtered out (default: contigs with coverage below 60). After filtering, the file whose combined nucleotide count is closest to the target mitochondrion size (default: 16500) is selected. The number of contigs in this file is counted. If the file contains 1 contig, perform annotation. If no or only very few contigs match the reference mitogenome, restart with a lower coverage cutoff. If the file contains 2+ contigs, perform reassembly with NOVOPlasty. If the first NOVOPlasty reassembly fails, try again with the file containing the highest average nucleotide count per contig instead. Afterward, perform annotation with the mitogenome, or if the reassembly is incomplete, with the largest reassembly contig.
 
#### Options to improve reassembly
Based on the extraction process described above, there are 3 ways to potentially improve the reassembly. It is advisable to check `output_dir/mitogenome_extraction/stats.txt` before estimating what may need to change.

##### 1) Changing coverage cutoff 
If the majority of blastn search files match a nucleotide count lower than your desired mitochondrion size it might be helpful to lower the coverage cutoff. E.g.: `--coverage_cutoff 3` will lower the cutoff from 60 to 30. Conversely, if your reassembly fails and all blast result files contain too many nucleotides it might be helpful to raise the cutoff to the maximum of 90.

##### 2) Lowering minimum blast word size
If lowering the coverage cutoff does not remove the issue of matching too few contigs, or some genes are still not annotated, lowering the bar required for contig matches might help, e.g.: `--min_blast_wordsize 8`. Be careful with these options, as lowering both the coverage cutoff and the minimum blast word size increases the risk of matching contigs not belonging to the mitochondrial genome.

##### 3) Changing the target nucleotide size
Increasing ( `--nucleotide_size 18000`) or decreasing ( `--nucleotide_size 15000`) the target nucleotide size provides an easy way to direct on which blast result file the reassembly/annotation is performed with.

 
#### Options if NOVOPlasty reassembly fails
If the NOVOPlasty reassembly fails after including options to improve reassembly (e.g. open `/output_dir/mitogenome_extraction/stats.txt` and check if the sum_len of `single_contig_mitogenome.fa` is listed as below 10000) it is recommended to try reassembly with [PRICE](https://pubmed.ncbi.nlm.nih.gov/23550143/) instead. Unlike NOVOPlasty, PRICE is a general reassembler not restricted to assembling organelles. It is, therefore, more robust to the seed file containing non-mitochondrion sequences. Warning: Price will take a lot longer to run than NOVOPlasty. PRICE needs to be [downloaded](https://sourceforge.net/projects/pricedenovo/) and extracted. PriceTI then needs to be placed inside the bin of the base directory of the pipeline: `mitogenome-extraction/bin/PriceTI`. The file extension of the reads needs to be either `.fastq` or `.fastq.gz`. To run PRICE instead of NOVOPlasty, start with the input: 
```
Nextflow run price.nf [commands]
```
