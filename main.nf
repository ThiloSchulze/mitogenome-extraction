#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

ch_contigs = Channel
  .fromPath( params.contigs, type: 'file' )
ch_mitogenome = Channel
  .fromPath( params.mitogenome, type: 'file' )
ch_rawReads = Channel
  .fromFilePairs( params.reads, size : 2, type: 'file' )
  .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
  .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.reads}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'\n" +
             "For single-end reads, specify '--single-end'" }

def helpMessage() {
    log.info"""
    Example of command to run:
    nextflow run main.nf --contigs 9GyTe_Gymnonereis_tenera/final_contigs.fasta --reads '9GyTe_Gymnonereis_tenera/raw_reads/NG-29255_9GyTe_lib572147_7903_2_{1,2}.fastq.gz' --mitogenome platynereis_dumerilii_complete_mito.fna --species_id 9_Gymnonereis_tenera --mitos_reference testfolder/  -resume
    """.stripIndent()
}

process extract_mitogenome {
    publishDir "${params.output}/mitogenome_extraction", mode: 'copy'
//    label 'process_low'

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    path(contigs)
    path(mitogenome)

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the host mitogenome
    path('mitogenome_candidates*')
    path('possible_mitogenomes.fa')
    path('mitogenome.fa'), emit: mitogenome optional true
    path('split_mitogenome.fa'), emit: split_mitgenome optional true
    path('warning.txt'), emit: no_mitogenome_match optional true    
    path('stats.txt')
    path('*.txt')
    path('*.fa')

    script:
    """
    touch mitogenome.fa
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa
    cat $contigs | bfg "cov_[6-9][0-9]{1,}\\.[0-9]+" > cov_60_to_99.fa
    cat $contigs | bfg "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > cov_100_plus.fa
    cat cov_60_to_99.fa cov_100_plus.fa > possible_mitogenomes.fa
    makeblastdb -in $contigs -title contig -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {17,25}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query $mitogenome -db db -outfmt "10 sseqid" -word_size \$i -num_threads ${task.cpus} > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        echo "made seqids unique"
        cat possible_mitogenomes.fa | bfg -f unique_seqid.txt > "mitogenome_candidates_wordsize_\$i.fa"
      done

    for file in *candidate*
    do
      grep -v  '^>' \$file | wc -m
    done > nucleotide_count.txt
    closest_match=\$( awk -v c=1 -v t=$params.nucleotide_size 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count.txt )
    echo "\$closest_match" > closest_match.txt
    for blast_result in *candidate*
    do
      if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
      then
        cat \$blast_result > identified_mitogenome.fa
        break
      fi
    done

    grep '^>' identified_mitogenome.fa | cat -n | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_mito_seqid.txt

    if [[ \$(wc -l unique_mito_seqid.txt) = "0 unique_mito_seqid.txt" ]]
    then
      echo "The mitogenome was not identified! It may have a coverage below 60." > warning.txt
    elif [[ \$(wc -l unique_mito_seqid.txt) = "1 unique_mito_seqid.txt" ]]
    then
      cat identified_mitogenome.fa > mitogenome.fa
    else
      cat identified_mitogenome.fa > split_mitogenome.fa
    fi

    seqkit stats *.fa > stats.txt
    """
}

process reassemble_mitogenome {
    publishDir "${params.output}/NOVOPlasty_reassembly", mode: 'copy'
    label 'process_low'
//    beforeScript 'ulimit -s unlimited'

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), path(rawreads)
    path(mitogenome_contigs)

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the host mitogenome
    path('*')
    path('assembled_mitogenome.fasta'), emit: assembled_mitogenome

    script:
    """
    echo "Project:
    -----------------------
    Project name          = Mitogenome
    Type                  = mito
    Genome Range          = $params.min_size-$params.max_size
    K-mer                 = $params.kmer_size
    Max memory            = ${task.memory.toGiga()}
    Extended log          = 0
    Save assembled reads  = no
    Seed Input            = $mitogenome_contigs
    Extend seed directly  = no
    Reference sequence    = 
    Variance detection    = 
    
    Dataset 1:
    -----------------------
    Read Length           = $params.read_length
    Insert size           = $params.insert_size
    Platform              = illumina
    Single/Paired         = PE
    Combined reads        = 
    Forward reads         = ${rawreads[0]}
    Reverse reads         = ${rawreads[1]}
    Store Hash            =
    
    Optional:
    -----------------------
    Insert size auto      = yes
    Use Quality Scores    = no
    Output path           = " > config.txt

    NOVOPlasty.pl -c config.txt
    
    if [[ -f "Circularized_assemblies_1_Mitogenome.fasta" ]]
    then
      cat Circularized_assemblies_1_Mitogenome.fasta > assembled_mitogenome.fasta
      break
    elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
    then
      cat Uncircularized_assemblies_1_Mitogenome.fasta > assembled_mitogenome.fasta
      break
    elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
    then
      bfg Contig01 Contigs_1_Mitogenome.fasta > assembled_mitogenome.fasta
      break
    fi
    """
}

process annotate_mitogenome {
    publishDir "${params.output}/MITOS_annotation", mode: 'copy'
//    cpus ${params.threads}
    label 'process_low'

    input:
    // Fasta file of assembled genome
    path mitogenome
    tuple val(name), path(rawreads)

    output:
    // Mitochondrial genome
    path "*"

    conda '/home/student/anaconda3/envs/mitos'
    script:
    """  
    mkdir -p mitos_output
    runmitos.py -i $mitogenome -o mitos_output -r $params.mitos_reference -R $baseDir -c 05

    mkdir -p individual_genes_nuc
    mkdir -p individual_genes_prot
    if [[ "$params.species_id" ]]
    then
      id=\$( echo "$params.species_id" )
    else
      id=\$( echo "${rawreads[0].simpleName}" )
    fi
    sed "s/^.*\\(; \\)/>\${id}@/g" mitos_output/result.fas | sed 's/(.*//' > individual_genes_nuc/result.fas
    sed "s/^.*\\(; \\)/>\${id}@/g" mitos_output/result.faa | sed 's/(.*//' > individual_genes_prot/result.faa

    cat individual_genes_nuc/result.fas | grep '^>' | sed 's/^.*@//' > individual_genes_nuc.txt
    while read -r line; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_nuc/result.fas > individual_genes_nuc/\$gene.fna; done < individual_genes_nuc.txt
    
    cat individual_genes_prot/result.faa | grep '^>' | sed 's/^.*@//' > individual_genes_prot.txt
    while read -r line; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_prot/result.faa > individual_genes_prot/\$gene.faa; done < individual_genes_prot.txt
    """  
}


workflow {
    extract_mitogenome(ch_contigs, ch_mitogenome)
    if ( extract_mitogenome.out.split_mitgenome ) {
      reassemble_mitogenome(ch_rawReads, extract_mitogenome.out.split_mitgenome ) 
      annotate_mitogenome(reassemble_mitogenome.out.assembled_mitogenome, ch_rawReads) }
    else { annotate_mitogenome(extract_mitogenome.out.mitogenome_contigs, ch_rawReads) }
}

workflow.onComplete {
    // Display complete message
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
    // Display error message
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}