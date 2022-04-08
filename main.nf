#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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


process extract_mitogenome {
    publishDir "${params.output}", mode: 'copy'
//    cpus ${params.threads}

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
    for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query $mitogenome -db db -outfmt "10 sseqid" -word_size \$i > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        echo "made seqids unique"
        cat possible_mitogenomes.fa | bfg -f unique_seqid.txt > "mitogenome_candidates_wordsize_\$i.fa"
      done

    for file in *candidate*
    do
      grep -v  '^>' \$file | wc -m
    done > nucleotide_count.txt
    closest_match=\$( awk -v c=1 -v t=15843 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count.txt )
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

    for file in *.fa
    do
      contigs=\$( grep '^>' \$file | wc -l )
      nucleotides=\$( grep -v '^>' \$file | wc -m )
      echo "\$contigs \$nucleotides \$file" > \${file}_stats.txt
    done
    stats=\$( cat *_stats.txt )
          echo "Contigs | Nucleotides | Filename
          \$stats" > stats.txt

    """
}

process reassemble_mitogenome {
    publishDir "${params.output}", mode: 'copy'
//    cpus ${params.threads}

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), path(rawreads)
    path(mitogenome_contigs)


    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the host mitogenome
    path('*'), emit: assembled_mitogenome

    script:
    """
    if [[ ${rawreads[0]} == *.gz ]]
    then
      gunzip -f ${rawreads[0]}
      gunzip -f ${rawreads[1]}
    fi
    
    echo "${rawreads[0].simpleName}.fastq ${rawreads[1].simpleName}.fastq $mitogenome_contigs" > test2.txt
    PriceTI -fpp ${rawreads[0].simpleName}.fastq ${rawreads[1].simpleName}.fastq 300 95 -icf $mitogenome_contigs 1 1 3 -nc 10 -dbmax 151 -dbk 91 -mol 30 -tol 20 -mpi 80 -target 90 2 1 1 -lenf 500 2 -a $params.threads -o assembly.fasta

    gzip ${rawreads[0].simpleName}.fastq
    gzip ${rawreads[1].simpleName}.fastq    
    """
}

process annotate_mitogenome {
    publishDir "${params.output}", mode: 'copy'
//    cpus ${params.threads}

    input:
    // Fasta file of assembled genome
    path mitogenome

    output:
    // Mitochondrial genome
    path "*"

    conda '/home/student/anaconda3/envs/mitos'
    script:
    """  
    mkdir -p mkdir mitos_output
    python2 /home/student/anaconda3/envs/mitos/bin/runmitos.py -i $mitogenome -o mitos_output -r 'refseq63m/' -R '/home/student/training_grounds/mitos/testfolder/' -c 05
    """  
}


workflow {
    extract_mitogenome(ch_contigs, ch_mitogenome)
    if ( extract_mitogenome.out.split_mitgenome ) {
        reassemble_mitogenome(ch_rawReads, extract_mitogenome.out.split_mitgenome) }
//       MITOS(reassemble_mitogenome.out.assembled_mitogenome) }
    else {}
//        MITOS(extract_mitogenome.out.mitogenome_contigs) }
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