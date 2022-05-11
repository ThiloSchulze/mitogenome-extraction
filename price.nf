#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

ch_contigs = Channel
  .fromPath( params.contigs, type: 'file' )
ch_mitogenome = Channel
  .fromPath( params.mitogenome, type: 'file' )
//  .first() //transform from queue to value channel for process strand_test
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
    label 'process_low'

    input:
    // Assembled contigs fasta file, reference mitogenome, forward and reverse read corresponding to contigs
    path(contigs)
    path(mitogenome)
    tuple val(name), path(rawreads)

    output:
    // Mitogenome (assembled if necessary), NOVOPlasty results, statistics
    path('assembled_mitogenome.fasta') optional true
    path('stats.txt')
    path('split_mitogenome.fa'), emit: split_mitogenome
    path('NOVOPlasty_out'), type: 'dir' optional true
    path('NOVOPlasty_out_highest_average'), type: 'dir' optional true
    path('unique_mito_seqid.txt') optional true
    path('config.txt') optional true
    path('warning.txt') optional true




    conda "${baseDir}/environment1.yml"

    script:
    """
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa
    cat $contigs > possible_mitogenomes.fa
    makeblastdb -in $contigs -title contig -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
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
        cat $contigs | bfg "cov_[1-9][0-9]{1,}\\.[0-9]+" > cov_10_plus.fa
      for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
        do
          echo "starting iteration with word size \$i"
          cat unique_seqid.txt > prev_seqid.txt
          blastn -query $mitogenome -db db -outfmt "10 sseqid" -word_size \$i -num_threads ${task.cpus} > seqid.txt
          echo "blastn complete"
          cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
          echo "made seqids unique"
          cat cov_10_plus.fa | bfg -f unique_seqid.txt > "mitogenome_candidates_wordsize_\$i.fa"
        done

      for file in *candidate*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count.txt
      closest_match=\$( awk -v c=1 -v t=$params.nucleotide_size 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count.txt )
      for blast_result in *candidate*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > identified_mitogenome.fa
          break
        fi
      done
    fi

    grep '^>' identified_mitogenome.fa | cat -n | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_mito_seqid.txt

    if [[ \$(wc -l unique_mito_seqid.txt) = "0 unique_mito_seqid.txt" ]]
    then
      echo "The mitogenome was not identified! It may have a coverage below 10." > warning.txt
    elif [[ \$(wc -l unique_mito_seqid.txt) = "1 unique_mito_seqid.txt" ]]
    then
      cat identified_mitogenome.fa > assembled_mitogenome.fasta
    else
      cat identified_mitogenome.fa > split_mitogenome.fa
    fi

    seqkit stats *.fa > stats.txt
    """
}


process reassemble_mitogenome {
    publishDir "${params.output}/Price", mode: 'copy'

    label 'process_high'

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
    if [[ ${rawreads[0]} == *.gz ]]
    then
      gunzip -f ${rawreads[0]}
      gunzip -f ${rawreads[1]}
    fi
    
    echo "${rawreads[0].simpleName}.fastq ${rawreads[1].simpleName}.fastq $mitogenome_contigs" > test2.txt
    PriceTI -fpp ${rawreads[0].simpleName}.fastq ${rawreads[1].simpleName}.fastq 300 95 -icf $mitogenome_contigs 1 1 3 -nc 5 -dbmax 151 -dbk 91 -mol 30 -tol 20 -mpi 80 -target 90 2 1 1 -lenf 300 0 -lenf 500 2 -a ${task.cpus} -o assembly.fasta

    if [[ -f assembly.cycle5.fasta ]]
    then 
      if [[ \$(grep '^>' -c assembly.cycle5.fasta) = '1' ]]
      then
        cat assembly.cycle5.fasta > assembled_mitogenome.fasta
      else
        bfg -F contig_1 assembly.cycle5.fasta > assembled_mitogenome.fasta
      fi
    fi

    seqkit stats *.fasta > stats_price.fa
    """
}

process strand_control {
    publishDir "${params.output}/strand_control", mode: 'copy'
    label 'process_low'

    input:
    // Fasta file of assembled genome
    path mitogenome_reference
    path assembled_mitogenome

    output:
    // Mitochondrial genome
    path('assembled_mitogenome.fasta'), emit: strand_tested_mitogenome
    path('original_assembled_mitogenome.fasta') optional true
    path('blast_output.txt')
    path('blast_strands.txt')


    conda "${baseDir}/environment1.yml"

    script:
    """
    if [[ \$( cat assembled_mitogenome.fasta | grep -v '^>' | grep -c -i -e [*] ) > '0' ]]
    then
      
      tr -d \\* < assembled_mitogenome.fasta > new_assembled_mitogenome.fasta
      cat new_assembled_mitogenome.fasta > assembled_mitogenome.fasta
      rm new_assembled_mitogenome.fasta
    fi

    makeblastdb -in $mitogenome_reference -dbtype nucl -out reference
    blastn -db reference -query $assembled_mitogenome -word_size 15 -out blast_output.txt
    cat blast_output.txt | grep 'Strand' > blast_strands.txt
    if [[ \$( head -n 1 blast_strands.txt ) == *'Strand=Plus/Minus'* ]]
    then
      cat assembled_mitogenome.fasta > original_assembled_mitogenome.fasta
      revseq -sequence assembled_mitogenome.fasta -reverse -complement -outseq assembled_mitogenome.fasta
    fi
    """
}

process annotate_mitogenome {
    publishDir "${params.output}/MITOS_annotation", mode: 'copy'
    label 'process_low'

    input:
    // Fasta file of assembled genome
    path mitogenome
    tuple val(name), path(rawreads)

    output:
    // Mitochondrial genome
    path "*"

    conda "${baseDir}/environment2.yml"

    script:
    """
    mkdir -p mitos_output
    runmitos.py -i $mitogenome -o mitos_output -r $params.mitos_reference -R $baseDir -c 05 > mitos_output.txt

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

    if [[ -f individual_genes_nuc/nad4.fna ]]
    then
    bfg -v -F nad4l individual_genes_nuc/nad4.fna > nad4.fna
    mv nad4.fna individual_genes_nuc/nad4.fna
    fi

    if [[ -f individual_genes_prot/nad4.faa ]]
    then
    bfg -v -F nad4l individual_genes_prot/nad4.faa > nad4.faa
    mv nad4.faa individual_genes_prot/nad4.faa
    fi

    if grep -q '\\-cox1' "mitos_output/result.geneorder"
    then
      cat mitos_output/result.geneorder > mitos_output/original_result.geneorder
      sed -i -e 's/-//g' mitos_output/result.geneorder
    fi
    if grep -q 'cox1' "mitos_output/result.geneorder"
    then
      grep -v '^>' mitos_output/result.geneorder > mitos_output/current_order.txt
      while read -r line; do
          if [[ \$( cat mitos_output/current_order.txt | awk '{print \$1;}' ) == *"cox1"* ]]
          then
              cat mitos_output/current_order.txt > mitos_output/adjusted_result.geneorder
          else
              last_gene=\$( awk '{ print \$NF }' mitos_output/current_order.txt )
              new_order=\$( sed "s/\\<\$last_gene\\>//" mitos_output/current_order.txt )
              echo "\$last_gene \$new_order" > mitos_output/current_order.txt
          fi
      done < mitos_output/current_order.txt
    fi
    """
}

workflow {
    extract_mitogenome(ch_contigs, ch_mitogenome, ch_rawReads)
    reassemble_mitogenome(ch_rawReads, extract_mitogenome.out.split_mitogenome)
    strand_control(ch_mitogenome, reassemble_mitogenome.out.assembled_mitogenome)
    annotate_mitogenome(strand_control.out.strand_tested_mitogenome, ch_rawReads)
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
