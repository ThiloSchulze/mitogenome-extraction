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
    path('single_contig_mitogenome.fa'), emit: mitogenome
    path('stats.txt')
    path('split_mitogenome.fa') optional true
    path('NOVOPlasty_out'), type: 'dir' optional true
    path('NOVOPlasty_out_highest_average'), type: 'dir' optional true
    path('unique_mito_seqid.txt') optional true
    path('config.txt') optional true
    path('warning.txt') optional true




    script:
    """
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa
    cat $contigs | bfg "cov_[${params.coverage_cutoff}-9][0-9]{1,}\\.[0-9]+" > cov_${params.coverage_cutoff}0_to_99.fa
    cat $contigs | bfg "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > cov_100_plus.fa
    cat cov_${params.coverage_cutoff}0_to_99.fa > possible_mitogenomes.fa
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

    if [[ \$(wc -l unique_mito_seqid.txt) = "0 unique_mito_seqid.txt" ]] || [[ \$(cat mitogenome_candidates_wordsize_2{3,4,5}.fa | wc -m) == '0' ]]
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
      cat identified_mitogenome.fa > single_contig_mitogenome.fa
    else
      cat identified_mitogenome.fa > split_mitogenome.fa
      seqkit stats *.fa > stats.txt

      echo "Project:
      -----------------------
      Project name          = Mitogenome
      Type                  = mito
      Genome Range          = $params.min_size-$params.max_size
      K-mer                 = $params.kmer_size
      Max memory            = ${task.memory.toGiga()}
      Extended log          = 0
      Save assembled reads  = no
      Seed Input            = split_mitogenome.fa
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
      mkdir -p NOVOPlasty_out
      mv contigs_tmp_Mitogenome.txt log_Mitogenome.txt NOVOPlasty_out
      if [[ -f "Merged_contigs_Mitogenome.txt" ]]
      then
        mv Merged_contigs_Mitogenome.txt NOVOPlasty_out
      fi
      if [[ -f "Circularized_assembly_1_Mitogenome.fasta" ]]
      then
        cat Circularized_assembly_1_Mitogenome.fasta > single_contig_mitogenome.fa
        mv Circularized_assembly_1_Mitogenome.fasta NOVOPlasty_out
      elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
      then
        cat Uncircularized_assemblies_1_Mitogenome.fasta > single_contig_mitogenome.fa
        mv Uncircularized_assemblies_1_Mitogenome.fasta NOVOPlasty_out
      elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
      then
        contig=\$( head -n 1 Contigs_1_Mitogenome.fasta )
        bfg -F \$contig Contigs_1_Mitogenome.fasta > single_contig_mitogenome.fa
        mv Contigs_1_Mitogenome.fasta NOVOPlasty_out
      fi

      if [[ -f "NOVOPlasty_out/Contigs_1_Mitogenome.fasta" ]] && [[ \$( bfg Contig01 NOVOPlasty_out/Contigs_1_Mitogenome.fasta --output-sequences | wc -m ) -lt '12000' ]]
      then
          for blastn_result in *candidates*.fa
          do
              if [[ \$( grep -v '^>' \$blastn_result | wc -m | awk '{print int(\$1+0.5)}' ) -gt '10000' ]]
              then
                  grep '^>' "\$blastn_result" > header_list.txt
                  while read -r header
                      do
                      bfg "\$header" "\$blastn_result" | grep -v '^>' | wc -m
                  done < header_list.txt > "\${blastn_result%.fa}_nuc_per_header.txt"
                  awk 'BEGIN{s=0;}{s+=\$1;}END{print s/NR;}' "\${blastn_result%.fa}_nuc_per_header.txt" > "\${blastn_result}_avg_len.txt"
              fi
          done

          cat *_avg_len.txt | sort -gr | head -1 | cut -d ' ' -f3 > highest_avg.txt


          for avg_len in *_avg_len.txt
          do
            if [[ \$(cat "\$avg_len") = \$(cat highest_avg.txt) ]]
            then
                novoplasty_seed="\${avg_len%_avg_len.txt}"
                break
            fi
          done
          if [[ \$( grep -c '^>' "\$novoplasty_seed" ) == '1' ]]
          then
              cat "\$novoplasty_seed" > single_contig_mitogenome.fa
          else
              cat "\$novoplasty_seed" > split_mitogenome.fa
              NOVOPlasty.pl -c config.txt
              mkdir -p NOVOPlasty_out_highest_average
              mv log_Mitogenome.txt NOVOPlasty_out_highest_average
              if [[ -f "Merged_contigs_Mitogenome.txt" ]]
              then
                mv Merged_contigs_Mitogenome.txt NOVOPlasty_out_highest_average
              fi
              if [[ -f "Circularized_assembly_1_Mitogenome.fasta" ]]
              then
                cat Circularized_assembly_1_Mitogenome.fasta > single_contig_mitogenome.fa
                mv Circularized_assembly_1_Mitogenome.fasta NOVOPlasty_out_highest_average
              elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
              then
                cat Uncircularized_assemblies_1_Mitogenome.fasta > single_contig_mitogenome.fa
                mv Uncircularized_assemblies_1_Mitogenome.fasta NOVOPlasty_out_highest_average
              elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
              then
                bfg -F Contig01 Contigs_1_Mitogenome.fasta > single_contig_mitogenome.fa
                mv Contigs_1_Mitogenome.fasta NOVOPlasty_out_highest_average
              fi
          fi
      fi
    fi

    seqkit stats *.fa > stats.txt
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
    path('single_contig_mitogenome.fa'), emit: strand_tested_mitogenome
    path('original_single_contig_mitogenome.fa') optional true
    path('blast_output.txt')
    path('blast_strands.txt')



    script:
    """
    if [[ \$( cat single_contig_mitogenome.fa | grep -v '^>' | grep -c -i -e [*] ) > '0' ]]
    then
      
      tr -d \\* < single_contig_mitogenome.fa > new_single_contig_mitogenome.fa
      cat new_single_contig_mitogenome.fa > single_contig_mitogenome.fa
      rm new_single_contig_mitogenome.fa
    fi

    makeblastdb -in $mitogenome_reference -dbtype nucl -out reference
    blastn -db reference -query $assembled_mitogenome -word_size 15 -out blast_output.txt
    cat blast_output.txt | grep 'Strand' > blast_strands.txt
    if [[ \$( head -n 1 blast_strands.txt ) == *'Strand=Plus/Minus'* ]]
    then
      cat single_contig_mitogenome.fa > original_single_contig_mitogenome.fa
      revseq -sequence single_contig_mitogenome.fa -reverse -complement -outseq single_contig_mitogenome.fa
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


    script:
    """
    mkdir -p mitos_output
    runmitos.py -i $mitogenome -o mitos_output -r $params.mitos_reference -R $baseDir -c $params.genetic_code > mitos_output.txt

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
    strand_control(ch_mitogenome, extract_mitogenome.out.mitogenome)
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
