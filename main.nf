#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

ch_contigs = Channel
  .fromPath( params.contigs, type: 'file' )
ch_mitogenome = Channel
  .fromPath( params.mitogenome, type: 'file' )
//  .first() //transform from queue to value channel for process strand_test
//ch_barcode = Channel
//  .fromPath( params.barcode, type: 'file' )
if ( params.barcode ) {
  ch_barcode = Channel
      .fromPath (params.barcode, type: 'file')
}

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

    output:
    // Mitogenome (assembled if necessary), NOVOPlasty results, statistics
    path("mito_candidate_*"), emit: mitogenome_candidates
    path('stats.txt') optional true

    conda "${baseDir}/environment1.yml"

    script:
    """
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa
    if [[ "$params.mito_min_size" = 'false' ]] || [[ "$params.mito_min_size" -gt "$params.mito_size" ]]
    then
      threshold_080=\$( echo "$params.mito_size*0.80" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
    else
      threshold_080="$params.mito_min_size"
    fi
    echo " \$threshold_080"
    threshold_100=\$( echo "$params.mito_size" | awk '{printf("%d\\n",\$1 + 0.5)}' )
    threshold_135=\$( echo "$params.mito_size*1.35" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
    threshold_200=\$( echo "$params.mito_size*2" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
    threshold_500=\$( echo "$params.mito_size*3" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
    if [[ "$params.assembler" = 'spades' ]]
    then
      cat $contigs | bfg "cov_[5-9][0-9]{1,}\\.[0-9]+" > cov_50_to_99.fa
      cat $contigs | bfg "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > cov_100_plus.fa
      cat cov_50_to_99.fa cov_100_plus.fa > cov_50_plus.fa
    fi
    cat $contigs > cov_0_plus.fa
    

    makeblastdb -in $contigs -title contig -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query ${params.mitogenome} -db db -outfmt "10 sseqid" -word_size \$i -num_threads ${task.cpus} > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        check_seqid=\$(head -n 1 unique_seqid.txt)
        if [[ \${check_seqid:0-1} = '|' ]]
        then 
          sed -i 's/\\(.*\\)|/\\1 /' unique_seqid.txt
        fi
        echo "made seqids unique"
        if [[ "\$i" = "${params.min_blast_wordsize}" ]]
        then
          create_top_hits_set () {
            sequences=\$((\$seq_limit + 1))
            touch mito_candidate_"\$queue_id"_covcut_0_top_"\$sequences"_match.fa
            while read -r top_hit || [ -n "\$top_hit" ]
            do
              if [[ \$( bfg "\$top_hit" cov_0_plus.fa | tr -d '\n' | wc -m ) -lt "\$threshold_200" ]]
              then
                bfg bfg "\$top_hit" cov_0_plus.fa >> mito_candidate_"\$queue_id"_covcut_0_top_"\$sequences"_match.fa
              fi
              if [[ \$(grep '^>' -c mito_candidate_"\$queue_id"_covcut_0_top_"\$sequences"_match.fa) -gt "\$seq_limit" ]] && [[ \$( cat mito_candidate_"\$queue_id"_covcut_0_top_"\$sequences"_match.fa | tr -d '\n' | wc -m ) -gt "$params.mito_size" ]]
              then 
                break
              fi
            done < unique_seqid.txt
          }
          queue_id='7'; seq_limit='4'; create_top_hits_set
          queue_id='8'; seq_limit='9'; create_top_hits_set
        fi
        if [[ "$params.assembler" = 'spades' ]]
        then
          cat cov_100_plus.fa | bfg -f unique_seqid.txt > "blastn_covcut_100_wordsize_\$i.fa"
          cat cov_50_plus.fa | bfg -f unique_seqid.txt > "blastn_covcut_50_wordsize_\$i.fa"
        fi
        cat cov_0_plus.fa | bfg -f unique_seqid.txt > "blastn_covcut_0_wordsize_\$i.fa"
    done
    for file in blastn_*
    do
      if [[ \$(grep -c  '^>' \$file) -eq '1' ]] && [[ \$(grep -v  '^>' \$file | wc -m) -gt "\$threshold_080" ]]
      then
        cat \$file > mito_candidate_mitogenome.fa
        echo "Found the mitogenome on a single contig."
      fi
    done
    size_match () {
      echo "Starting search for closest mitogenome size match."
      for file in blastn_covcut_\${covcut}_*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count_covcut_\${covcut}.txt
      closest_match=\$( awk -v c=1 -v t=\$threshold 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count_covcut_\${covcut}.txt )
      for blast_result in blastn_covcut_\${covcut}_*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > mito_candidate_\${counter}_covcut_\${covcut}_size_match.fa
          break
        fi
      done
      echo "Finished search for closest mitogenome size match (cov \${covcut})."
    }
    
    if [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      echo "Start size script"
      if [[ "$params.assembler" = 'spades' ]]
      then
        covcut='100'; threshold=\$( echo "\$threshold_100" ); counter='1'; size_match
        covcut='100'; threshold=\$( echo "\$threshold_135" ); counter='extra_1'; size_match
        covcut='50'; threshold=\$( echo "\$threshold_100" ); counter='3'; size_match
        covcut='50'; threshold=\$( echo "\$threshold_135" ); counter='extra_2'; size_match
      fi
      covcut='0'; threshold=\$( echo "\$threshold_100" ); counter='5'; size_match
      covcut='0'; threshold=\$( echo "\$threshold_200" ); counter='extra_3'; size_match
      echo "End size script"

      contig_match () {
      for blastn_result in blastn_covcut_\${covcut}_*
        do
                if [[ \$(cat \$blastn_result | wc -m) != '0' ]]
                then
                grep '^>' "\$blastn_result" > covcut_\${covcut}_header_list.txt
                while read -r header || [ -n "\$header" ]
                    do
                    bfg "\$header" "\$blastn_result" | grep -v '^>' | wc -m
                done < covcut_\${covcut}_header_list.txt > "\${blastn_result%.fa}_covcut_\${covcut}_nuc_per_header.txt"
                awk 'BEGIN{s=0;}{s+=\$1;}END{print s/NR;}' "\${blastn_result%.fa}_covcut_\${covcut}_nuc_per_header.txt" > "\${blastn_result}_covcut_\${covcut}_avg_len.txt"
                fi
        done
        echo "Determined the average nucleotide size per contig for each blast result (cov \${covcut})."
        cat *_covcut_\${covcut}_avg_len.txt | sort -gr | head -1 | cut -d ' ' -f3 > covcut_\${covcut}_highest_avg.txt
        echo "Saved the highest average to the file covcut_\${covcut}_highest_avg.txt (cov 100)."
        for avg_len in *_covcut_\${covcut}_avg_len.txt
        do
          if [[ \$(cat "\$avg_len") = \$(cat covcut_\${covcut}_highest_avg.txt) ]]
          then
              novoplasty_seed="\${avg_len%_covcut_\${covcut}_avg_len.txt}"
              cat \$novoplasty_seed > mito_candidate_\${counter}_covcut_\${covcut}_contig_match.fa
          fi
        done
      }
      echo "Start contig script"
      if [[ "$params.assembler" = 'spades' ]]
      then
        covcut='100'; threshold=\$( echo "\$threshold_100" ); counter='2'; contig_match
        covcut='50'; threshold=\$( echo "\$threshold_100" ); counter='4'; contig_match
      fi
      covcut='0'; threshold=\$( echo "\$threshold_100" ); counter='6'; contig_match
      echo "End contig script"

      seqkit stats *.fa > stats.txt
      echo '0' > candidate_size_list.txt
      for candidate in mito_candidate_*
      do
        nucleotide_count=\$( grep -v '^>' \$candidate | wc -m)
        if grep -Fxq "\$nucleotide_count" candidate_size_list.txt
        then
          rm \$candidate
          echo "Removed candidate \$candidate. A file with \$nucleotide_count nucleotides is already included."
        else
          echo "\$nucleotide_count" >> candidate_size_list.txt
        fi
      done
    fi
    """
}

process reassemble_mitogenome {
    publishDir "${params.output}/mitogenome_reassembly", mode: 'copy'
    label 'process_low'

    input:
    // Assembled contigs fasta file, reference mitogenome, forward and reverse read corresponding to contigs
    path(mitogenomes)
    tuple val(name), path(rawreads)

    output:
    path('single_contig_mitogenome.fa'), emit: mitogenome
    path("NOVOPlasty_run_*"), type: 'dir' optional true
    path('stats.txt') optional true
    path('assembly.log') optional true

    conda "${baseDir}/environment1.yml"

    script:
    """
    touch "assembly.log"
    echo "Initiating reassembly..." >> assembly.log
    if [[ -f mito_candidate_mitogenome.fa ]]
    then
      cat mito_candidate_mitogenome.fa > single_contig_mitogenome.fa
      echo "The mitogenome was identified." >> assembly.log
    elif [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      if [[ "$params.mito_min_size" = 'false' ]] || [[ "$params.mito_min_size" -gt "$params.mito_size" ]]
      then
        threshold_080=\$( echo "$params.mito_size*0.80" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
        threshold_070=\$( echo "$params.mito_size*0.7" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )        
      else
        threshold_080="$params.mito_min_size"
        threshold_070="$params.mito_min_size"
      fi
      threshold_100=\$( echo "$params.mito_size" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
      threshold_200=\$( echo "$params.mito_size*2" | bc | awk '{printf("%d\\n",\$1 + 0.5)}' )
      calc_threshold_500=\$(( "$params.mito_size*5" ))
      threshold_500=\$( echo \$calc_threshold_500 | awk '{printf("%d\\n",\$1 + 0.5)}' )
      
      echo "Target mitochondrion size: \$threshold_100 nucleotides \nMinimum mitochondrion size: \$threshold_080 nucleotides" >> assembly.log

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

      counter='0'
      candidate_list=($mitogenomes)    
      for i in "\${candidate_list[@]}"
      do
        counter=\$((\$counter + 1))
        echo "\n Staring reassembly run \$counter \n =========================" >> assembly.log
        if [[ \$(grep -c '^>' \$i) -eq '1' ]]
        then 
          echo "Seed file \$i contains only a single contig (nucleitode size: \$(grep -v '^>' \$i | tr -d '\n' | wc -m) ). Skipping reassembly." >> assembly.log
          mkdir -p NOVOPlasty_run_\$counter
          cat \$i > largest_single_contig.fa
          mv \$i largest_single_contig.fa NOVOPlasty_run_\$counter
        else
          if [[ \$(grep -v '^>' \$i | wc -m) -eq '0' ]]
          then
            echo "Seed file \$i is empty." >> assembly.log
            continue
          fi
          cat \$i > split_mitogenome.fa
          echo "Starting to run NOVOPlasty" >> assembly.log
          NOVOPlasty.pl -c config.txt
          echo "Finished running NOVOPlasty" >> assembly.log
          mkdir NOVOPlasty_run_\$counter
          mv contigs_tmp_Mitogenome.txt log_Mitogenome.txt NOVOPlasty_run_\$counter
          if [[ -f "Merged_contigs_Mitogenome.txt" ]]
          then
            mv Merged_contigs_Mitogenome.txt NOVOPlasty_run_\$counter
          fi

          separate_contigs () {
          COUNT="0"
          grep "^>" \$input | sort | uniq | while read -r header || [ -n "\$header" ]
          do
            COUNT=\$((\$COUNT + 1))
            echo \$header > contig_name_\${COUNT}.txt
          done
          COUNT="0"
          for header in contig_name_*.txt
          do
              COUNT=\$((\$COUNT + 1))
              search=\$( cat "\$header" )
              PRINT="0"
              while read line || [ -n "\$line" ]
              do
                  if [[ "\$line" = "\$search" ]]
                  then
                      PRINT="1"
                      echo \$line
                      search='re_set_variable'
                      continue
                  fi
                  if [[ \$PRINT = "1" ]] && [[ \${line:0:1} != ">" ]]
                  then
                      echo \$line
                  else
                      PRINT='0'
                  fi
              done < \$input > \${step}_NOVOPlasty_contig_\${COUNT}.fa
          done
          rm contig_name_*.txt
          }

          check_for_mitogenome () {
          echo "Evaluating output..." >> assembly.log
          for contig in post_NOVOPlasty_contig_*.fa
          do
            grep -v "^>" \$contig | wc -m
          done > contig_sizes.txt
          largest_contig=\$( cat contig_sizes.txt | sort -gr | uniq | head -n 1 )
          rm contig_sizes.txt
          for contig in post_NOVOPlasty_contig_*.fa
          do
            if [[ \$(grep -v "^>" \$contig | wc -m) = "\$largest_contig" ]] && [[ \$(grep -v "^>" \$contig | wc -m) -gt "\$threshold_070" ]]
            then
              cat \$contig > largest_single_contig.fa
              echo "The mitochondrion has been successfully reassembled!\n It encompasses \$(grep -v '^>' \$contig | tr -d '\n' | wc -m) nucleotides and was assembled from a selection of \$(grep -c '^>' split_mitogenome.fa) different contigs." >> assembly.log
              break
            fi
          done

          if [[ ! -f "largest_single_contig.fa" ]]
          then
            if [[ "$params.high_coverage" = 'false' ]]
            then
              echo "Mitochondrion reassembly was unsuccessful. Searching seed file for the largest pre-reassembly contig." >> assembly.log
              for contig in pre_NOVOPlasty_contig_*.fa
              do
                if [[ \$(grep -v "^>" \$contig | tr -d '\n' | wc -m) -lt "\$threshold_200" ]]
                then
                  grep -v "^>" \$contig | tr -d '\n' | wc -m
                fi
              done > contig_sizes.txt

              if [[ \$(cat contig_sizes.txt | wc -l) -gt '0' ]]
              then
                closest_match=\$( awk -v c=1 -v t=\$threshold_100 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' contig_sizes.txt )
                for contig in pre_NOVOPlasty_contig_*.fa
                do
                  if [[ \$(grep -v "^>" \$contig | tr -d '\n' | wc -m) = "\$closest_match" ]]
                  then
                    cat \$contig > largest_single_contig.fa
                    echo "The contig closest to the target mitochondrium size of \$threshold_100 is \$contig with a nucleotide size of \$(grep -v "^>" largest_single_contig.fa | tr -d '\n' | wc -m)." >> assembly.log                   
                  fi
                done
              fi
              rm contig_sizes.txt
              
            else
              echo "Mitochondrion reassembly was unsuccessful. Searching seed file for the mitogenome." >> assembly.log
              touch mito_size_range.txt
              for contig in pre_NOVOPlasty_contig_*.fa
              do
                nuc_size=\$(grep -v "^>" \$contig | tr -d '\n' | wc -m)
                if [[ "\$nuc_size" -gt "\$threshold_080" ]] && [[ "\$nuc_size" -lt "\$threshold_200" ]]
                then
                  echo \$nuc_size >> mito_size_range.txt
                fi
              done
              if [[ \$(cat mito_size_range.txt | wc -l) -gt '0' ]]
              then
                cat mito_size_range.txt | sort -gr | uniq > uniq_mito_size_range.txt
                echo "\$(cat uniq_mito_size_range.txt | wc -l) potential mitochondrium candidate(s) found." >> assembly.log
                if [[ ! -f "cox1_archive.pin" ]]
                then 
                  echo "creating cox1 database"
                  makeblastdb -in "$projectDir/refseq63m/featureProt/cox1.fas" -dbtype prot -out cox1_archive
                  echo "creating cob database"
                  makeblastdb -in "$projectDir/refseq63m/featureProt/cob.fas" -dbtype prot -out cob_archive
                fi
  
                while read -r suspect_size || [ -n "\$suspect_size" ]
                do
                  for suspect in pre_NOVOPlasty_contig_*.fa
                  do
                    if [[ \$(grep -v "^>" \$suspect | tr -d '\n' | wc -m) = "\$suspect_size" ]]
                    then
                      echo "start checking contig \$suspect for the mitogenome"
                      blastx -db cox1_archive -query \$suspect -word_size 5 -evalue "1e-100" -query_gencode 5 -outfmt "10 sseqid evalue pident" -num_threads ${task.cpus} -out \${suspect%.fa}_cox1_output.txt
                      if [[ \$(cat "\${suspect}_cox1_output.txt" | wc -l) -eq '0' ]]
                      then
                        blastx -db cob_archive -query \$suspect -word_size 5 -evalue "1e-100" -query_gencode 5 -outfmt "10 sseqid evalue pident" -num_threads ${task.cpus} -out \${suspect%.fa}_cob_output.txt
                      fi
                      if [[ \$(cat "\${suspect}_cox1_output.txt" | wc -l) -gt '0' ]] || [[ \$(cat "\${suspect}_cob_output.txt" | wc -l) -gt '0' ]]
                      then
                        echo "The assembled mitogenome has been found! \nIt is contig \$suspect and encompasses \$(grep -v '^>' \$suspect | tr -d '\n' | wc -m) nucleotides." >> assembly.log
                        cat "\$suspect" > largest_single_contig.fa
                        break 2
                      else
                        rm \${suspect%.fa}_cox1_output.txt \${suspect%.fa}_cob_output.txt
                      fi
                    fi
                  done                 
                done < uniq_mito_size_range.txt
                rm uniq_mito_size_range.txt
              fi
              rm mito_size_range.txt
            fi
          fi
          echo "Finished searching for the mitochondrion." >> assembly.log
          }

          if [[ -f "Circularized_assembly_1_Mitogenome.fasta" ]]
          then
            input="\$i"; step='pre'; separate_contigs
            seqkit stats *.fa* > stats.txt
            if [[ \$(grep -v '^>' Circularized_assembly_1_Mitogenome.fasta | wc -m) -gt '0' ]]
            then
              cat Circularized_assembly_1_Mitogenome.fasta > largest_single_contig.fa
              mv largest_single_contig.fa NOVOPlasty_run_\$counter
              echo " The mitochondrion has been fully circularized!" >> assembly.log
            fi
            mv \$i *_NOVOPlasty_contig_*.fa Circularized_assembly_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter

          elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
          then
            input="\$i"; step='pre'; separate_contigs
            seqkit stats *.fa* > stats.txt
            if [[ \$(grep -v '^>' Uncircularized_assemblies_1_Mitogenome.fasta | wc -m) -gt '0' ]]
            then
              cat Uncircularized_assemblies_1_Mitogenome.fasta > largest_single_contig.fa
              mv largest_single_contig.fa NOVOPlasty_run_\$counter
              echo "The uncircularized mitochondrion has been reassembled!" >> assembly.log
            fi
            mv \$i *_NOVOPlasty_contig_*.fa Uncircularized_assemblies_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter

          elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
          then
              echo "Mitogenome was not circularized."
              input="\$i"; step='pre'; separate_contigs
              input='Contigs_1_Mitogenome.fasta'; step='post'; separate_contigs
              check_for_mitogenome
              seqkit stats *.fa* > stats.txt
              mv \$i *_NOVOPlasty_contig_*.fa Contigs_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter
              if [[ -f largest_single_contig.fa ]]
              then
                mv largest_single_contig.fa NOVOPlasty_run_\$counter
              fi

          fi
          if [[ -f NOVOPlasty_run_\${counter}/largest_single_contig.fa ]] && [[ \$(grep -v '^>' NOVOPlasty_run_\${counter}/largest_single_contig.fa | wc -m) -gt "\$threshold_070" ]]
          then
            cat NOVOPlasty_run_\${counter}/largest_single_contig.fa > single_contig_mitogenome.fa
            cp single_contig_mitogenome.fa NOVOPlasty_run_\${counter}
            echo "Starting annotation." >> assembly.log
            break
          fi
        fi
      done

      check_contigs=( NOVOPlasty_run_*/largest_single_contig.fa )
      if [[ ! -f single_contig_mitogenome.fa ]] && [[ -f "\$check_contigs" ]]
      then
        touch largest_contigs_list.txt
        for dir in NOVOPlasty_run_*/
        do
          if [[ -f "\${dir}"largest_single_contig.fa ]]
          then
            grep -v "^>" \${dir}largest_single_contig.fa | wc -m
          fi
        done > contig_sizes.txt
        largest_contig=\$( cat contig_sizes.txt | sort -gr | uniq | head -n 1 )
        rm contig_sizes.txt
        for dir in NOVOPlasty_run_*/
        do
          if [[ \$(grep -v "^>" "\${dir}"largest_single_contig.fa | wc -m) = "\$largest_contig" ]]
          then
            cat "\${dir}"largest_single_contig.fa > single_contig_mitogenome.fa
          fi
        done
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
    path('mitogenome.fa'), emit: strand_tested_mitogenome
    path('original_single_contig_mitogenome.fa') optional true
    path('blast_output.txt')
    path('blast_strands.txt')


    conda "${baseDir}/environment1.yml"

    script:
    """
    if [[ \$( cat $assembled_mitogenome | grep -v '^>' | grep -c -i -e [*] ) > '0' ]]
    then
      
      tr -d \\* < $assembled_mitogenome > new_single_contig_mitogenome.fa
      cat new_single_contig_mitogenome.fa > single_contig_mitogenome.fa
      rm new_single_contig_mitogenome.fa
    fi

    makeblastdb -in $mitogenome_reference -dbtype nucl -out reference
    blastn -db reference -query $assembled_mitogenome -word_size 15 -out blast_output.txt
    cat blast_output.txt | grep 'Strand' > blast_strands.txt
    if [[ \$( head -n 1 blast_strands.txt ) == *'Strand=Plus/Minus'* ]]
    then
      cat "${assembled_mitogenome}" > original_single_contig_mitogenome.fa
      revseq -sequence single_contig_mitogenome.fa -reverse -complement -outseq mitogenome.fa
    else
      cat "${assembled_mitogenome}" > mitogenome.fa
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
    path('individual_genes_nuc/cox1.fna'), emit: cox1 optional true

    conda "${baseDir}/environment2.yml"

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
    while read -r line || [ -n "\$line" ]; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_nuc/result.fas > individual_genes_nuc/\$gene.fna; done < individual_genes_nuc.txt

    cat individual_genes_prot/result.faa | grep '^>' | sed 's/^.*@//' > individual_genes_prot.txt
    while read -r line || [ -n "\$line" ]; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_prot/result.faa > individual_genes_prot/\$gene.faa; done < individual_genes_prot.txt

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

    if grep -q '\\-' "mitos_output/result.geneorder"
    then
      cat mitos_output/result.geneorder > mitos_output/original_result.geneorder
      sed -i -e 's/-/_/g' mitos_output/result.geneorder
    fi
    if grep -q 'cox1' "mitos_output/result.geneorder"
    then
      grep -v '^>' mitos_output/result.geneorder > mitos_output/current_order.txt
      while read -r line || [ -n "\$line" ]; do
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

process get_barcode {
    publishDir "${params.output}/MITOS_annotation", mode: 'copy'
    label 'process_low'

    input:
    // Fasta file of assembled genome
    path barcode_ref
    path cox1

    output:
    // Mitochondrial genome
    path('barcode'), type: 'dir' optional true

    conda "${baseDir}/environment1.yml"

    script:
    """
    if [[ -f $cox1 ]]
    then
      makeblastdb -in $barcode_ref -dbtype nucl -out platynereis_ref
      blastn -query $cox1 -db platynereis_ref -word_size 13 -outfmt "10 qseq" > blast_out.fa
      tr -d \\- < blast_out.fa > cox1_subunit_seq.fna
      head -n 1 $cox1 > header.txt
      sed 's/@cox1/@cox1_subunit/' header.txt > header_edit.txt
      cat header_edit.txt cox1_subunit_seq.fna > cox1_subunit.fna
      mkdir -p barcode
      mv cox1_subunit.fna barcode/
    fi
    """
}

workflow {
    extract_mitogenome(ch_contigs, ch_mitogenome)
    reassemble_mitogenome(extract_mitogenome.out.mitogenome_candidates, ch_rawReads)
    strand_control(ch_mitogenome, reassemble_mitogenome.out.mitogenome)
    annotate_mitogenome(strand_control.out.strand_tested_mitogenome, ch_rawReads)
    if (!params.barcode){
    }
    else {
      get_barcode(ch_barcode, annotate_mitogenome.out.cox1)
    }
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
