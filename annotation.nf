#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()


ch_mitogenome = Channel
  .fromPath( params.mitogenome, type: 'file' )
//  .first() //transform from queue to value channel for process strand_test

def helpMessage() {
    log.info"""
    Example of command to run:
    nextflow run main.nf --contigs 9GyTe_Gymnonereis_tenera/final_contigs.fasta --reads '9GyTe_Gymnonereis_tenera/raw_reads/NG-29255_9GyTe_lib572147_7903_2_{1,2}.fastq.gz' --mitogenome platynereis_dumerilii_complete_mito.fna --species_id 9_Gymnonereis_tenera --mitos_reference testfolder/  -resume
    """.stripIndent()
}

process annotate_mitogenome {
    publishDir "${params.output}", mode: 'copy'
    label 'process_low'

    input:
    // Fasta file of assembled genome
    path mitogenome

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

    id=\$( echo "${mitogenome.simpleName}" )

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
    

    mkdir -p MITOS_annotation
    mv individual_genes_nuc/ individual_genes_prot/ mitos_output/ *.txt MITOS_annotation
    name=\$( echo "${mitogenome.simpleName}_out" )
    mkdir -p \$name
    mv MITOS_annotation/ \$name
    """
}

workflow {
    annotate_mitogenome(ch_mitogenome)
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
