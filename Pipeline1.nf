#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.out = "${projectDir}/output"
params.store = "${projectDir}/downloads"
params.downloadurl = "https://gitlab.com/dabrowskiw/cq-examples/-/raw/master/data/sequences.sam?inline=false"
params.indir = "${projectDir}/input"



process downloadFile {
  storeDir params.store
  publishDir params.out, mode: "copy", overwrite: true

  output:
    path "sequences.sam"

  script:
  """
  wget ${params.downloadurl} -O sequences.sam
  """
}

process convertToFasta {
  input:
    path "sequences.sam"

  output:
    path "fasta_files"

  script:
  """
  mkdir -p fasta_files

  # Extract each sequence and write to a separate FASTA file
  awk 'BEGIN {i=0} !/^@/ {i++; print ">"\$1"\\n"\$10 > "fasta_files/seq"i".fasta"}' sequences.sam
  """
}



process countStopCodons {
  publishDir params.out, mode: "copy", overwrite: true

  input:
    path "fasta_files"

  output:
    path "stop_codon_counts.txt"

  script:
  """
  > stop_codon_counts.txt
  for file in fasta_files/*.fasta; do
    taa=\$(grep -o TAA "\$file" | wc -l)
    tag=\$(grep -o TAG "\$file" | wc -l)
    tga=\$(grep -o TGA "\$file" | wc -l)
    total=\$((taa + tag + tga))
    echo "\$(basename \$file): \$total stop codons (TAA=\$taa, TAG=\$tag, TGA=\$tga)" >> stop_codon_counts.txt
  done
  """
}


process makeSummaryCSV {
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    path "stop_codon_counts.txt"

    output:
    path "codon_summary.csv"

    script:
    """
    echo "filename,TAA,TAG,TGA" > codon_summary.csv

    while read line; do
        filename=\$(echo "\$line" | cut -d' ' -f1)
        taa=\$(echo "\$line" | grep -o 'TAA=[0-9]*' | cut -d'=' -f2)
        tag=\$(echo "\$line" | grep -o 'TAG=[0-9]*' | cut -d'=' -f2)
        tga=\$(echo "\$line" | grep -o 'TGA=[0-9]*' | cut -d'=' -f2)
        echo "\$filename,\$taa,\$tag,\$tga" >> codon_summary.csv
    done < stop_codon_counts.txt
    """
}


workflow {
  downloadFile()
  convertToFasta(downloadFile.out)
  countStopCodons(convertToFasta.out)
  makeSummaryCSV(countStopCodons.out)
}


