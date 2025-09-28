//parameters
params.out="${projectDir}/output"
params.inDir= null
params.inURL= null
params.temp="${projectDir}/downloads"

//processes
process downloadFile{
    storeDir params.temp
    output:
        path "sequences.sam"
    """
    wget ${params.inURL} -O sequences.sam
    """
}

process splitSeqs {
    input:
        path inputfile
    output:
        path "x*.txt"
    """
    tail -n +3 ${inputfile}| split -l 1 --additional-suffix ".txt" -d
    """
    //alternative
    //grep -ve "^@" ${inputfile} | split -l 1 --additional-suffix ".txt" -d
}  

process makeFasta {
	publishDir params.out, mode: "copy", overwrite:true
    input:
        path inputfiles
    output:
        path "*.fasta"
    """
    echo -n ">" > ${inputfiles.getSimpleName()}.fasta
    cat ${inputfiles} | cut -f 1 >> ${inputfiles.getSimpleName()}.fasta
    cat ${inputfiles} | cut -f 10 >> ${inputfiles.getSimpleName()}.fasta 
    """
    
}
process getSeqName {
    //publishDir params.out, mode: "copy", overwrite:true
    input:
        path inputfiles
    output:
        path "x*_name.txt"
    """
    cat ${inputfiles} | cut -f 1 > ${inputfiles.getSimpleName()}_name.txt
    """
}


process countStart{
	//publishDir params.out, mode: "copy", overwrite:true
	input:
		path fastafile
	output:
		path "*_startcount.txt"
		// .getSimpleName() removes the extension of a file
	"""
    grep -o "ATG" ${fastafile} | wc -l > ${fastafile.getSimpleName()}_startcount.txt
	"""
}
process countStop{
	//publishDir params.out, mode: "copy", overwrite:true
	input:
		path fastafile
	output:
		path "*_stopcount.txt"
		// .getSimpleName() removes the extension of a file
	"""
    grep -o -E "TAA|TAG|TGA" ${fastafile} | wc -l > ${fastafile.getSimpleName()}_stopcount.txt
	"""
}

process summaryRep{
    publishDir params.out, mode: "copy", overwrite:true
	input:
		path inputfiles
	output:
		path "summary.csv"
	"""
	for f in \$(ls *.txt); do 
        cat \$f |tr -d "\n";
        echo -n ", ";
	done > summary.csv
	"""
    //if you want to print the summary with linebreaks use this command : cat summary.csv | sed -e 's/, S/\nS/g'
}

//workflow
workflow {
   if(params.inURL != null && params.inDir == null){
        c_download=downloadFile()
    } 
    else if(params.inDir != null && params.inURL == null){
        c_download=channel.fromPath("{params.inDir}/*.sam")
    } else {
        print "Error: Please provide either --inDir or --inURL"
        System.exit(0)
    }
    c_flat = splitSeqs(c_download) | flatten
    c_fasta= makeFasta(c_flat) | flatten

    /*
    c_seq= getSeqName(c_flat) | collect
    c_start= countStart(c_fasta) | collect
    c_stop = countStop(c_fasta) | collect
    c_report=c_seq.combine(c_start.combine(c_stop))| summaryRep
    */

    c_seq= getSeqName(c_flat) 
    c_start= countStart(c_fasta)
    c_stop = countStop(c_fasta) 
    c_report=c_seq.concat(c_start, c_stop)|collect|summaryRep
    
    }