#!/bin/bash




function run_bbduk() {
	for prefix in $(ls rawdata/ | rev | cut -c 13- | rev | uniq)
	do 
		bbduk.sh in1=rawdata/"$prefix"R1_001.fastq in2=rawdata/"$prefix"R2_001.fastq out=trim/"$prefix"R1_001_trim.fastq out2=trim/"$prefix"R2_001_trim.fastq minlength=50 qtrim=rl trimq=20
	done
}


function run_bbmerge() {
	for prefix in $(ls rawdata/ | rev | cut -c 13- | rev | uniq)
	do
		bbmerge.sh in1=trim/"$prefix"R1_001_trim.fastq in2=trim/"$prefix"R2_001_trim.fastq out=output_bbmerge/"$prefix"merged_reads.fastq outu1=output_bbmerge/"$prefix"unmerged_R1.fastq outu2=output_bbmerge/"$prefix"unmerged_R2.fastq ouq=t 
	done
}

function rename_file_names() {
	for suffix in $(ls output_bbmerge_trim_adapt | cut -d "_" -f2-)
	do
		mv output_bbmerge_trim_adapt/tmp_"$suffix" output_bbmerge_trim_adapt/"$suffix"
	done
}


function sub_sampling() {
	for file in $(ls output_bbmerge/*fasta)
	do
		suf=$(echo "$file" | cut -f2 -d"/")
		mothur "#sub.sample(fasta=$file, size=50000)"
	done
}

function sub_sampling_groups_file() {
	for file in $(ls output_bbmerge_trim_adapt/*fastq.groups)
	do

		head -500 "$file" > "$file"_1000.groups
	done
}

function fastq_to_fasta {
	for file in $(ls output_bbmerge/*fastq)
	do
		NAME=`echo "$file" | cut -d'.' -f1`
		awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' "$file" > "$NAME".fasta
	done
}


function rename_file() {
	for file in $(ls output_bbmerge/*.fasta)
	do
		sample=$(echo "$file" | cut -f2 -d"/" | cut -f1 -d"_")

		mv "$file" output_bbmerge/"$sample"_contigs.fasta
	done
}

rename_file

function bbmerge_to_fasta() {
	for f in $(ls output_bbmerge/*fastq)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		mothur "#fastq.info(fastq=$f); make.group(fasta=current, groups=$f, output=$suf.groups)"
	done
}



function make_group() {
	for f in $(ls output_bbmerge/*fasta)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		echo $suf
		mothur "#make.group(fasta=$f, groups=$f, output=$suf.groups)"
	done
}

function rename_group_files() {
	for f in $(ls output_bbmerge/*groups)
	do
		cut -f1-8 -d"_" "$f" | sed 's/output_bbmerge\///g' > "$f".tmp && mv "$f".tmp "$f"
	done

}

function to_trim_adapt() {
	for f in $(ls output_bbmerge_trim_adapt/*fastq)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		python3 trim_adapt.py "$f" > output_bbmerge_trim_adapt/tmp_"$suf"
	done
}

function run_ITSX() {

	count = 0
	for f in $(ls output_bbmerge/FILTRE*fasta )
	do
		((count++))
		ITSx -i "$f" -o ITSX/FILTRE"$count"
	done
}


run_ITSX



#make.group(fasta=output_bbmerge/11VP10_lib176564_5287_merged_reads.fasta, groups=output_bbmerge/11VP10_lib176564_5287_2_merged_reads.fastq, output=output_bbmerge/11VP10_lib176564_5287_2_merged_reads.fastq.groups)

# egrep '^CGCACCTGGACTGGAC.+TTCTGGCGCTGGACCATGGG$' output_bbmerge/10CF278a_lib176536_5287_2_merged_reads.fasta




#supprimer score qualités associés

# regex name groups

# sed '{s/output_bbmerge_trim_adapt\///;s/_lib[0-9]\{6\}_5287_2_merged_reads.fastq//}' test/final.good.groups





