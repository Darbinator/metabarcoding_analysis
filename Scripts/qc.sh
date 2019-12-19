#!/usr/bin/env bash

function run_bbduk() {
	for prefix in $(ls rawdata/ | rev | cut -c 8- | rev | uniq)
	do 
		bbduk.sh in1=rawdata/"$prefix"1.fastq in2=rawdata/"$prefix"2.fastq out=trim/"$prefix"trim_1.fastq out2=trim/"prefix"trim_2.fastq minlength=5
	done
}


function run_bbmerge() {
	for prefix in $(ls rawdata/ | rev | cut -c 8- | rev | uniq)
	do
		bbmerge.sh in1=rawdata/"$prefix"1.fastq in2=rawdata/"$prefix"2.fastq out=output_bbmerge/"$prefix"merged_reads.fastq outu1=output_bbmerge/"$prefix"unmerged_R1.fastq outu2=output_bbmerge/"$prefix"unmerged_R2.fastq ouq=t 
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


function rename_file() {
	for file in $(ls output_bbmerge/*subsample.fasta)
	do
		sample=$(echo "$file" | cut -f2 -d"/" | cut -f1 -d"_")

		mv "$file" subsample_rawdata/"$sample"_contigs.fasta
	done
}

function bbmerge_to_fasta() {
	for f in $(ls output_bbmerge/*fastq)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		mothur "#fastq.info(fastq=$f); make.group(fasta=current, groups=$f, output=$suf.groups)"
	done
}

function make_group() {
	for f in $(ls subsample_rawdata/*fasta)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		mothur "#make.group(fasta=$f, groups=$f, output=$suf.groups)"
	done
}

function to_trim_adapt() {
	for f in $(ls output_bbmerge_trim_adapt/*fastq)
	do
		suf=$(echo "$f" | cut -f2 -d"/")
		python3 trim_adapt.py "$f" > output_bbmerge_trim_adapt/tmp_"$suf"
	done
}



#make.group(fasta=output_bbmerge/11VP10_lib176564_5287_merged_reads.fasta, groups=output_bbmerge/11VP10_lib176564_5287_2_merged_reads.fastq, output=output_bbmerge/11VP10_lib176564_5287_2_merged_reads.fastq.groups)

# egrep '^CGCACCTGGACTGGAC.+TTCTGGCGCTGGACCATGGG$' output_bbmerge/10CF278a_lib176536_5287_2_merged_reads.fasta

run_bbduk


#supprimer score qualités associés

# regex name groups

# sed '{s/output_bbmerge_trim_adapt\///;s/_lib[0-9]\{6\}_5287_2_merged_reads.fastq//}' test/final.good.groups