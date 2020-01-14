#!/usr/bin/env python

import sys
import re
import os
import pandas



def read_otu_table(file):

	try:
		otu_table = pandas.read_table(file,sep="\t",header=0)
		return otu_table
	except:
		exit(0)


def nameToSeq(file):

	''' Dictionnaire associant le nom d'une séquence et la séquence,
	à partir d'un fichier FASTA multilines '''

	lignes = file.readlines()

	sequences = {}

	for ligne in lignes:

		if ligne.startswith(">"):
			seq = ligne.rstrip().lstrip('>')
			sequences[seq] = ''

		else:
			sequences[seq] += ligne.rstrip()

	return sequences


df = read_otu_table(sys.argv[1])

total = 0
for i in df.columns[3:]:
	total+=df[i].sum()



conserve = 0
for i in df.columns[3:]:

	if (df[i].sum()/total) * 100 < 0.1:
		del df[i]
	else:

		conserve+= df[i].sum()

print(conserve/total*100)

otus_conserved = list(df.columns[3:])

df.to_csv(path_or_buf="stats/shared_new_new_01.tsv",sep="\t", header=True,index=False)

def read_list_file(file):
	try:
		liste_file = pandas.read_table(file,sep="\t",header=0)
		return liste_file
	except:
		exit(0)


liste = read_list_file(sys.argv[2])

otu_to_ref = {}
for i in liste.columns:
	if i in otus_conserved:
		print(i)
		otu_to_ref[i] = str(liste[i].values[0]).split(',')[0]

print(otu_to_ref)

file_3 = open(sys.argv[3],'r')

fasta = nameToSeq(file_3)

to_write = open('80.fasta','w')


for otu,seq in otu_to_ref.items():
	if seq in fasta:
		to_write.write('>'+otu+'\n')
		to_write.write(fasta[seq]+'\n')



