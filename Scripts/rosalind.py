#!/usr/bin/env python

import sys
import os
import re
from collections import defaultdict


def nameToSeq(file):

	''' Dictionnaire associant le nom d'une séquence et la séquence,
	à partir d'un fichier FASTA multilines '''

	lignes = file.readlines()

	sequences = {}

	for ligne in lignes:

		if ligne.startswith(">"):
			seq = ligne.rstrip()
			sequences[seq] = ''

		else:
			sequences[seq] += ligne.rstrip()

	return sequences



def rnaToProt(seq):
	codonTable = ''' UUU F      CUU L      AUU I      GUU V
		UUC F      CUC L      AUC I      GUC V
		UUA L      CUA L      AUA I      GUA V
		UUG L      CUG L      AUG M      GUG V
		UCU S      CCU P      ACU T      GCU A
		UCC S      CCC P      ACC T      GCC A
		UCA S      CCA P      ACA T      GCA A
		UCG S      CCG P      ACG T      GCG A
		UAU Y      CAU H      AAU N      GAU D
		UAC Y      CAC H      AAC N      GAC D
		UAA Stop   CAA Q      AAA K      GAA E
		UAG Stop   CAG Q      AAG K      GAG E
		UGU C      CGU R      AGU S      GGU G
		UGC C      CGC R      AGC S      GGC G
		UGA Stop   CGA R      AGA R      GGA G
		UGG W      CGG R      AGG R      GGG G '''

	prot=''

	codon_to_aa = {}

	for i in range(0,len(codonTable.split()),2):
		codon_to_aa[codonTable.split()[i]] = codonTable.split()[i+1]

	for j in range(0,len(seq),3):

		if codon_to_aa[seq[j:j+3]] != "Stop":
			prot += codon_to_aa[seq[j:j+3]]

		else:
			break

	return prot


def reverseComplement(seq):
	complement = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

	reverseSeq = seq[::-1]
	revComplementSeq = ''

	for i in range(len(reverseSeq)):
		revComplementSeq += complement[reverseSeq[i]]

	return revComplementSeq



