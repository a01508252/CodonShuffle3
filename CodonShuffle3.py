import os
import numpy as np
import scipy
import matplotlib
import pandas as pd
import statsmodels
import patsy
import sys
import argparse
import matplotlib.pyplot as plt
import re
from random import shuffle, random, randint, choice, seed
from collections import Counter
from subprocess import call
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from plotnine import ggplot, aes, geom_density, labs, geom_vline, ggtitle, geom_point, ggsave

# Codon translation table
tt = {
    "TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser","TCA":"S|Ser","TCG":"S|Ser", 
    "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp","TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp",
    "CTT":"L|Leu","CTC":"L|Leu","CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
    "CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg","CGA":"R|Arg","CGG":"R|Arg",
    "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met","ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr",
    "AAT":"N|Asn","AAC":"N|Asn","AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
    "GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala","GCA":"A|Ala","GCG":"A|Ala",
    "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu","GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"
}

nts = ['A','C','G','T']

def gc3(seq):
    """
    Analyzes and shuffles a DNA sequence based on GC and AT content at the third codon position.
    This function calculates the GC and AT content at the third position of each codon in the given DNA sequence.
    It then shuffles the sequence while preserving the overall GC and AT content at the third codon position.
    Args:
        seq (str): A DNA sequence composed of the characters 'A', 'T', 'G', and 'C'.
    Returns:
        str: A shuffled DNA sequence with preserved GC and AT content at the third codon position.
    """
    gc=0
    at=0
    for num in range(2,len(seq),3):
        if seq[num]=='A'or seq[num]=='T':
            at+=1
        elif seq[num]=='G'or seq[num]=='C':
            gc+=1
    at=at/(len(seq)/3.)
    gc=gc/(len(seq)/3.)
    
    seq1=[]
    for num in range(2,len(seq),3):
        seq1+=seq[num-2:num],
        if (seq[num]=='T'or seq[num]=='C')and(seq[num-2:num]=='TT'or seq[num-2:num]=='TA'or seq[num-2:num]=='TG'or seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GA'):
            seq1+='_Y_',
        elif (seq[num]=='A'or seq[num]=='G')and(seq[num-2:num]=='TT'or seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GA'):
            seq1+='_R_',
        elif seq[num-2:num+1]=='ATT'or seq[num-2:num+1]=='ATC'or seq[num-2:num+1]=='ATA':
            seq1+='_H_',
        elif (seq[num]=='A'or seq[num]=='G'or seq[num]=='T'or seq[num]=='C')and(seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            seq1+='_N_',
        else: seq1+=seq[num],
    seq2=''
    for i in seq1:
        if i == '_Y_':
            x=random()
            if x<=gc:
                seq2+='C'
            elif gc<x<=gc+at:
                seq2+='T'
            else: seq2+=choice('TC')
        elif i == '_R_':
            x=random()
            if x<=gc:
                seq2+='G'
            elif gc<x<=gc+at:
                seq2+='A'
            else: seq2+=choice('AG')
        elif i == '_H_':
            x=random()
            if x<=gc:
                seq2+='C'
            elif gc<x<=gc+at:
                seq2+=choice('AT')
            else: seq2+=choice('ATC')
        elif i == '_N_':
            x=random()
            if x<=gc:
                seq2+=choice('GC')
            elif gc<x<=gc+at:
                seq2+=choice('AT')
            else: seq2+=choice('AGTC')
        else: seq2+=i
    seq=seq2
    return seq

def third_simple(seq):
    """
    Shuffles specific nucleotides in a DNA sequence based on certain conditions.

    This function processes the input DNA sequence in four stages, each targeting specific nucleotides
    (Y, R, H, N) and shuffling them according to predefined conditions.

    Args:
        seq (str): The input DNA sequence.

    Returns:
        str: The DNA sequence after shuffling specific nucleotides.

    Notes:
        - Y nucleotides (T or C) are shuffled if they appear in specific codons.
        - R nucleotides (A or G) are shuffled if they appear in specific codons.
        - H nucleotides (A, C, or T) are shuffled if they appear in specific codons.
        - N nucleotides (A, C, T, or G) are shuffled if they appear in specific codons.
    """
    Y=[]
    seq1=[]
    for num in range(2,len(seq),3):
        if (seq[num]=='T' or seq[num]=='C')and(seq[num-2:num]=='TT'or seq[num-2:num]=='TC'or seq[num-2:num]=='TA'or seq[num-2:num]=='TG'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CA'or seq[num-2:num]=='CG'or seq[num-2:num]=='AT'or seq[num-2:num]=='AC'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GU' or seq[num-2:num]=='GC'or seq[num-2:num]=='GA'or seq[num-2:num]=='GG'):
            Y+=seq[num],
            seq1+=seq[num-2:num],'_Y_',
        else:seq1+=seq[num-2:num+1],
    shuffle(Y)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='_Y_':seq2+=Y.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    R=[]
    seq1=[]
    for num in range(2,len(seq),3):
        if (seq[num]=='A' or seq[num]=='G')and(seq[num-2:num]=='TT'or seq[num-2:num]=='CA'or seq[num-2:num]=='AA'or seq[num-2:num]=='AG'or seq[num-2:num]=='GA'):
            R+=seq[num],
            seq1+=seq[num-2:num],'_R_',
        else:seq1+=seq[num-2:num+1],
    shuffle(R)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='_R_':seq2+=R.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    H=[]
    seq1=[]
    for num in range(2,len(seq),3):
        if (seq[num]=='A'or seq[num]=='C'or seq[num]=='T')and(seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AT'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT' or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            H+=seq[num],
            seq1+=seq[num-2:num],'_H_',
        else:seq1+=seq[num-2:num+1],
    shuffle(H)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='_H_':seq2+=H.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    N=[]
    seq1=[]
    for num in range(2,len(seq),3):
        if (seq[num] in ['A','C','T','G']) and (seq[num-2:num]=='TC'or seq[num-2:num]=='CT'or seq[num-2:num]=='CC'or seq[num-2:num]=='CG'or seq[num-2:num]=='AC'or seq[num-2:num]=='GT'or seq[num-2:num]=='GC'or seq[num-2:num]=='GG'):
            N+=seq[num],
            seq1+=seq[num-2:num],'_N_',
        else:seq1+=seq[num-2:num+1],
    shuffle(N)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='_N_':seq2+=N.pop(0)
        else:seq2+=seq1[i]
    seq=seq2        
    return seq

def dn23(seq):
    """
    Analyzes and shuffles codon sequences in a given DNA sequence.

    This function takes a DNA sequence as input, counts the occurrences of 
    each possible dinucleotide pair, and then generates a new sequence by 
    shuffling the codons based on the observed frequencies.

    Args:
        seq (str): A string representing the DNA sequence to be analyzed and shuffled.

    Returns:
        str: A new DNA sequence with shuffled codons based on the observed dinucleotide frequencies.

    Notes:
        - The input sequence should be a multiple of 3 in length.
        - The function avoids creating stop codons (TAA, TAG, TGA) in the new sequence.
    """
    aa=ag=ac=at=ga=gg=gc=gt=ca=cg=cc=ct=ta=tg=tc=tt=0
    for num in range(2,len(seq),3):
        if seq[num-1]=='A':
            if seq[num]=='A':
                aa+=1
            elif seq[num]=='G':
                ag+=1
            elif seq[num]=='C':
                ac+=1
            elif seq[num]=='T':
                at+=1
        elif seq[num-1]=='G':
            if seq[num]=='A':
                ga+=1
            elif seq[num]=='G':
                gg+=1
            elif seq[num]=='C':
                gc+=1
            elif seq[num]=='T':
                gt+=1
        elif seq[num-1]=='C':
            if seq[num]=='A':
                ca+=1
            elif seq[num]=='G':
                cg+=1
            elif seq[num]=='C':
                cc+=1
            elif seq[num]=='T':
                ct+=1
        elif seq[num-1]=='T':
            if seq[num]=='A':
                ta+=1
            elif seq[num]=='G':
                tg+=1
            elif seq[num]=='C':
                tc+=1
            elif seq[num]=='T':
                tt+=1
    aa,ag,ac,at,ga,gg,gc,gt,ca,cg,cc,ct,ta,tg,tc,tt=aa/(len(seq)/3.),ag/(len(seq)/3.),ac/(len(seq)/3.),at/(len(seq)/3.),ga/(len(seq)/3.),gg/(len(seq)/3.),gc/(len(seq)/3.),gt/(len(seq)/3.),ca/(len(seq)/3.),cg/(len(seq)/3.),cc/(len(seq)/3.),ct/(len(seq)/3.),ta/(len(seq)/3.),tg/(len(seq)/3.),tc/(len(seq)/3.),tt/(len(seq)/3.)
    seq2=''
    for num in range(2,len(seq),3):
        seq2+=seq[num-2:num]
        if seq[num-1]=='A'and seq[num-2:num+1]!='TAA'and seq[num-2:num+1]!='TAG':
            if seq[num]=='T'or seq[num]=='C':
                space=at+ac
                AT,AC=at/space,ac/space
                x=random()
                if x<=AT:
                    seq2+='T'
                elif AT<x<=AT+AC:
                    seq2+='C'
            elif seq[num]=='A'or seq[num]=='G':
                space=aa+ag
                AA,AG=aa/space,ag/space
                x=random()
                if x<=AA:
                    seq2+='A'
                elif AA<x<=AA+AG:
                    seq2+='G'
            else:seq2+=seq[num]
        elif seq[num-1]=='G'and seq[num-2:num+1]!='TGA'and seq[num-2:num+1]!='TGG':
            if (seq[num-2]=='T'or seq[num-2]=='A')and(seq[num]=='C'or seq[num]=='T'):
                space = gt+gc
                GT,GC=gt/space,gc/space
                x=random()
                if x<=GT:
                    seq2+='T'
                elif GT<x<=GT+GC:
                    seq2+='C'
            elif seq[num-2:num+1]=='AGA'or seq[num-2:num+1]=='AGG':
                space=ga+gg
                GA,GG=ga/space,gg/space
                x=random()
                if x<=GA:
                    seq2+='A'
                elif GA<x<=GA+GG:
                    seq2+='G'
            elif seq[num-2]=='C'or seq[num-2]=='G':
                space=ga+gg+gc+gt
                GA,GG,GC,GT=ga/space,gg/space,gc/space,gt/space
                x=random()
                if x<=GA:seq2+='A'
                elif GA<x<=GA+GG:seq2+='G'
                elif GA+GG<x<=GA+GG+GC:seq2+='C'
                elif GA+GG+GC<x<=GA+GG+GC+GT:seq2+='T'
            else:seq2+=seq[num]
        elif seq[num-1]=='C':
            space=ca+cg+cc+ct
            CA,CG,CC,CT=ca/space,cg/space,cc/space,ct/space
            x=random()
            if x<=CA:seq2+='A'
            elif CA<x<=CA+CG:seq2+='G'
            elif CA+CG<x<=CA+CG+CC:seq2+='C'
            elif CA+CG+CC<x<=CA+CG+CC+CT:seq2+='T'
        elif seq[num-1]=='T':
            if seq[num-2:num+1]=='TTT'or seq[num-2:num+1]=='TTC':
                space = tt+tc
                TT,TC=tt/space,tc/space
                x=random()
                if x<=TT:
                    seq2+='T'
                elif TT<x<=TT+TC:
                    seq2+='C'
            elif seq[num-2:num+1]=='TTA'or seq[num-2:num+1]=='TTG':
                space = ta+tg
                TA,TG=ta/space,tg/space
                x=random()
                if x<=TA:
                    seq2+='A'
                elif TA<x<=TA+TG:
                    seq2+='G'
            elif seq[num-2:num+1]=='ATT'or seq[num-2:num+1]=='ATC'or seq[num-2:num+1]=='ATA':
                space=tt+tc+ta
                TT,TC,TA=tt/space,tc/space,ta/space
                x=random()
                if x<=TA:seq2+='A'
                elif TA<x<=TA+TC:seq2+='C'
                elif TA+TC<x<=TA+TC+TT:seq2+='T'
            elif seq[num-2]=='C'or seq[num-2]=='G':
                space=ta+tg+tc+tt
                TA,TG,TC,TT=ta/space,tg/space,tc/space,tt/space
                x=random()
                if x<=TA:seq2+='A'
                elif TA<x<=TA+TG:seq2+='G'
                elif TA+TG<x<=TA+TG+TC:seq2+='C'
                elif TA+TG+TC<x<=TA+TG+TC+TT:seq2+='T'
            else:seq2+=seq[num]
        else:seq2+=seq[num]
    seq=seq2
    return seq

def third(seq: str) -> str:
    """Shuffle nucleotides in the given sequence to preserve certain dinucleotide frequencies.

    This function creates a synonymously mutated sequence maintaining dinucleotide content in 
    codon positions 2–3 and 3–1 identical to that of the input sequence. It applies multiple 
    rounds of codon position shuffling within amino acid groups (like four-codon sets, six-codon 
    sets, etc.) to preserve the underlying dinucleotide structure. 

    Args:
        seq (str): The input nucleotide sequence (typically coding DNA).

    Returns:
        str: The synonymously mutated sequence with preserved dinucleotide frequencies at specific 
        codon positions.
    """
    seq1 = []
    TNT,TNC,TNA,TNG,GNT,GNA,GNC,GNG,CNG,CNA,CNT,CNC = [],[],[],[],[],[],[],[],[],[],[],[]
    for num in range(2,len(seq)-3,3):
        seq1 += seq[num-2:num]
        if seq[num] in ['T','C','A','G']:
            if seq[num-2:num] in ['CT','GT']:  # LEU4 or VAL
                if seq[num+1]=='T':
                    seq1 += 'TNT',
                    TNT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1 += 'TNC',
                    TNC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1 += 'TNA',
                    TNA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1 += 'TNG',
                    TNG.append(seq[num])
                else:
                    seq1 += seq[num]
            elif seq[num-1]=='C':  # SER4 or PRO or THR or ALA
                if seq[num+1]=='T':
                    seq1 += 'CNT',
                    CNT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1 += 'CNC',
                    CNC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1 += 'CNA',
                    CNA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1 += 'CNG',
                    CNG.append(seq[num])
                else:
                    seq1 += seq[num]
            elif seq[num-2:num] in ['CG','GG']:  # ARG4 or GLY
                if seq[num+1]=='T':
                    seq1 += 'GNT',
                    GNT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1 += 'GNC',
                    GNC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1 += 'GNA',
                    GNA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1 += 'GNG',
                    GNG.append(seq[num])
                else:
                    seq1 += seq[num]
            else:
                seq1 += seq[num]
        else:
            seq1 += seq[num]
    seq1 += seq[-3:]

    shuffle(TNT),shuffle(TNC),shuffle(TNA),shuffle(TNG),shuffle(GNG),shuffle(GNA),shuffle(GNT),shuffle(GNC),shuffle(CNT),shuffle(CNC),shuffle(CNA),shuffle(CNG)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='TNT':seq2+=TNT.pop(0)
        elif seq1[i]=='TNC':seq2+=TNC.pop(0)
        elif seq1[i]=='TNG':seq2+=TNG.pop(0)
        elif seq1[i]=='TNA':seq2+=TNA.pop(0)
        elif seq1[i]=='GNT':seq2+=GNT.pop(0)
        elif seq1[i]=='GNA':seq2+=GNA.pop(0)
        elif seq1[i]=='GNC':seq2+=GNC.pop(0)
        elif seq1[i]=='GNG':seq2+=GNG.pop(0)
        elif seq1[i]=='CNT':seq2+=CNT.pop(0)
        elif seq1[i]=='CNC':seq2+=CNC.pop(0)
        elif seq1[i]=='CNG':seq2+=CNG.pop(0)
        elif seq1[i]=='CNA':seq2+=CNA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    seq1=[]
    THT,THC,THA,THG,GHT,GHA,GHC,GHG,CHG,CHA,CHT,CHC=[],[],[],[],[],[],[],[],[],[],[],[]
    for num in range(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num] in ['T','C','A']:
            if seq[num-2:num] in ['CT','GT','AT']: # ILE3 or LEU4 or VAL
                if seq[num+1]=='T':
                    seq1+='THT',
                    THT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='THC',
                    THC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='THA',
                    THA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='THG',
                    THG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C': # SER4 or PRO or THR or ALA
                if seq[num+1]=='T':
                    seq1+='CHT',
                    CHT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='CHC',
                    CHC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='CHA',
                    CHA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='CHG',
                    CHG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-2:num] in ['CG','GG']: # ARG4 or GLY
                if seq[num+1]=='T':
                    seq1+='GHT',
                    GHT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='GHC',
                    GHC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='GHA',
                    GHA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='GHG',
                    GHG.append(seq[num])
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]

    shuffle(THT),shuffle(THC),shuffle(THA),shuffle(THG),shuffle(GHG),shuffle(GHA),shuffle(GHT),shuffle(GHC),shuffle(CHT),shuffle(CHC),shuffle(CHA),shuffle(CHG)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='THT':seq2+=THT.pop(0)
        elif seq1[i]=='THC':seq2+=THC.pop(0)
        elif seq1[i]=='THA':seq2+=THA.pop(0)
        elif seq1[i]=='THG':seq2+=THG.pop(0)
        elif seq1[i]=='GHT':seq2+=GHT.pop(0)
        elif seq1[i]=='GHA':seq2+=GHA.pop(0)
        elif seq1[i]=='GHC':seq2+=GHC.pop(0)
        elif seq1[i]=='GHG':seq2+=GHG.pop(0)
        elif seq1[i]=='CHT':seq2+=CHT.pop(0)
        elif seq1[i]=='CHC':seq2+=CHC.pop(0)
        elif seq1[i]=='CHG':seq2+=CHG.pop(0)
        elif seq1[i]=='CHA':seq2+=CHA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    seq1=[]
    TRT,TRC,TRA,TRG,ART,ARC,ARG,ARA,GRT,GRA,GRC,GRG,CRG,CRA,CRT,CRC=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for num in range(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num] in ['A','G']:
            if seq[num-1]=='T' and seq[num-2:num]!='AT': # not MET
                if seq[num+1]=='T':
                    seq1+='TRT',
                    TRT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='TRC',
                    TRC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='TRA',
                    TRA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='TRG',
                    TRG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':
                if seq[num+1]=='T':
                    seq1+='CRT',
                    CRT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='CRC',
                    CRC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='CRA',
                    CRA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='CRG',
                    CRG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='A' and seq[num-2:num]!='TA': # not Amber, not Ochre
                if seq[num+1]=='T':
                    seq1+='ART',
                    ART.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='ARC',
                    ARC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='ARA',
                    ARA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='ARG',
                    ARG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='G' and seq[num-2:num]!='TG': # not TRP, not Opal
                if seq[num+1]=='T':
                    seq1+='GRT',
                    GRT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='GRC',
                    GRC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='GRA',
                    GRA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='GRG',
                    GRG.append(seq[num])
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]

    shuffle(TRT),shuffle(TRC),shuffle(TRA),shuffle(TRG),shuffle(GRG),shuffle(GRA),shuffle(GRT),shuffle(GRC),shuffle(ARG),shuffle(ARC),shuffle(ART),shuffle(ARA),shuffle(CRT),shuffle(CRC),shuffle(CRA),shuffle(CRG)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='TRT':seq2+=TRT.pop(0)
        elif seq1[i]=='TRC':seq2+=TRC.pop(0)
        elif seq1[i]=='TRA':seq2+=TRA.pop(0)
        elif seq1[i]=='TRG':seq2+=TRG.pop(0)
        elif seq1[i]=='ART':seq2+=ART.pop(0)
        elif seq1[i]=='ARC':seq2+=ARC.pop(0)
        elif seq1[i]=='ARG':seq2+=ARG.pop(0)
        elif seq1[i]=='ARA':seq2+=ARA.pop(0)
        elif seq1[i]=='GRT':seq2+=GRT.pop(0)
        elif seq1[i]=='GRA':seq2+=GRA.pop(0)
        elif seq1[i]=='GRC':seq2+=GRC.pop(0)
        elif seq1[i]=='GRG':seq2+=GRG.pop(0)
        elif seq1[i]=='CRT':seq2+=CRT.pop(0)
        elif seq1[i]=='CRC':seq2+=CRC.pop(0)
        elif seq1[i]=='CRG':seq2+=CRG.pop(0)
        elif seq1[i]=='CRA':seq2+=CRA.pop(0)
        else:seq2+=seq1[i]
    seq=seq2

    seq1=[]
    TYT,TYC,TYA,TYG,AYT,AYC,AYG,AYA,GYT,GYA,GYC,GYG,CYG,CYA,CYT,CYC=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    for num in range(2,len(seq)-3,3):
        seq1+=seq[num-2:num]
        if seq[num] in ['T','C']:
            if seq[num-1]=='T':
                if seq[num+1]=='T':
                    seq1+='TYT',
                    TYT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='TYC',
                    TYC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='TYA',
                    TYA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='TYG',
                    TYG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='C':
                if seq[num+1]=='T':
                    seq1+='CYT',
                    CYT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='CYC',
                    CYC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='CYA',
                    CYA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='CYG',
                    CYG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='A':
                if seq[num+1]=='T':
                    seq1+='AYT',
                    AYT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='AYC',
                    AYC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='AYA',
                    AYA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='AYG',
                    AYG.append(seq[num])
                else:
                    seq1+=seq[num]
            elif seq[num-1]=='G':
                if seq[num+1]=='T':
                    seq1+='GYT',
                    GYT.append(seq[num])
                elif seq[num+1]=='C':
                    seq1+='GYC',
                    GYC.append(seq[num])
                elif seq[num+1]=='A':
                    seq1+='GYA',
                    GYA.append(seq[num])
                elif seq[num+1]=='G':
                    seq1+='GYG',
                    GYG.append(seq[num])
                else:
                    seq1+=seq[num]
            else:
                seq1+=seq[num]
        else:
            seq1+=seq[num]
    seq1+=seq[-3:]

    shuffle(TYT),shuffle(TYC),shuffle(TYA),shuffle(TYG),shuffle(GYG),shuffle(GYA),shuffle(GYT),shuffle(GYC),shuffle(AYG),shuffle(AYC),shuffle(AYT),shuffle(AYA),shuffle(CYT),shuffle(CYC),shuffle(CYA),shuffle(CYG)
    seq2=''
    for i in range(len(seq1)):
        if seq1[i]=='TYT':seq2+=TYT.pop(0)
        elif seq1[i]=='TYC':seq2+=TYC.pop(0)
        elif seq1[i]=='TYA':seq2+=TYA.pop(0)
        elif seq1[i]=='TYG':seq2+=TYG.pop(0)
        elif seq1[i]=='AYT':seq2+=AYT.pop(0)
        elif seq1[i]=='AYC':seq2+=AYC.pop(0)
        elif seq1[i]=='AYG':seq2+=AYG.pop(0)
        elif seq1[i]=='AYA':seq2+=AYA.pop(0)
        elif seq1[i]=='GYT':seq2+=GYT.pop(0)
        elif seq1[i]=='GYA':seq2+=GYA.pop(0)
        elif seq1[i]=='GYC':seq2+=GYC.pop(0)
        elif seq1[i]=='GYG':seq2+=GYG.pop(0)
        elif seq1[i]=='CYT':seq2+=CYT.pop(0)
        elif seq1[i]=='CYC':seq2+=CYC.pop(0)
        elif seq1[i]=='CYG':seq2+=CYG.pop(0)
        elif seq1[i]=='CYA':seq2+=CYA.pop(0)
        else:seq2+=seq1[i]

    return seq2


def exchange6deg(seq):
    """
    Shuffles the first nucleotide for two six-codon amino acids (LEU and ARG) with the third codon position,
    preserving overall dinucleotide content (but not position-specific dinucleotide content).
    It also handles SER (TCN + AGY) in a more complicated manner.

    After this shuffling, the 'third'-function (or similar) should be used as described in the original code.

    This is a direct Python-3 translation of the old Python-2 code. The code is intentionally
    very large and specialized. Use caution if you plan to modify it.

    Args:
        seq (str): A DNA sequence (typically coding) whose codons will be shuffled according to specific rules.

    Returns:
        str: The shuffled DNA sequence following the exchange6deg procedure.
    """

    # --------------------------------------------------------------------------
    # 1) LEU section
    # --------------------------------------------------------------------------
    seq1 = []
    # LEU has TTR and CTN codons.
    # First conversion of CTY to CTR is required to make them compatible with TTR codons
    # and improve the shuffling efficiency. CTY LEU codons might be re-introduced later
    # upon the third-position shuffling.

    CTYT, CTYA, CTTG, CTCG, CTYC = [], [], [], [], []
    TAT, TAA, TAG, TAC = [], [], [], []
    TGT, TGA, TGG, TGC = [], [], [], []
    TAG_arg = []

    # Add the first two nucleotides of the sequence as-is
    seq1.append(seq[:2])

    # Parse the sequence from index=2 to len(seq)-2
    for num in range(2, len(seq) - 2):
        # Condition 1
        if (num % 3 == 2 and seq[num - 2:num] == 'CT' and seq[num] in ('C', 'T')):
            # We have a codon starting with CT followed by C/T
            if seq[num + 1] == 'A':
                CTYA.append(seq[num])
                seq1.append('CTYA')
            elif seq[num + 1] == 'G':
                if seq[num] == 'C':
                    CTCG.append(seq[num])
                    seq1.append('CTYG')
                else:  # seq[num] == 'T'
                    CTTG.append(seq[num])
                    seq1.append('CTYG')
            elif seq[num + 1] == 'C':
                CTYC.append(seq[num])
                seq1.append('CTYC')
            elif seq[num + 1] == 'T':
                CTYT.append(seq[num])
                seq1.append('CTYT')
            else:
                seq1.append(seq[num])

        # Condition 2
        elif (num % 3 == 2 and (seq[num - 2:num + 1] == 'GTA' or seq[num - 2:num + 1] == 'ATA')):
            # e.g. GTA? / ATA?
            if seq[num + 1] == 'A':
                TAA.append(seq[num])
                seq1.append('TAA')
            elif seq[num + 1] == 'G':
                TAG.append(seq[num])
                seq1.append('TAG')
            elif seq[num + 1] == 'C':
                TAC.append(seq[num])
                seq1.append('TAC')
            elif seq[num + 1] == 'T':
                TAT.append(seq[num])
                seq1.append('TAT')
            else:
                seq1.append(seq[num])

        # Condition 3
        elif (num % 3 == 2 and seq[num - 2:num + 1] == 'GTG'):
            if seq[num + 1] == 'A':
                TGA.append(seq[num])
                seq1.append('TGA')
            elif seq[num + 1] == 'G':
                TGG.append(seq[num])
                seq1.append('TGG')
            elif seq[num + 1] == 'C':
                TGC.append(seq[num])
                seq1.append('TGC')
            elif seq[num + 1] == 'T':
                TGT.append(seq[num])
                seq1.append('TGT')
            else:
                seq1.append(seq[num])

        # Condition 4
        elif (num % 3 == 0 and (seq[num:num + 3] in ('AGA', 'AGG'))
              and seq[num - 1] == 'T' and seq[num - 2:num] != 'TT'
              and seq[num - 3:num] != 'CTT'):
            TAG_arg.append(seq[num])
            seq1.append('TAG_arg')

        else:
            seq1.append(seq[num])

    seq1.append(seq[-2:])  # Add the last two nucleotides

    # Flatten seq1
    seq1 = ''.join(str(x) for x in seq1)

    # Now replacing the third position Y with R in LEU
    # We'll do some partial logic to swap CTCG / TAG_arg, etc.
    CTAG = []
    change_num = min(len(CTCG), len(TAG_arg))
    CTAG = TAG_arg[:change_num]
    TAG_arg = TAG_arg[change_num:]
    tmp = CTCG[:change_num]
    CTCG = CTCG[change_num:]
    # leftover in CTCG
    CTYG = CTCG + CTTG

    shuffle(CTYG)
    change_num = min(len(CTYG), len(TAG))
    CTYG[:change_num], TAG[:change_num] = TAG[:change_num], CTYG[:change_num]
    CTYG.reverse()
    change_num = min(len(CTYG), len(TGG))
    CTYG[:change_num], TGG[:change_num] = TGG[:change_num], CTYG[:change_num]
    CTYG += CTAG

    # More LEU swaps
    change_num = min(len(CTYT), len(TAT))
    CTYT[:change_num], TAT[:change_num] = TAT[:change_num], CTYT[:change_num]
    CTYT.reverse()
    change_num = min(len(CTYT), len(TGT))
    CTYT[:change_num], TGT[:change_num] = TGT[:change_num], CTYT[:change_num]

    change_num = min(len(CTYA), len(TAA))
    CTYA[:change_num], TAA[:change_num] = TAA[:change_num], CTYA[:change_num]
    CTYA.reverse()
    change_num = min(len(CTYA), len(TGA))
    CTYA[:change_num], TGA[:change_num] = TGA[:change_num], CTYA[:change_num]

    change_num = min(len(CTYC), len(TAC))
    CTYC[:change_num], TAC[:change_num] = TAC[:change_num], CTYC[:change_num]
    CTYC.reverse()
    change_num = min(len(CTYC), len(TGC))
    CTYC[:change_num], TGC[:change_num] = TGC[:change_num], CTYC[:change_num]

    # Shuffle them all
    shuffle(CTYT); shuffle(CTYA); shuffle(CTYG); shuffle(CTYC)
    shuffle(TAT);  shuffle(TAA);  shuffle(TAG);  shuffle(TAC)
    shuffle(TGT);  shuffle(TGA);  shuffle(TGG);  shuffle(TGC)
    shuffle(TAG_arg)

    # Replace placeholders
    seq2 = []
    idx = 0
    while idx < len(seq1):
        sub4 = seq1[idx:idx + 4]

        if sub4 == 'CTYT':
            seq2.append(CTYT.pop(0))
            idx += 4
        elif sub4 == 'CTYC':
            seq2.append(CTYC.pop(0))
            idx += 4
        elif sub4 == 'CTYA':
            seq2.append(CTYA.pop(0))
            idx += 4
        elif sub4 == 'CTYG':
            seq2.append(CTYG.pop(0))
            idx += 4
        elif seq1[idx:idx + 3] == 'TAG' and len(TAG) > 0:
            seq2.append(TAG.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TAA':
            seq2.append(TAA.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TAC':
            seq2.append(TAC.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TAT':
            seq2.append(TAT.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TGG':
            seq2.append(TGG.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TGA':
            seq2.append(TGA.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TGC':
            seq2.append(TGC.pop(0))
            idx += 3
        elif seq1[idx:idx + 3] == 'TGT':
            seq2.append(TGT.pop(0))
            idx += 3
        elif seq1[idx:idx + 7] == 'TAG_arg':
            seq2.append(TAG_arg.pop(0))
            idx += 7
        else:
            seq2.append(seq1[idx])
            idx += 1

    seq = ''.join(seq2)

    # --------------------------------------------------------------------------
    # 2) Shuffle the first nucleotide of TTR and CTR LEU codons, etc.
    # --------------------------------------------------------------------------
    seq1 = []
    TYT, AYT, GYT, CYT = [], [], [], []
    seq1.append(seq[:2])

    for num in range(2, len(seq) - 2):
        if (
            (
                (seq[num] in ('T', 'C')) and seq[num + 1] == 'T' and (num % 3 == 2)
                and seq[num + 1:num + 4] not in ('TTA', 'TTG', 'CTA', 'CTG')
            ) or (
                (num % 3 == 0)
                and (seq[num:num + 3] in ('TTA', 'TTG', 'CTA', 'CTG'))
            )
        ):
            if seq[num - 1] == 'T':
                seq1.append('TYT')
                TYT.append(seq[num])
            elif seq[num - 1] == 'C':
                seq1.append('CYT')
                CYT.append(seq[num])
            elif seq[num - 1] == 'A':
                seq1.append('AYT')
                AYT.append(seq[num])
            elif seq[num - 1] == 'G':
                seq1.append('GYT')
                GYT.append(seq[num])
            else:
                seq1.append(seq[num])
        else:
            seq1.append(seq[num])

    seq1.append(seq[-2:])
    seq1 = ''.join(str(x) for x in seq1)

    shuffle(TYT)
    shuffle(GYT)
    shuffle(AYT)
    shuffle(CYT)

    seq2 = []
    idx = 0
    while idx < len(seq1):
        sub3 = seq1[idx:idx + 3]
        if sub3 == 'TYT':
            seq2.append(TYT.pop(0))
            idx += 3
        elif sub3 == 'AYT':
            seq2.append(AYT.pop(0))
            idx += 3
        elif sub3 == 'GYT':
            seq2.append(GYT.pop(0))
            idx += 3
        elif sub3 == 'CYT':
            seq2.append(CYT.pop(0))
            idx += 3
        else:
            seq2.append(seq1[idx])
            idx += 1

    seq = ''.join(seq2)

    # --------------------------------------------------------------------------
    # 3) SER section
    # --------------------------------------------------------------------------
    seq1 = []
    TCRC, TCRA, TCRG, TCRT = [], [], [], []
    CYC, CYA, CYG, CYT = [], [], [], []

    for num in range(2, len(seq) - 2, 3):
        seq1.append(seq[num - 2:num])
        maybe = seq[num - 2:num + 1]
        if maybe in ('TCA', 'TCG'):
            if seq[num + 1] == 'C':
                TCRC.append(seq[num])
                seq1.append('TCRC')
            elif seq[num + 1] == 'G':
                TCRG.append(seq[num])
                seq1.append('TCRG')
            elif seq[num + 1] == 'A':
                TCRA.append(seq[num])
                seq1.append('TCRA')
            elif seq[num + 1] == 'T':
                TCRT.append(seq[num])
                seq1.append('TCRT')
            else:
                seq1.append(seq[num])
        elif (
            seq[num - 2:num] in ('CC', 'AC', 'GC')
            and seq[num] in ('T', 'C')
        ):
            if seq[num + 1] == 'C':
                CYC.append(seq[num])
                seq1.append('CYC')
            elif seq[num + 1] == 'G':
                CYG.append(seq[num])
                seq1.append('CYG')
            elif seq[num + 1] == 'A':
                CYA.append(seq[num])
                seq1.append('CYA')
            elif seq[num + 1] == 'T':
                CYT.append(seq[num])
                seq1.append('CYT')
            else:
                seq1.append(seq[num])
        else:
            seq1.append(seq[num])

    seq1.append(seq[-3:])
    seq1 = ''.join(str(x) for x in seq1)

    # Swapping
    change_num = min(len(TCRC), len(CYC))
    TCRC[:change_num], CYC[:change_num] = CYC[:change_num], TCRC[:change_num]

    change_num = min(len(TCRA), len(CYA))
    TCRA[:change_num], CYA[:change_num] = CYA[:change_num], TCRA[:change_num]

    change_num = min(len(TCRG), len(CYG))
    TCRG[:change_num], CYG[:change_num] = CYG[:change_num], TCRG[:change_num]

    change_num = min(len(TCRT), len(CYT))
    TCRT[:change_num], CYT[:change_num] = CYT[:change_num], TCRT[:change_num]

    shuffle(TCRC); shuffle(TCRA); shuffle(TCRG); shuffle(TCRT)
    shuffle(CYC);  shuffle(CYA);  shuffle(CYG);  shuffle(CYT)

    seq2 = []
    idx = 0
    while idx < len(seq1):
        sub4 = seq1[idx:idx + 4]
        if sub4 == 'TCRC':
            seq2.append(TCRC.pop(0))
            idx += 4
        elif sub4 == 'TCRA':
            seq2.append(TCRA.pop(0))
            idx += 4
        elif sub4 == 'TCRG':
            seq2.append(TCRG.pop(0))
            idx += 4
        elif sub4 == 'TCRT':
            seq2.append(TCRT.pop(0))
            idx += 4
        elif sub4 == 'CYT':
            seq2.append(CYT.pop(0))
            idx += 3
        elif sub4 == 'CYA':
            seq2.append(CYA.pop(0))
            idx += 3
        elif sub4 == 'CYG':
            seq2.append(CYG.pop(0))
            idx += 3
        elif sub4 == 'CYC':
            seq2.append(CYC.pop(0))
            idx += 3
        else:
            seq2.append(seq1[idx])
            idx += 1

    seq = ''.join(seq2)

    # --------------------------------------------------------------------------
    # 4) ATCY, AAGY conversion to BTCY, BAGY
    # --------------------------------------------------------------------------
    seq1 = []
    AAA_R, GAA_N, GAA_R, CAA_N = [], [], [], []
    TAA_R, TAA_N, TAA_H = [], [], []
    AAT_R, GAT_N, GAT_R, CAT_N = [], [], [], []
    TAT_R, TAT_N, TAT_H = [], [], []
    AGA_R, GBA_N, GGA_R = [], [], []
    TGA_R, TBA_N, TYA_H = [], [], []
    AGT_R, GBT_N, GGT_R, CBT_N = [], [], [], []
    TGT_R, TBT_N, TYT_H, CBA_N = [], [], [], []

    for num in range(2, len(seq) - 2, 3):
        seq1.append(seq[num - 2:num])
        four_mer = seq[num:num + 4]

        # AAGC / AAGT
        if four_mer in ('AAGC', 'AAGT'):
            cod2 = seq[num - 2:num]
            if cod2 in ('CA', 'AA', 'GA'):
                AAA_R.append(seq[num])
                seq1.append('AAA_R')
            elif cod2 in ('CG', 'GG'):
                GAA_N.append(seq[num])
                seq1.append('GAA_N')
            elif cod2 == 'AG':
                GAA_R.append(seq[num])
                seq1.append('GAA_R')
            elif cod2 in ('CC', 'AC', 'GC'):
                CAA_N.append(seq[num])
                seq1.append('CAA_N')
            elif cod2 == 'TT':
                TAA_R.append(seq[num])
                seq1.append('TAA_R')
            elif cod2 in ('CT', 'GT'):
                TAA_N.append(seq[num])
                seq1.append('TAA_N')
            elif cod2 == 'AT':
                TAA_H.append(seq[num])
                seq1.append('TAA_H')
            else:
                seq1.append(seq[num])

        # ATCC / ATCT
        elif four_mer in ('ATCC', 'ATCT'):
            cod2 = seq[num - 2:num]
            if cod2 in ('CA', 'AA', 'GA'):
                AAT_R.append(seq[num])
                seq1.append('AAT_R')
            elif cod2 in ('CG', 'GG'):
                GAT_N.append(seq[num])
                seq1.append('GAT_N')
            elif cod2 == 'AG':
                GAT_R.append(seq[num])
                seq1.append('GAT_R')
            elif cod2 in ('CC', 'AC', 'GC'):
                CAT_N.append(seq[num])
                seq1.append('CAT_N')
            elif cod2 == 'TT':
                TAT_R.append(seq[num])
                seq1.append('TAT_R')
            elif cod2 in ('CT', 'GT'):
                TAT_N.append(seq[num])
                seq1.append('TAT_N')
            elif cod2 == 'AT':
                TAT_H.append(seq[num])
                seq1.append('TAT_H')
            else:
                seq1.append(seq[num])

        # AGA / AGG region
        elif (
            seq[num + 1:num + 4] not in ('AGC', 'AGT')
            and seq[num + 1] == 'A'
            and seq[num] != 'A'
        ):
            cod3 = seq[num - 2:num + 1]
            cod2 = seq[num - 2:num]
            if cod3 in ('CAG', 'AAG', 'GAG'):
                AGA_R.append(seq[num])
                seq1.append('AGA_R')
            elif cod2 in ('CG', 'GG'):
                GBA_N.append(seq[num])
                seq1.append('GBA_N')
            elif cod3 == 'AGG':
                GGA_R.append(seq[num])
                seq1.append('GGA_R')
            elif cod2 in ('CC', 'AC', 'GC'):
                CBA_N.append(seq[num])
                seq1.append('CBA_N')
            elif cod3 == 'TTG':
                TGA_R.append(seq[num])
                seq1.append('TGA_R')
            elif cod2 in ('CT', 'GT'):
                TBA_N.append(seq[num])
                seq1.append('TBA_N')
            elif cod2 == 'AT' and seq[num] != 'G':
                TYA_H.append(seq[num])
                seq1.append('TYA_H')
            else:
                seq1.append(seq[num])

        # T... region
        elif (
            seq[num + 1:num + 4] not in ('TCC', 'TCT')
            and seq[num + 1] == 'T'
            and seq[num] != 'A'
        ):
            cod3 = seq[num - 2:num + 1]
            cod2 = seq[num - 2:num]
            if cod3 in ('CAG', 'AAG', 'GAG'):
                AGT_R.append(seq[num])
                seq1.append('AGT_R')
            elif cod2 in ('CG', 'GG'):
                GBT_N.append(seq[num])
                seq1.append('GBT_N')
            elif cod3 == 'AGG':
                GGT_R.append(seq[num])
                seq1.append('GGT_R')
            elif cod2 in ('CC', 'AC', 'GC'):
                CBT_N.append(seq[num])
                seq1.append('CBT_N')
            elif cod3 == 'TTG':
                TGT_R.append(seq[num])
                seq1.append('TGT_R')
            elif cod2 in ('CT', 'GT'):
                TBT_N.append(seq[num])
                seq1.append('TBT_N')
            elif cod2 == 'AT' and seq[num] != 'G':
                TYT_H.append(seq[num])
                seq1.append('TYT_H')
            else:
                seq1.append(seq[num])
        else:
            seq1.append(seq[num])

    seq1.append(seq[-3:])
    seq1 = ''.join(str(x) for x in seq1)

    def swap_min(a, b):
        n = min(len(a), len(b))
        a[:n], b[:n] = b[:n], a[:n]

    swap_min(AAA_R, AGA_R)
    swap_min(GAA_R, GGA_R)
    swap_min(TAA_R, TGA_R)
    swap_min(AAT_R, AGT_R)
    swap_min(GAT_R, GGT_R)
    swap_min(TAT_R, TGT_R)
    swap_min(TAT_H, TYT_H)
    swap_min(TAA_H, TYA_H)
    swap_min(CAA_N, CBA_N)
    swap_min(GAA_N, GBA_N)
    swap_min(TAA_N, TBA_N)
    swap_min(GAT_N, GBT_N)
    swap_min(CAT_N, CBT_N)
    swap_min(TAT_N, TBT_N)

    all_lists = [
        AAA_R, GAA_N, GAA_R, CAA_N, TAA_R, TAA_N, TAA_H,
        AAT_R, GAT_N, GAT_R, CAT_N, TAT_R, TAT_N, TAT_H,
        AGA_R, GBA_N, GGA_R, TGA_R, TBA_N, TYA_H,
        AGT_R, GBT_N, GGT_R, CBT_N, TGT_R, TBT_N, TYT_H, CBA_N
    ]
    for lst in all_lists:
        shuffle(lst)

    seq2 = []
    idx = 0
    while idx < len(seq1):
        chunk5 = seq1[idx:idx + 5]

        replaced = False
        if chunk5.startswith('AAA_R'):
            seq2.append(AAA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GAA_N'):
            seq2.append(GAA_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GAA_R'):
            seq2.append(GAA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('CAA_N'):
            seq2.append(CAA_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAA_R'):
            seq2.append(TAA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAA_N'):
            seq2.append(TAA_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAA_H'):
            seq2.append(TAA_H.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('AAT_R'):
            seq2.append(AAT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GAT_N'):
            seq2.append(GAT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GAT_R'):
            seq2.append(GAT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('CAT_N'):
            seq2.append(CAT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAT_R'):
            seq2.append(TAT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAT_N'):
            seq2.append(TAT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TAT_H'):
            seq2.append(TAT_H.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('AGA_R'):
            seq2.append(AGA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GBA_N'):
            seq2.append(GBA_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GGA_R'):
            seq2.append(GGA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TGA_R'):
            seq2.append(TGA_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TBA_N'):
            seq2.append(TBA_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TYA_H'):
            seq2.append(TYA_H.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('AGT_R'):
            seq2.append(AGT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GBT_N'):
            seq2.append(GBT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('GGT_R'):
            seq2.append(GGT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('CBT_N'):
            seq2.append(CBT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TGT_R'):
            seq2.append(TGT_R.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TBT_N'):
            seq2.append(TBT_N.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('TYT_H'):
            seq2.append(TYT_H.pop(0))
            idx += 5
            replaced = True
        elif chunk5.startswith('CBA_N'):
            seq2.append(CBA_N.pop(0))
            idx += 5
            replaced = True
        else:
            seq2.append(seq1[idx])
            idx += 1

    seq = ''.join(seq2)

    # --------------------------------------------------------------------------
    # 5) G_SER (some further conversions)
    # --------------------------------------------------------------------------
    seq1 = []
    GTG, GAG = [], []
    CTG, CAG = [], []
    TTG, TAG = [], []
    TCC, TGC = [], []
    TCT, TGT = [], []
    GAGC, CAGC, TAGC = [], [], []
    GAGT, CAGT, TAGT = [], [], []
    GTCC, CTCC, TTCC = [], [], []
    GTCT, CTCT, TTCT = [], [], []

    i = 2
    while i < len(seq) - 3:
        chunk_2 = seq[i - 2:i]
        chunk_3 = seq[i - 3:i]
        # Big checks, as in original code
        if (chunk_2 == 'AG' and seq[i] in ('C', 'T') and chunk_3[-1] in ('G', 'C', 'T')):
            # GAGC, GAGT, ...
            if chunk_3[-1] == 'G':
                if seq[i] == 'C':
                    seq1.append('GAGC')
                    seq1.append(seq[i + 1])
                    GAGC.append(chunk_2)
                else:
                    seq1.append('GAGT')
                    seq1.append(seq[i + 1])
                    GAGT.append(chunk_2)
            elif chunk_3[-1] == 'C':
                if seq[i] == 'C':
                    seq1.append('CAGC')
                    seq1.append(seq[i + 1])
                    CAGC.append(chunk_2)
                else:
                    seq1.append('CAGT')
                    seq1.append(seq[i + 1])
                    CAGT.append(chunk_2)
            elif chunk_3[-1] == 'T':
                if seq[i] == 'C':
                    seq1.append('TAGC')
                    seq1.append(seq[i + 1])
                    TAGC.append(chunk_2)
                else:
                    seq1.append('TAGT')
                    seq1.append(seq[i + 1])
                    TAGT.append(chunk_2)
            i += 2

        elif (chunk_2 == 'TC' and seq[i] in ('C', 'T') and chunk_3[-1] in ('G', 'C', 'T')):
            if chunk_3[-1] == 'G':
                if seq[i] == 'C':
                    seq1.append('GTCC')
                    seq1.append(seq[i + 1])
                    GTCC.append(chunk_2)
                else:
                    seq1.append('GTCT')
                    seq1.append(seq[i + 1])
                    GTCT.append(chunk_2)
            elif chunk_3[-1] == 'C':
                if seq[i] == 'C':
                    seq1.append('CTCC')
                    seq1.append(seq[i + 1])
                    CTCC.append(chunk_2)
                else:
                    seq1.append('CTCT')
                    seq1.append(seq[i + 1])
                    CTCT.append(chunk_2)
            elif chunk_3[-1] == 'T':
                if seq[i] == 'C':
                    seq1.append('TTCC')
                    seq1.append(seq[i + 1])
                    TTCC.append(chunk_2)
                else:
                    seq1.append('TTCT')
                    seq1.append(seq[i + 1])
                    TTCT.append(chunk_2)
            i += 2

        elif (seq[i - 1:i + 2] == 'GTG' and seq[i - 2] != 'C' and i % 3 == 2):
            seq1.append(seq[i - 2:i])
            seq1.append('GTG')
            GTG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'GAG' and seq[i - 2] != 'C' and i % 3 == 2):
            seq1.append(seq[i - 2:i])
            seq1.append('GAG')
            GAG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'CTG' and chunk_2 in ('AC', 'GC', 'CC')):
            seq1.append(chunk_2)
            seq1.append('CTG')
            CTG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'CAG' and chunk_2 in ('AC', 'GC', 'CC')):
            seq1.append(chunk_2)
            seq1.append('CAG')
            CAG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'TTG' and chunk_2 in ('AT', 'GT', 'CT')):
            seq1.append(chunk_2)
            seq1.append('TTG')
            TTG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'TAG' and chunk_2 in ('AT', 'GT', 'CT')):
            seq1.append(chunk_2)
            seq1.append('TAG')
            TAG.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'TCC' and chunk_2 in ('GT', 'CT')):
            seq1.append(chunk_2)
            seq1.append('TCC')
            TCC.append(seq[i])
            i += 2
        elif (seq[i - 1:i + 2] == 'TGC' and chunk_2 in ('GT', 'CT')):
            seq1.append(chunk_2)
            seq1.append('TGC')
            TGC.append(seq[i])
            i += 2
        elif (
            seq[i - 1:i + 2] == 'TCT'
            and chunk_2 in ('GT', 'CT')
            and seq[i + 1:i + 4] not in ('TCC', 'TCT')
        ):
            seq1.append(chunk_2)
            seq1.append('TCT')
            TCT.append(seq[i])
            i += 2
        elif (
            seq[i - 1:i + 2] == 'TGT'
            and chunk_2 in ('GT', 'CT')
            and seq[i + 1:i + 4] not in ('TCC', 'TCT')
        ):
            seq1.append(chunk_2)
            seq1.append('TGT')
            TGT.append(seq[i])
            i += 2
        else:
            seq1.append(seq[i - 2:i + 1])
            i += 1

    seq1.append(seq[i - 2:])
    seq1 = ''.join(str(x) for x in seq1)

    # Next, define the min-lens for the swapping logic
    gAGc = min(len(GAGC), len(GTG), len(TCC))
    gTCc = min(len(GTCC), len(GAG), len(TGC))

    gAGt = min(len(GAGT), len(GTG), len(TCT))
    gTCt = min(len(GTCT), len(GAG), len(TGT))

    cAGc = min(len(CAGC), len(CTG), len(TCC))
    cTCc = min(len(CTCC), len(CAG), len(TGC))

    cAGt = min(len(CAGT), len(CTG), len(TCT))
    cTCt = min(len(CTCT), len(CAG), len(TGT))

    tAGc = min(len(TAGC), len(TTG), len(TCC))
    tTCc = min(len(TTCC), len(TAG), len(TGC))

    tAGt = min(len(TAGC), len(TTG), len(TCT))
    tTCt = min(len(TTCC), len(TAG), len(TGT))

    # Adjust distributions if the sum is larger than the length of the target list
    if gAGc + gAGt > len(GTG):
        ratio = len(GTG) / float(gAGc + gAGt)
        gAGc = int(round(gAGc * ratio))
        gAGt = int(round(gAGt * ratio))
    if gTCc + gTCt > len(GAG):
        ratio = len(GAG) / float(gTCc + gTCt)
        gTCc = int(round(gTCc * ratio))
        gTCt = int(round(gTCt * ratio))

    if cAGc + cAGt > len(CTG):
        ratio = len(CTG) / float(cAGc + cAGt)
        cAGc = int(round(cAGc * ratio))
        cAGt = int(round(cAGt * ratio))
    if cTCc + cTCt > len(CAG):
        ratio = len(CAG) / float(cTCc + cTCt)
        cTCc = int(round(cTCc * ratio))
        cTCt = int(round(cTCt * ratio))

    if tAGc + tAGt > len(TTG):
        ratio = len(TTG) / float(tAGc + tAGt)
        tAGc = int(round(tAGc * ratio))
        tAGt = int(round(tAGt * ratio))
    if tTCc + tTCt > len(TAG):
        ratio = len(TAG) / float(tTCc + tTCt)
        tTCc = int(round(tTCc * ratio))
        tTCt = int(round(tTCt * ratio))

    if (gAGc + cAGc + tAGc) > len(TCC):
        s = float(gAGc + cAGc + tAGc)
        ratio = len(TCC) / s
        gAGc = int(round(gAGc * ratio))
        cAGc = int(round(cAGc * ratio))
        tAGc = int(round(tAGc * ratio))

    if (gAGt + cAGt + tAGt) > len(TCT):
        s = float(gAGt + cAGt + tAGt)
        ratio = len(TCT) / s
        gAGt = int(round(gAGt * ratio))
        cAGt = int(round(cAGt * ratio))
        tAGt = int(round(tAGt * ratio))

    if (gTCc + cTCc + tTCc) > len(TGC):
        s = float(gTCc + cTCc + tTCc)
        ratio = len(TGC) / s
        gTCc = int(round(gTCc * ratio))
        cTCc = int(round(cTCc * ratio))
        tTCc = int(round(tTCc * ratio))

    if (gTCt + cTCt + tTCt) > len(TGT):
        s = float(gTCt + cTCt + tTCt)
        ratio = len(TGT) / s
        gTCt = int(round(gTCt * ratio))
        cTCt = int(round(cTCt * ratio))
        tTCt = int(round(tTCt * ratio))

    # Now do random partial swaps
    change_num = randint(0, gAGc)
    for _ in range(change_num):
        GAGC[0] = 'TC'
        GTG.pop(0); GTG.append('A')
        TCC.pop(0); TCC.append('G')

    change_num = randint(0, gTCc)
    for _ in range(change_num):
        GTCC[0] = 'AG'
        GAG.pop(0); GAG.append('T')
        TGC.pop(0); TGC.append('C')

    change_num = randint(0, gAGt)
    for _ in range(change_num):
        GAGT[0] = 'TC'
        GTG.pop(0); GTG.append('A')
        TCT.pop(0); TCT.append('G')

    change_num = randint(0, gTCt)
    for _ in range(change_num):
        GTCT[0] = 'AG'
        GAG.pop(0); GAG.append('T')
        TGT.pop(0); TGT.append('C')

    change_num = randint(0, cAGc)
    for _ in range(change_num):
        CAGC[0] = 'TC'
        CTG.pop(0); CTG.append('A')
        TCC.pop(0); TCC.append('G')

    change_num = randint(0, cTCc)
    for _ in range(change_num):
        CTCC[0] = 'AG'
        CAG.pop(0); CAG.append('T')
        TGC.pop(0); TGC.append('C')

    change_num = randint(0, cAGt)
    for _ in range(change_num):
        CAGT[0] = 'TC'
        CTG.pop(0); CTG.append('A')
        TCT.pop(0); TCT.append('G')

    change_num = randint(0, cTCt)
    for _ in range(change_num):
        CTCT[0] = 'AG'
        CAG.pop(0); CAG.append('T')
        TGT.pop(0); TGT.append('C')

    change_num = randint(0, tAGc)
    for _ in range(change_num):
        TAGC[0] = 'TC'
        TTG.pop(0); TTG.append('A')
        TCC.pop(0); TCC.append('G')

    change_num = randint(0, tTCc)
    for _ in range(change_num):
        TTCC[0] = 'AG'
        TAG.pop(0); TAG.append('T')
        TGC.pop(0); TGC.append('C')

    change_num = randint(0, tAGt)
    for _ in range(change_num):
        TAGC[0] = 'TC'  # Original code reuses TAGC, potentially a logic quirk
        TTG.pop(0); TTG.append('A')
        TCT.pop(0); TCT.append('G')

    change_num = randint(0, tTCt)
    for _ in range(change_num):
        TTCC[0] = 'AG'
        TAG.pop(0); TAG.append('T')
        TGT.pop(0); TGT.append('C')

    # Shuffle all involved lists
    shuffle(GTG); shuffle(GAG); shuffle(CTG); shuffle(CAG)
    shuffle(TTG); shuffle(TAG); shuffle(TCC); shuffle(TGC)
    shuffle(TCT); shuffle(TGT)
    shuffle(GAGC); shuffle(CAGC); shuffle(TAGC); shuffle(GAGT)
    shuffle(CAGT); shuffle(TAGT)
    shuffle(GTCC); shuffle(CTCC); shuffle(TTCC); shuffle(GTCT)
    shuffle(CTCT); shuffle(TTCT)

    seq2 = []
    i = 0
    while i < len(seq1):
        # Replace placeholders with shuffled items
        if seq1[i:i+3] == 'GTG':
            seq2.append(GTG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'GAG':
            seq2.append(GAG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'CTG':
            seq2.append(CTG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'CAG':
            seq2.append(CAG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TTG':
            seq2.append(TTG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TAG':
            seq2.append(TAG.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TCC':
            seq2.append(TCC.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TGC':
            seq2.append(TGC.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TCT':
            seq2.append(TCT.pop(0))
            i += 3
        elif seq1[i:i+3] == 'TGT':
            seq2.append(TGT.pop(0))
            i += 3
        elif seq1[i:i+4] == 'GAGC':
            seq2.append(GAGC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'CAGC':
            seq2.append(CAGC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'TAGC':
            seq2.append(TAGC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'GAGT':
            seq2.append(GAGT.pop(0))
            i += 4
        elif seq1[i:i+4] == 'CAGT':
            seq2.append(CAGT.pop(0))
            i += 4
        elif seq1[i:i+4] == 'TAGT':
            seq2.append(TAGT.pop(0))
            i += 4
        elif seq1[i:i+4] == 'GTCC':
            seq2.append(GTCC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'CTCC':
            seq2.append(CTCC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'TTCC':
            seq2.append(TTCC.pop(0))
            i += 4
        elif seq1[i:i+4] == 'GTCT':
            seq2.append(GTCT.pop(0))
            i += 4
        elif seq1[i:i+4] == 'CTCT':
            seq2.append(CTCT.pop(0))
            i += 4
        elif seq1[i:i+4] == 'TTCT':
            seq2.append(TTCT.pop(0))
            i += 4
        else:
            seq2.append(seq1[i])
            i += 1

    seq = ''.join(seq2)

    # ARG
    # Shuffling the first codon position of arginine codons requires conversion
    # of CGY to CGR. Then the first codon positions of AGR and CGR are shuffled
    # with the third codon position of other codons. Overall dinucleotide content
    # and amino acid sequence are maintained.

    seq1 = []
    GYA, GYT, GYC, GYG = [], [], [], []
    GRA = []  # only gly
    GRT, GRC, GRG = [], [], []

    for num in range(2, len(seq) - 3, 3):
        seq1.append(seq[num - 2:num])
        if seq[num - 2:num] == 'CG' and seq[num] in ('C', 'T'):
            # CGY
            if seq[num + 1] == 'A':
                GYA.append(seq[num])
                seq1.append('GYA')
            elif seq[num + 1] == 'T':
                GYT.append(seq[num])
                seq1.append('GYT')
            elif seq[num + 1] == 'C':
                GYC.append(seq[num])
                seq1.append('GYC')
            elif seq[num + 1] == 'G':
                GYG.append(seq[num])
                seq1.append('GYG')
            else:
                seq1.append(seq[num])
        elif (seq[num - 2:num] == 'GG'
              and seq[num] in ('A', 'G')):  # only gly
            if seq[num + 1] == 'A':
                GRA.append(seq[num])
                seq1.append('GRA')
            elif seq[num + 1] == 'T':
                GRT.append(seq[num])
                seq1.append('GRT')
            elif seq[num + 1] == 'C':
                GRC.append(seq[num])
                seq1.append('GRC')
            elif seq[num + 1] == 'G':
                GRG.append(seq[num])
                seq1.append('GRG')
            else:
                seq1.append(seq[num])
        else:
            seq1.append(seq[num])

    seq1.append(seq[-3:])
    seq1 = ''.join(str(x) for x in seq1)

    change_num = min(len(GYA), len(GRA))
    GYA[:change_num], GRA[:change_num] = GRA[:change_num], GYA[:change_num]

    change_num = min(len(GYT), len(GRT))
    GYT[:change_num], GRT[:change_num] = GRT[:change_num], GYT[:change_num]

    change_num = min(len(GYC), len(GRC))
    GYC[:change_num], GRC[:change_num] = GRC[:change_num], GYC[:change_num]

    change_num = min(len(GYG), len(GRG))
    GYG[:change_num], GRG[:change_num] = GRG[:change_num], GYG[:change_num]

    shuffle(GYA); shuffle(GYT); shuffle(GYC); shuffle(GYG)
    shuffle(GRA); shuffle(GRG); shuffle(GRC); shuffle(GRT)

    seq2_list = []
    i = 0
    while i < len(seq1):
        placeholder = seq1[i:i+3]
        if placeholder == 'GYA':
            seq2_list.append(GYA.pop(0))
            i += 3
        elif placeholder == 'GRA':
            seq2_list.append(GRA.pop(0))
            i += 3
        elif placeholder == 'GYT':
            seq2_list.append(GYT.pop(0))
            i += 3
        elif placeholder == 'GRT':
            seq2_list.append(GRT.pop(0))
            i += 3
        elif placeholder == 'GYC':
            seq2_list.append(GYC.pop(0))
            i += 3
        elif placeholder == 'GRC':
            seq2_list.append(GRC.pop(0))
            i += 3
        elif placeholder == 'GYG':
            seq2_list.append(GYG.pop(0))
            i += 3
        elif placeholder == 'GRG':
            seq2_list.append(GRG.pop(0))
            i += 3
        else:
            seq2_list.append(seq1[i])
            i += 1

    seq = ''.join(seq2_list)

    # Next mini-block
    seq1 = []
    CGTG, CGCG, GTGp, GCGp, GAGp = [], [], [], [], []
    seq1.append(seq[:2])

    for num in range(2, len(seq) - 2):
        four_sub = seq[num - 2:num + 2]
        if four_sub == 'CGCG' and num % 3 == 2:
            CGCG.append(seq[num])
            seq1.append('CGYG')
        elif four_sub == 'CGTG' and num % 3 == 2:
            CGTG.append(seq[num])
            seq1.append('CGYG')
        elif seq[num - 1:num + 2] == 'GTG' and seq[num - 2] != 'C' and num % 3 == 2:
            GTGp.append(seq[num])
            seq1.append('GTG')
        elif seq[num - 1:num + 2] == 'GCG' and seq[num - 2] != 'C' and num % 3 == 2:
            GCGp.append(seq[num])
            seq1.append('GCG')
        elif (seq[num - 1:num + 3] in ('GAGA', 'GAGG') and num % 3 == 0):
            GAGp.append(seq[num])
            seq1.append('GAG')
        else:
            seq1.append(seq[num])

    seq1.append(seq[-2:])
    seq1 = ''.join(str(x) for x in seq1)

    CGYG = []
    change_num = min(len(CGTG), len(GCGp))
    # We swap some items, then unify them into CGYG
    for _ in range(change_num):
        if CGTG and GCGp:
            v = CGTG.pop(0)
            w = GCGp.pop(0)
            CGYG.append(w)
            CGYG.append(v)

    # Add CGCG to CGYG
    CGYG.extend(CGCG)

    change_num2 = min(len(CGYG), len(GAGp))
    for _ in range(change_num2):
        if CGYG and GAGp:
            v2 = CGYG.pop(0)
            w2 = GAGp.pop(0)
            CGYG.append(w2)
            GAGp.append(v2)

    # Finally extend with leftover CGTG if needed
    # (The original code does something like CGYG += CGTG[CGTG_num:])
    # We'll approximate it:
    CGYG.extend(CGTG)

    shuffle(CGYG)
    shuffle(GTGp)
    shuffle(GCGp)
    shuffle(GAGp)

    seq2_list = []
    i = 0
    while i < len(seq1):
        chunk4 = seq1[i:i+4]
        chunk3 = seq1[i:i+3]
        if chunk4 == 'CGYG':
            if CGYG:
                seq2_list.append(CGYG.pop(0))
            i += 4
        elif chunk3 == 'GTG':
            if GTGp:
                seq2_list.append(GTGp.pop(0))
            i += 3
        elif chunk3 == 'GCG':
            if GCGp:
                seq2_list.append(GCGp.pop(0))
            i += 3
        elif chunk3 == 'GAG':
            if GAGp:
                seq2_list.append(GAGp.pop(0))
            i += 3
        else:
            seq2_list.append(seq1[i])
            i += 1

    seq = ''.join(seq2_list)

    # Last section
    seq1 = []
    TMG, AMG, GMG, CMG = [], [], [], []
    # The code appends the first two nucleotides
    seq1.append(seq[:2])

    for num in range(2, len(seq) - 2):
        # If we see CG/AG patterns at certain positions, we do placeholders
        condition1 = (num % 3 == 0 and seq[num:num+3] in ('CGA', 'CGG', 'AGA', 'AGG'))
        condition2 = (num % 3 == 2 and seq[num:num+2] in ('CG', 'AG')
                      and ((seq[num-1] == 'T' and seq[num-2] != 'T')
                           or seq[num-1] == 'C'
                           or seq[num-2:num] == 'GG'))
        if condition1 or condition2:
            if seq[num-1] == 'T':
                seq1.append('TMG')
                TMG.append(seq[num])
            elif seq[num-1] == 'C':
                seq1.append('CMG')
                CMG.append(seq[num])
            elif seq[num-1] == 'A':
                seq1.append('AMG')
                AMG.append(seq[num])
            elif seq[num-1] == 'G':
                seq1.append('GMG')
                GMG.append(seq[num])
            else:
                seq1.append(seq[num])
        else:
            seq1.append(seq[num])

    seq1.append(seq[-2:])
    # Convert seq1 from list to string
    seq1 = ''.join(str(x) for x in seq1)

    shuffle(TMG)
    shuffle(GMG)
    shuffle(AMG)
    shuffle(CMG)

    seq2_list = []
    i = 0
    while i < len(seq1):
        placeholder = seq1[i:i+3]
        if placeholder == 'TMG':
            seq2_list.append(TMG.pop(0))
            i += 3
        elif placeholder == 'AMG':
            seq2_list.append(AMG.pop(0))
            i += 3
        elif placeholder == 'GMG':
            seq2_list.append(GMG.pop(0))
            i += 3
        elif placeholder == 'CMG':
            seq2_list.append(CMG.pop(0))
            i += 3
        else:
            seq2_list.append(seq1[i])
            i += 1

    return ''.join(seq2_list)


def get_difference(seq1,seq2):
    assert len(seq1) == len(seq2)
    return sum(seq1c != seq2c for seq1c, seq2c in zip(seq1,seq2))

def make_protein_record(nuc_record):
    return SeqRecord(seq = nuc_record.seq.translate(to_stop=True),
                     id = "trans_" + nuc_record.id,
                     description = "translation of CDS, using default table")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CodonShuffle.')
    parser.add_argument('-i', nargs='?', help='Input Filename', required=True, dest="input_file_name")
    parser.add_argument('-s', choices=['dn23', 'dn31', 'dn231', 'n3','gc3'], nargs='?', help='Type of shuffle', default="dn23", dest="random_type")
    parser.add_argument('-r', nargs='?', help='Number of replications (int)', default='1000', dest="reps", type=int)
    parser.add_argument('-m', choices=['CAI', 'CPB', 'DN', 'ENC', 'VFOLD', 'UFOLD', 'all'], nargs='*', help='Control Features', default='all', dest="modules")
    parser.add_argument('-g', dest="graphics", help='Generate Feature Plots', action="store_true")
    parser.add_argument('--seed', type=int, nargs='?', dest='randomseed', help='Optional integer for random seed', const=99)
    args = parser.parse_args()

    if args.randomseed is not None:
        seed(args.randomseed)

    types_of_rnd = args.random_type

    infile = open(args.input_file_name,'r')
    names_list = []
    data = {}
    for line in infile:
        if line[0]=='>':
            strain=line
            data[strain]=''
            names_list.append(strain)
        else:
            for liter in line:
                if liter not in ['\n','-','~']:
                    data[strain]+=liter
    infile.close()

    for strain in names_list:
        seq_name=''
        for liter in strain:
            if liter in ['\\','/',' ','-',',','|',':']:
                seq_name+='_'
            else:
                seq_name+=liter
        inseq_file=open(seq_name[1:-1]+'.fas','w')
        inseq_file.write(strain+data[strain]+'\n')
        inseq_file.close()

        outfile=open(seq_name[1:-1]+'_'+args.random_type+'.fas','w')
        outfile.write(strain+data[strain]+'\n')
        outfile.close()

        outfile=open(seq_name[1:-1]+'_'+args.random_type+'.fas','a')
        for i in range(args.reps):
            outseq=data[strain]
            if args.random_type =='gc3':
                outseq=gc3(outseq)
            elif args.random_type == 'dn23':
                outseq=dn23(outseq)
            elif args.random_type=='dn31':
                outseq=third(outseq)
            elif args.random_type=='dn231':
                outseq=third(exchange6deg(outseq))
            elif args.random_type=='n3':
                outseq=third_simple(outseq)
            outfile.write('>replicate_'+str(i+1)+'\n'+outseq+'\n')
        outfile.close()

    final_table = pd.DataFrame()

    filename = seq_name[1:-1]+'_'+args.random_type+'.fas'

    print("Calculating Hamming distance")
    seq_records=[]
    # 'rU' is deprecated, use 'r' instead
    inseq_file=open(filename, newline='')
    outnt_file=open(filename+'.hamming', 'w')
    for seq_record in SeqIO.parse(inseq_file, "fasta"):
        seq_records.append(seq_record.seq)
    inseq_file.close()

    n = len(seq_records)
    my_array = np.zeros((n,n))
    for i in range(0,n):
        for j in range(i+1,n):
            difference = get_difference(seq_records[i], seq_records[j])
            outnt_file.write(str(difference)+'\n')
            my_array[i, j] = difference
            my_array[j, i] = difference
    outnt_file.close()

    if args.graphics:
        hamming_graphname = filename+'.hamming.pdf'
        hamming_table = pd.read_csv(filename+'.hamming', sep='\t', names=['distance'])
        hamming_graph = (ggplot(hamming_table, aes('distance'))
                         + geom_density()
                         + labs(x="Hamming distance", y="Frequency")
                         + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (Hamming)'))
        ggsave(hamming_graph, filename=hamming_graphname)

    least_squares = pd.DataFrame(np.zeros(n))

    # CAI
    if 'CAI' in args.modules or 'all' in args.modules:
        print("Calculating CAI")
        cainame = seq_name[1:-1]+'_'+args.random_type+'.cai'
        bulkname = seq_name[1:-1]+'_'+args.random_type+'.bulk'
        
        call([
            "cai",
            filename,
            cainame,
            bulkname,
            "-cai",
            "-cfile", "Ehuman.cut",
            "-nomenu",
            "-nowarn",
            "-silent"
        ])
        
        cai_value = None
        with open(cainame, 'r') as f:
            for line in f:
                if "Codon Adaptation Index" in line:
                    parts = line.strip().split()
                    try:
                        cai_value = float(parts[-1])
                    except ValueError:
                        print(f"Error parsing CAI value: {parts[-1]}")
                        sys.exit(1)
                    break
        
        if cai_value is None:
            print(f"CAI value not found in {cainame}")
            sys.exit(1)
        
        cai_table = pd.Series([cai_value])
        cai_table_z =  ((cai_table - cai_table.mean()) / cai_table.std())
        cai_table_wt_z = cai_table_z.iloc[0]
        cai_table_wt_z_rep = np.repeat(cai_table_wt_z, len(cai_table_z))
        cai_table_z_ls = (cai_table_wt_z_rep - cai_table_z)**2
        least_squares = least_squares.add(cai_table_z_ls, axis=0)
        final_table.insert(0, "Cai", cai_table)
        
        if args.graphics:
            cai_graphname = cainame +'.pdf'
            cai_graph = (
                ggplot(pd.DataFrame({'cai': cai_table}), aes('cai'))
                + geom_density()
                + labs(x="CAI", y="Frequency")
                + geom_vline(xintercept=[cai_table.iloc[0]], colour="red", linetype="dashed")
                + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (CAI)')
            )
            ggsave(cai_graph, filename=cai_graphname)

    # ENC
    if 'ENC' in args.modules or 'all' in args.modules:
        # Run ENC calculation using codonw
        print("Calculating ENC using codonw")
        # Modify the path to codonw as needed, or ensure it's in your PATH
        call(["codonw", filename, "-enc", "-nomenu", "-nowarn", "-silent"])

        enc_filename = seq_name[1:-1] + '_' + args.random_type + '.out'
        enc_table = pd.read_csv(enc_filename, sep='\t')

        # Drop unwanted column if it exists
        if 'Unnamed: 2' in enc_table.columns:
            enc_table = enc_table.drop('Unnamed: 2', axis=1)

        # Z-score normalization of ENC values
        enc_table_z = (enc_table['Nc'] - enc_table['Nc'].mean()) / enc_table['Nc'].std()
        enc_table_wt_z = enc_table_z.iloc[0]
        enc_table_wt_z_rep = np.repeat(enc_table_wt_z, len(enc_table_z))
        enc_table_z_ls = (enc_table_wt_z_rep - enc_table_z)**2

        # Update least_squares DataFrame
        least_squares = least_squares.add(enc_table_z_ls, axis=0)

        # Insert ENC values into final_table
        final_table.insert(0, "ENC", enc_table['Nc'])

        # If graphics option is enabled, generate a density plot
        if args.graphics:
            print("Generating ENC density plot")
            # Re-read the enc_table if needed (not strictly necessary if enc_table is unchanged)
            enc_table = pd.read_csv(enc_filename, sep='\t')
            if 'Unnamed: 2' in enc_table.columns:
                enc_table = enc_table.drop('Unnamed: 2', axis=1)

            enc_graphname = enc_filename + '.enc.pdf'
            enc_graph = (
                ggplot(enc_table, aes('Nc'))
                + geom_density()
                + labs(x="ENC", y="Frequency")
                + geom_vline(xintercept=[enc_table['Nc'].iloc[0]], colour="red", linetype="dashed")
                + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (ENC)')
            )
            ggsave(enc_graph, filename=enc_graphname)


    if 'VFOLD' in args.modules or 'all' in args.modules:
        print("Calculating VFOLD using ViennaRNA Python interface")

        import RNA

        mfe_values = []

        for seq_record in SeqIO.parse(filename, "fasta"):
            seq_str = str(seq_record.seq)
            structure, mfe = RNA.fold(seq_str)
            mfe_values.append(mfe)

        fold_table = pd.DataFrame({'mfe': mfe_values})

        fold_table_z = (fold_table['mfe'] - fold_table['mfe'].mean()) / fold_table['mfe'].std()
        fold_table_wt_z = fold_table_z.iloc[0]
        fold_table_wt_z_rep = np.repeat(fold_table_wt_z, len(fold_table_z))
        fold_table_z_ls = (fold_table_wt_z_rep - fold_table_z)**2

        least_squares = least_squares.add(fold_table_z_ls, axis=0)

        final_table.insert(0, "VFOLD (mfe)", fold_table['mfe'])

        if args.graphics:
            fold_graphname = seq_name[1:-1] + '_' + args.random_type + '.fold.pdf'
            fold_graph = (
                ggplot(fold_table, aes('mfe'))
                + geom_density()
                + labs(x="MFE", y="Frequency")
                + geom_vline(xintercept=[fold_table['mfe'].iloc[0]], colour="red", linetype="dashed")
                + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (FOLD)')
            )
            ggsave(fold_graph, filename=fold_graphname)


    # UFOLD
    if 'UFOLD' in args.modules:
        # Needs implementation
        pass

    # DN
    if 'DN' in args.modules or 'all' in args.modules:
        print("Calculating DN")
        dnname = seq_name[1:-1]+'_'+args.random_type+'.dn'
        dn_file=open(dnname, 'w')
        dn_file.write("id")
        for nt1 in nts:
            for nt2 in nts:
                dinut = nt1+nt2
                dn_file.write("\t"+dinut)                
        dn_file.write("\n") 

        for nuc_rec in SeqIO.parse(filename, "fasta"):
            nucs = [str(nuc_rec.seq[i:i+1]) for i in range(len(nuc_rec.seq))]
            dinucs = [str(nuc_rec.seq[i:i+2]) for i in range(len(nuc_rec.seq)-1)]  
            nuc_counts = Counter(nucs)
            dinuc_counts = Counter(dinucs)
            seq_len = len(nuc_rec.seq)
            dn_file.write(nuc_rec.id)
            for nt1 in nts:
                for nt2 in nts:
                    dinut = nt1+nt2
                    if (dinut in dinuc_counts):
                        freq = (dinuc_counts[dinut] / (seq_len - 1)) / ((nuc_counts[nt1]/seq_len)*(nuc_counts[nt2]/seq_len))
                        dn_file.write("\t"+str(freq))
                    else:
                        dn_file.write("\t0")
            dn_file.write("\n")
        dn_file.close()

        dnlsname = seq_name[1:-1]+'_'+args.random_type+'.dnls'
        dn_table = pd.read_csv(dnname, sep='\t')

        dn_table_least = np.sqrt(((dn_table.iloc[:,1:] - dn_table.iloc[0,1:])**2).sum(axis=1))
        dn_table_least_z =  ((dn_table_least - dn_table_least.mean()) / dn_table_least.std())
        dn_table_least_z_wt = dn_table_least_z.iloc[0]
        dn_table_least_z_wt_rep = np.repeat(dn_table_least_z_wt, len(dn_table_least_z))
        dn_table_least_z_ls = (dn_table_least_z_wt_rep - dn_table_least_z)**2
        least_squares = least_squares.add(dn_table_least_z_ls, axis=0)
        dn_table_ls = pd.DataFrame({'Replication':range(len(dn_table_least)),'DN_least_square':dn_table_least})
        dn_table_ls.to_csv(dnlsname, sep="\t",index=False)
        final_table.insert(0, "DN_least_square", dn_table_ls['DN_least_square'])

        if args.graphics:
            dnls_graphname = dnlsname + '.pdf'
            dn_graph = (ggplot(dn_table_ls, aes('DN_least_square'))
                        + geom_density()
                        + labs(x="Dinucleotide", y="Dinucleotide Least Squares")
                        + geom_vline(xintercept=[dn_table_ls['DN_least_square'].iloc[0]], colour="red", linetype="dashed")
                        + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (DN)'))
            ggsave(dn_graph, filename=dnls_graphname)

    # CPB
    if 'CPB' in args.modules or 'all' in args.modules:
        print("Calculating CPB")
        cpbname = seq_name[1:-1]+'_'+args.random_type+'.cpb'
        cpb_file=open(cpbname,'w')
        cps_human = pd.read_csv("Coleman_CPS.csv", sep=';')
        cps_human = cps_human.drop(['Aapair','Expected','Observed','Observed/Expected'], axis=1)
        cps_human = cps_human.sort_values(['CodonPair'], ascending=[True])

        for nuc_rec in SeqIO.parse(filename, "fasta"):
            prot_rec = make_protein_record(nuc_rec)
            codons = [str(nuc_rec.seq[i:i+3]) for i in range(0,len(nuc_rec.seq)-3,3)]
            dicodons = [str(nuc_rec.seq[i:i+6]) for i in range(0,len(nuc_rec.seq)-3,3)]
            aas = [str(prot_rec.seq[i:i+1]) for i in range(len(prot_rec.seq))]
            diaas = [str(prot_rec.seq[i:i+2]) for i in range(len(prot_rec.seq)-1)]
            codon_counts = Counter(codons)
            dicodon_counts = Counter(dicodons)
            aa_counts = Counter(aas)
            diaa_counts = Counter(diaas)
            dicodon_df = pd.DataFrame.from_dict(dicodon_counts, orient='index').reset_index()
            dicodon_df = dicodon_df.sort_values(['index'], ascending=[True])
            dicodon_df.columns = ['CodonPair', 'Obs']
            cps_tb_final = pd.merge(cps_human, dicodon_df, on='CodonPair', how='inner')
            cps_tb_final['CPS'] = cps_tb_final['CPS'].replace({',':'.'}, regex=True).astype(float)
            cps_tb_final['result'] = cps_tb_final.CPS * cps_tb_final.Obs
            cpb = sum(cps_tb_final['result'])/sum(cps_tb_final['Obs'])
            cpb_file.write(str(cpb)+"\n")
        cpb_file.close()

        cpb_table = pd.read_csv(cpbname, sep=' ', names=['cpb'])
        cpb_table_z =  ((cpb_table['cpb'] - cpb_table['cpb'].mean()) / cpb_table['cpb'].std())
        cpb_table_wt_z = cpb_table_z.iloc[0]
        cpb_table_wt_z_rep = np.repeat(cpb_table_wt_z, len(cpb_table_z))
        cpb_table_z_ls = (cpb_table_wt_z_rep-cpb_table_z)**2
        least_squares = least_squares.add(cpb_table_z_ls, axis=0)
        final_table.insert(0, "CPB", cpb_table['cpb'])

        if args.graphics:
            cpb_graphname = cpbname +'.pdf'
            cpb_graph = (ggplot(cpb_table, aes('cpb'))
                         + geom_density()
                         + labs(x="CPB", y="Frequency")
                         + geom_vline(xintercept=[cpb_table['cpb'].iloc[0]], colour="red", linetype="dashed")
                         + ggtitle(seq_name[1:-1]+'_'+args.random_type+' (CPB)'))
            ggsave(cpb_graph, filename=cpb_graphname)

    print("Calculating Least Squares")
    least_squares = np.sqrt(least_squares)
    least_squares.columns = ['distance']
    least_table_name = seq_name[1:-1]+'_'+args.random_type+'_least_square.txt'
    least_squares.to_csv(least_table_name, sep="\t")
    final_table.insert(0, "Distance(ls)", least_squares['distance'])

    print("Making final table")
    nuc_distance_name=filename+'_distance_table.txt'
    nuc_distance_file=open(nuc_distance_name,'w')
    for j in range(1,n):
        difference = get_difference(seq_records[0], seq_records[j])
        nuc_distance_file.write(str(difference)+'\n')
    nuc_distance_file.close()
    col_name = ['Nucleotide_difference']
    new_nuc_distance_table = pd.read_csv(nuc_distance_name, sep=' ', names=col_name)
    new_nuc_distance_table.loc[-1]=[0]
    new_nuc_distance_table.index = new_nuc_distance_table.index + 1

    new_nuc_distance_table = new_nuc_distance_table.sort_index()

    final_table.insert(1, "Nucleotide_difference", new_nuc_distance_table['Nucleotide_difference'])

    final_tb_name = seq_name[1:-1]+'_'+args.random_type+'_final_table.txt'
    final_table.to_csv(final_tb_name, sep='\t')

    if args.graphics:
        new_table=pd.DataFrame()
        new_table.insert(0, "Distance", least_squares['distance'])
        new_table.insert(1, "Nucleotide_difference", new_nuc_distance_table['Nucleotide_difference'])
        final_graphname = filename +'final_graph.pdf'
        final_graph = (ggplot(new_table, aes('Distance', 'Nucleotide_difference'))
                       + geom_point()
                       + labs(x="Least Square Distance", y="Hamming Distance (nt)")
                       + ggtitle(seq_name[1:-1]+'_'+args.random_type))
        ggsave(final_graph, filename=final_graphname)
