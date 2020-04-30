#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 19:48:19 2018

@author: amina
"""


from Bio import SeqIO
import numpy as np

from Bio.SeqIO import FastaIO



import os
import random
import sklearn

from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from scipy.stats import rankdata
import math
import pickle
from scipy.stats import rankdata
#import time
from sklearn.externals import joblib
from sklearn.preprocessing import normalize
def prot_feats_seq(seq):

    aa=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


    f=[]



    X = ProteinAnalysis(str(seq))

    X.molecular_weight() #throws an error if 'X' in sequence. we skip such sequences
    p=X.get_amino_acids_percent()

    dp=[]
    for a in aa:
        dp.append(p[a])
    dp=np.array(dp)
    dp=normalize(np.atleast_2d(dp), norm='l2', copy=True, axis=1, return_norm=False)
    f.extend(dp[0])
    tm=np.array(twomerFromSeq(str(seq)))
    tm=normalize(np.atleast_2d(tm), norm='l2', copy=True, axis=1,return_norm=False)

    f.extend(tm[0])
    thm=np.array(threemerFromSeq(str(seq)))
    thm=normalize(np.atleast_2d(thm), norm='l2', copy=True, axis=1,return_norm=False)
    f.extend(thm[0])


    return np.array(f)



def prot_feats(filename):
    XX=[]
    ids=[]


    for rec in SeqIO.parse(filename, "fasta"):
        f=[]
        X = ProteinAnalysis(str(rec.seq))
#        import pdb; pdb.set_trace()
        try:
            X.molecular_weight() #throws an error if 'X' in sequence. we skip such sequences
            f=list(prot_feats_seq(str(rec.seq)))
    #
            XX.append(f)
            ids.append(rec.id)
            
        except:
#            print ("exception")
            continue






    XX=np.array(XX)
#    import pdb; pdb.set_trace()

    return XX,ids






def twomerFromSeq(s):
    k=2
    groups={'A':'1','V':'1','G':'1','I':'2','L':'2','F':'2','P':'2','Y':'3',
            'M':'3','T':'3','S':'3','H':'4','N':'4','Q':'4','W':'4',
            'R':'5','K':'5','D':'6','E':'6','C':'7'}
    crossproduct=[''.join (i) for i in product("1234567",repeat=k)]
    for i in range (0,len(crossproduct)): crossproduct[i]=int(crossproduct[i])
    ind=[]
    for i in range (0,len(crossproduct)): ind.append(i)
    combinations=dict(zip(crossproduct,ind))

    V=np.zeros(int((math.pow(7,k))))      #defines a vector of 343 length with zero entries
    try:
        for j in range (0,len(s)-k+1):
            kmer=s[j:j+k]
            c=''
            for l in range(0,k):
                c+=groups[kmer[l]]
                V[combinations[int(c)]]+=1
    except:
        count={'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
        for q in range(0,len(s)):
            if s[q]=='A' or s[q]=='V' or s[q]=='G':
                count['1']+=1
            if s[q]=='I' or s[q]=='L'or s[q]=='F' or s[q]=='P':
                count['2']+=1
            if s[q]=='Y' or s[q]=='M'or s[q]=='T' or s[q]=='S':
                count['3']+=1
            if s[q]=='H' or s[q]=='N'or s[q]=='Q' or s[q]=='W':
                count['4']+=1
            if s[q]=='R' or s[q]=='K':
                count['5']+=1
            if s[q]=='D' or s[q]=='E':
                count['6']+=1
            if s[q]=='C':
                count['7']+=1
        val=list(count.values())              #[ 0,0,0,0,0,0,0]
        key=list(count.keys())                #['1', '2', '3', '4', '5', '6', '7']
        m=0
        ind=0
        for t in range(0,len(val)):     #find maximum value from val
            if m<val[t]:
                m=val[t]
                ind=t
        m=key [ind]                     # m=group number of maximum occuring group alphabets in protein
        for j in range (0,len(s)-k+1):
            kmer=s[j:j+k]
            c=''
            for l in range(0,k):
                if kmer[l] not in groups:
                    c+=m
                else:
                    c+=groups[kmer[l]]
            V[combinations[int(c)]]+=1

    V=V/(len(s)-1)
    return np.array(V)


def twomerFromFile(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    X=[]
    for rec in records:
        s=str(rec.seq)
        X.append(twomerFromSeq(s))
    return X

def threemerFromSeq(s):
    k=3
    groups={'A':'1','V':'1','G':'1','I':'2','L':'2','F':'2','P':'2','Y':'3',
            'M':'3','T':'3','S':'3','H':'4','N':'4','Q':'4','W':'4',
            'R':'5','K':'5','D':'6','E':'6','C':'7'}
    crossproduct=[''.join (i) for i in product("1234567",repeat=k)]
    for i in range (0,len(crossproduct)): crossproduct[i]=int(crossproduct[i])
    ind=[]
    for i in range (0,len(crossproduct)): ind.append(i)
    combinations=dict(zip(crossproduct,ind))

    V=np.zeros(int((math.pow(7,k))))      #defines a vector of 343 length with zero entries
    try:
        for j in range (0,len(s)-k+1):
            kmer=s[j:j+k]
            c=''
            for l in range(0,k):
                c+=groups[kmer[l]]
                V[combinations[int(c)]]+=1
    except:
        count={'1':0,'2':0,'3':0,'4':0,'5':0,'6':0,'7':0}
        for q in range(0,len(s)):
            if s[q]=='A' or s[q]=='V' or s[q]=='G':
                count['1']+=1
            if s[q]=='I' or s[q]=='L'or s[q]=='F' or s[q]=='P':
                count['2']+=1
            if s[q]=='Y' or s[q]=='M'or s[q]=='T' or s[q]=='S':
                count['3']+=1
            if s[q]=='H' or s[q]=='N'or s[q]=='Q' or s[q]=='W':
                count['4']+=1
            if s[q]=='R' or s[q]=='K':
                count['5']+=1
            if s[q]=='D' or s[q]=='E':
                count['6']+=1
            if s[q]=='C':
                count['7']+=1
        val=list(count.values())              #[ 0,0,0,0,0,0,0]
        key=list(count.keys())                #['1', '2', '3', '4', '5', '6', '7']
        m=0
        ind=0
        for t in range(0,len(val)):     #find maximum value from val
            if m<val[t]:
                m=val[t]
                ind=t
        m=key [ind]                     # m=group number of maximum occuring group alphabets in protein
        for j in range (0,len(s)-k+1):
            kmer=s[j:j+k]
            c=''
            for l in range(0,k):
                if kmer[l] not in groups:
                    c+=m
                else:
                    c+=groups[kmer[l]]
            V[combinations[int(c)]]+=1

    V=V/(len(s)-1)
    return np.array(V)



def threemerFromFile(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    X=[]
    for rec in records:
        s=str(rec.seq)
        X.append(threemerFromSeq(s))
    return X



def analyzePhage(filename, output):
    X,ids=prot_feats(filename)

    model=pickle.load(open('xgb_rank.pickle', 'rb'))
    #model=joblib.load('/home/acranker/mysite/xgb_rank.pkl')

    scores =  model.predict(X)

    ranks=len(scores) - rankdata(scores, method='ordinal')+1

    ranks, scores, ids= zip(*sorted(zip(ranks, scores, ids)))
    res_fn=output+'.csv'
    with open(res_fn,'w') as f:
        f.write('protein ID,Rank,Score\n')

        for i in range(len(ids)):
            #f.write(str(X[i]))
            f.write(str(ids[i])+','+str(ranks[i])+','+str(scores[i])+'\n')
    f.close()
