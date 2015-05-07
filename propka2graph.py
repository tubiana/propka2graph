#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION

ARGUMENTS : 

USAGE : 

"""

__author__ = "Thibault TUBIANA"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "20../.."
#==============================================================================
#                     MODULES                            
#==============================================================================
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
from operator import sub

#==============================================================================
#                     GLOBAL VARIABLES                            
#==============================================================================
PKAMODEL={}

#==============================================================================
# TOOL FONCTION
#==============================================================================

def multi_dimensions(n, type):
  """
  Creates an n-dimension dictionary where the n-th dimension is of type 'type'
  """  
  if n<=1:
    return type()
  return defaultdict(lambda:multi_dimensions(n-1, type))


def FileToDict(filePath): 
    pkaDict=multi_dimensions(3,dict)
    myFile=open(filePath,'r')
    find=False
    for line in myFile:
        if find:
            if "---------------------------" in line:
                break
            resname=line[:6].split()[0]
            number=int(line[6:10])
            chain=line[11]
            pka=float(line[15:20])
            pkaModel=float(line[25:30])
            if not PKAMODEL.has_key(resname):
                PKAMODEL[resname]=pkaModel
            pkaDict[chain][resname][number]=pka
        elif "     RESIDUE    pKa   pKmodel   ligand atom-type" in line or \
             "       Group      pKa  model-pKa   ligand atom-type" in line:
            find=True
    myFile.close()
    return pkaDict
            

        
        
#==============================================================================
# FONCTIONS
#==============================================================================


def parseArg():
    """
    This fonction will the list of pdb files and the distance
    @return: dictionnary of arguments
    """
    arguments=argparse.ArgumentParser(description="\
            Description\
            ")
    arguments.add_argument('-f', "--files", help="propka file(s)",
                           nargs="*", required=True)
    arguments.add_argument('-c',"--chains", help="chain to select (separated\
                            by a coma)", default=None)
    arguments.add_argument('-p',"--comparaison", help="Compare 2 propka file (Y/N)",
                           default="N")
    arguments.add_argument('-b',"--basic", help="basic graph?",
                           default="Y")
    args = vars(arguments.parse_args())
    return(args)
    

def DictToGraph(pkaDict, filename, chains_selected):
    resi=["ASP","GLU", "HIS","TYR","LYS", "ARG"]
    if not chains_selected:
        chains=pkaDict.keys()
        chains.sort()
    else:    
        chains=chains_selected.split(",")
        
        
    for res in resi:
        labels=[]
        fig, ax = plt.subplots()
        pkaTheo=PKAMODEL[res]
        numbers=[]
        pkaCalc=[]
        for chain in chains:
            numbers_temp=pkaDict[chain][res].keys()
            numbers_temp.sort()
            numbers.extend(numbers_temp)
            pkaCalc_temp=[pkaDict[chain][res][x] for x in numbers_temp]
            pkaCalc.extend(pkaCalc_temp)
            labels.extend(["%i-%s" %(n,chain) for n in numbers_temp])
            
        pkaGraph=[x-pkaTheo for x in pkaCalc]    
        colors=[]
        for c in pkaGraph:
            if c>=0:
                colors.append("lightcoral")
            else:
                colors.append("lightskyblue")
        
        
        x=range(len(pkaGraph))
        ax.bar(x,pkaGraph, color=colors, alpha=0.5,edgecolor="black")
        ax.margins(0.05)
        ax.set_xlabel("Res numbers")
        ax.set_ylabel("pKa normalized")
        plt.xticks(x, labels, rotation=80)
        ax.tick_params("x",labelsize=6)    
        plt.title("%s -- pKa Theo : %.2f" %(res, pkaTheo))
        fig.savefig('%s-%s.png' %(filename,res),dpi=500, transparent=True)
            

def compare_graph(pkaDictList, chains_selected, filenames):
    pkaDict0=pkaDictList[0]
    pkaDict1=pkaDictList[1]
    
    resi=["ASP","GLU", "HIS","TYR","LYS", "ARG"]
    if not chains_selected:
        chains=pkaDict.keys()
        chains.sort()
        print chains
    else:    
        chains=chains_selected.split(",")
        
    for res in resi:
        labels=[]
        fig, ax = plt.subplots()
        pkaTheo=PKAMODEL[res]
        numbers=[]
        pkaGraph=[]
        for chain in chains:
            num0=pkaDict0[chain][res].keys()
            num1=pkaDict1[chain][res].keys()
            numbers_temp=list(set(num0).intersection(num1))
            numbers_temp.sort()
            numbers.extend(numbers_temp)
            pka0=[pkaDict0[chain][res][x] for x in numbers_temp]
            pka1=[pkaDict1[chain][res][x] for x in numbers_temp]
            pka=map(sub,pka0,pka1)
            pkaGraph.extend(pka)
            labels.extend(["%i-%s" %(n,chain) for n in numbers_temp])

        colors=[]
        for c in pkaGraph:
            if c>=0:
                colors.append("lightcoral")
            else:
                colors.append("lightskyblue")
                        
        x=range(len(pkaGraph))
        ax.bar(x,pkaGraph, color=colors, alpha=0.5,edgecolor="black")
        ax.margins(0.05)
        ax.set_xlabel("Res numbers")
        ax.set_ylabel("Diff PKA")
        plt.xticks(x, labels, rotation=80)
        ax.tick_params("x",labelsize=6)    
        plt.title("Comparaison PKA %s et %s Res = %s-- pKa Theo : %.2f" %(filenames[0],filenames[1],res, pkaTheo))
        fig.savefig('Comp-%s-%s_%s.png' %(filenames[0],filenames[1],res),dpi=500, transparent=True)
        
    
    

###############################################################################
#####                               MAIN                                 ######
###############################################################################

if __name__ == "__main__":
    print "***********************************************"
    print "**********  PROPKA 2 GRAPH V1.0  **************"
    print "***********************************************"
    print ""
    #We get all arguments
    myArgs=parseArg()
    
    files=myArgs["files"]
    chains=myArgs["chains"]    
    compare=myArgs["comparaison"].upper()
    basic_graph=myArgs["basic"].upper()
    pkaDictList=[]
    filenames=[]
    for filePath in files:
        filename=filePath.split("/")[-1].split(".")[0]
        filenames.append(filename)
        pkaDict=FileToDict(filePath)
        if basic_graph=="Y":
            DictToGraph(pkaDict, filename, chains)
        pkaDictList.append(pkaDict)
       
    if compare=="Y":
        if len(pkaDictList)==2:
            compare_graph(pkaDictList, chains, filenames)
        else:
            "please give only 2 proka files if you want a comparaison"
        
           
















