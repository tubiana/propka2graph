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
#==============================================================================
#                     GLOBAL VARIABLES                            
#==============================================================================
         

#==============================================================================
# TOOL FONCTION
#==============================================================================

def FileToDict(filePath): 
    pkaDict={}
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
            
            if pkaDict.has_key(resname):
                pkaDict[resname][1][0].append(number)
                pkaDict[resname][1][1].append(chain)
                pkaDict[resname][1][2].append(pka)
            else:
                pkaDict[resname]=[pkaModel,[[number],[chain],[pka]]]
        elif "     RESIDUE    pKa   pKmodel   ligand atom-type" in line:
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
    
    args = vars(arguments.parse_args())
    return(args)
    

def DictToGraph(pkaDict, filename):
    resi=["ASP","GLU", "HIS","TYR","LYS"]
    for res in resi:
        fig, ax = plt.subplots()
        pkaTheo=pkaDict[res][0]
        pkaCalc=pkaDict[res][1][2]
        pkaGraph=[x-pkaTheo for x in pkaCalc]
        x=range(len(pkaGraph))
        colors=[]
        for c in pkaGraph:
            if c>=0:
                colors.append("red")
            else:
                colors.append("blue")
        ax.bar(x,pkaGraph, color=colors, alpha=0.5)
        ax.margins(0.05)
        ax.set_xlabel("Res numbers")
        ax.set_ylabel("pKa normalized")
        labels=[]
        for i in xrange(len(pkaDict[res][1][0])):
            labels.append("%s-%s" %(pkaDict[res][1][0][i],pkaDict[res][1][1][i]))
        plt.xticks(x, labels, rotation=80)
        ax.tick_params("x",labelsize=6)    
        plt.title("%s -- pKa Theo : %.2f" %(res, pkaTheo))
        fig.savefig('%s-%s.png' %(filename,res),dpi=500, transparent=True)

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
    
    for filePath in files:
        filename=filePath.split("/")[-1].split(".propka")[0]
        pkaDict=FileToDict(filePath)
        DictToGraph(pkaDict, filename)

















