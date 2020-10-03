#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import scanpy as sc
import numpy as np
import logomaker as lm
import matplotlib.pyplot as plt
import warnings
import re
from myLogo import showlogo

warnings.filterwarnings('ignore')


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  default="test.txt",
                                    help='The input should be a list of peptides with counts, '
                                    'specificity and frequency, by default test.txt is read. ')
    my_parser.add_argument('-wOption', nargs='?',  default="no",
                                    help='Do you want to set weight for each position?')   
    
    my_parser.add_argument('-b', nargs='?', default="BLOSUM62.txt",
                                    help='Blosum matrix, by default it is BLOSUM 62. ')
    my_parser.add_argument('-colorClass', nargs='*', default="groups",
                                    help='PCA UMAP color by groups, rounds or cellTypes. ')
    
    my_parser.add_argument('-tsne', nargs='?', default="0",
                                    help='Please type an integer if you want to plot tsne. '
                                    '1: TSNE will be generated. Warning: it can take up to 10 mins.'
                                    '0: No TSNE plot. By default no tsne plot.' )  
    
    my_parser.add_argument('-o', nargs='?', default="1",
                                    help='Please type an integer if you want to add weight the Bolsum matrix. '
                                    '1: no weight. 2: peptide frequency used as weight. '
                                    '3: peptide frequency used as an additional feature. '
                                    'By default it is option 1. ') 
    my_parser.add_argument('-w', nargs='?',
                                    help='Please insert weight string.')

    return my_parser.parse_args()



def read_blossum_matrix(blosum_matrix_file):
    
    AA = ""
    with open(blosum_matrix_file) as fileBlo:

        refBlo = {}
        index=20

        blosumAA = list("ARNDCQEGHILKMFPSTWYV")
        while True:
            line = fileBlo.readline().rstrip()
            if (len(line) == 0):
                break
            if ("#" not in line):
                arr = line.split()
                countAA = 0
                for i in range(0,len(arr)):
                    if (arr[i] in blosumAA):
                        countAA = countAA + 1
                if (countAA > 1):
                    for i in range(0,len(arr)):
                        if (arr[i]  in "BZX*"):
                            index = i
                            break
                        else:
                            AA = AA + arr[i]

                else:
                    if (arr[0] not in "BZX*"):
                        score1 = []
                        for i in range(1,index + 1):                    
                            score = arr[i]                    
                            score1.append(score)
                        refBlo[arr[0]] = " ".join(score1)
    AAs = list(AA)
    return (refBlo,AAs)

def run_analysis(fileTmp,wOption, refBlo, AAs,colorClasses,weights,blosum_option):  
                
    refScore = {}        
    refPep = {}
    
    with open(fileTmp, 'r') as file0:                
        for line0 in file0:            
            if (len(line0.rstrip()) == 0):
                break
        
            arr = line0.rstrip().split("\t")
            
            
            
            if (len(arr[0]) == 0):
                refPep["Peptides"]= arr[1:] 
                
                classes = ["groups", arr[2], arr[3]]
                refScore["Peps"] = classes
                
            else:
                classes =   [arr[1], arr[2], arr[3]]
                
                pep=arr[0]               
                
                freq = float(arr[-2])
                peps = list(pep)

                if (wOption == "no"):
                    weights = [1] * len(peps)

                if (len(peps) == len(weights)):
                    refPep[pep]= arr[1:] 
                    score=[]
                    for i in range(0,len(peps)):
                        aa = peps[i]
                        weight = float(weights[i])
                        score1 = refBlo.get(aa).split(" ")
                        for tmp in score1:
                            tmp2 = ""
                            if (blosum_option == "blosum" or blosum_option == "blosumFreqAddfreq"):
                                tmp2 = float(tmp) * weight             ##option1: no frequency ##option3: add freq as a feature ##add weight to each position
                            elif (blosum_option == "blosumMultiplyWeight"):
                                tmp2 = float(tmp) * weight * freq ##option2: add freq as weight into matrix                     ##add weight to each position
                            else:
                                print ("Blosum option not correct.")

                            score.append(str(tmp2)) 
                    if (blosum_option == "blosum" or blosum_option=="blosumMultiplyWeight"):
                        refScore[pep] = classes + score          ##option1: no frequency ##option2: add freq as weight into matrix
                    elif (blosum_option == "blosumFreqAddfreq"):
                        refScore[pep] = classes + score + [freq] ##option3: add freq as a feature
                    else:
                        print ("Blosum option is not correct.")

                else:
                    print ("Your weights: ", weights, ", length is " + str(len(weights)) + ". Length of peptides: " + str(len(peps)) )
                    sys.exit("The number of weights does not match the number of positions.")
    arrTitle = []        
    for tmp in AAs:
        for leng in range(1,len(weights) + 1):
            arrTitle.append("P" + tmp + str(leng))
   
    refScore["Peps"] = refScore.get("Peps") + arrTitle
    
    print ("Your weights on each AA position:", ", ".join (map(str,weights)))
    return (refScore, refPep)


def figurePlots(refScore, figDIR, tSNE, colorClasses):
   
    sc.settings.autosave = False
    sc.settings.figdir = figDIR
    
    df = pd.DataFrame.from_dict(refScore, orient="index")
    new_header = df.iloc[0] 
    df = df[1:] 
    df.columns = new_header
    
    ####dataframe should contain at least 5 peptides:
    if (df.shape[0] < 5):
        sys.exit("There should be at least 5 peptides for clustering analysis!")
    else:
          
        new_cols = list(set(df.columns) - set(["cells","groups","rounds"]))
        new_cols.sort()
        rawdf = df[new_cols]
        #rawdf.to_csv("test_xli", sep="\t", index = True, header=True)

        adata = sc.AnnData(rawdf) 
        
        myColors = ["#490092","#24FF24","#FF6DB6", "#009292", "#396AB1","#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#e1f5fe","#1A4756","#33918D","#F6D82D","#E39517","#D4342D"]
        
        adata.obs[colorClasses] = df[colorClasses]     
        
        sc.tl.pca(adata, svd_solver='arpack', n_comps=49)
        sc.pl.pca(adata,color=colorClasses,title=["PCA " + colorClass for colorClass in colorClasses], palette = myColors,components=['1,2'],hspace=0.2, ncols=1)
        plt.savefig(figDIR + "PCA_1vs2.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca(adata,color=colorClasses,title=["PCA " + colorClass for colorClass in colorClasses], palette = myColors,components=['1,3'],hspace=0.2, ncols=1)
        plt.savefig(figDIR + "PCA_1vs3.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca(adata,color=colorClasses,title=["PCA " + colorClass for colorClass in colorClasses], palette = myColors,components=['2,3'],hspace=0.2, ncols=1)
        plt.savefig(figDIR + "PCA_2vs3.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca_variance_ratio(adata,n_pcs=10)
        plt.savefig(figDIR + "PCA_variance_ratio.pdf",bbox_inches='tight')
        plt.close("all")

        if (tSNE == "1"):
            
            sc.tl.tsne(adata, n_pcs = 49)
            sc.pl.tsne(adata, color=colorClasses,title=["tSNE " + colorClass for colorClass in colorClasses],hspace=0.2, ncols=1)
            plt.savefig(figDIR + "tSNE.pdf",bbox_inches='tight')
            plt.close("all")

        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=49)
        sc.tl.umap(adata)        
        sc.tl.louvain(adata)
        sc.pl.umap(adata, color=colorClasses + ['louvain'],title=["UMAP " + colorClass for colorClass in colorClasses] + ['louvain'],palette = myColors,hspace=0.2, ncols=1)
        plt.savefig(figDIR + "UMAP_Louvain.pdf",bbox_inches='tight')
        plt.close("all")

        
        return adata

    #########Output the clusters:
def write_data_out(adata, figDIR, louvainClusterDIR, sumDIR, refPep, colorClasses):
    
    L = list(set(adata.obs['louvain'].values))
    L.sort()

    
    for i in range(0,len(colorClasses)):
        colorClass = colorClasses[i]
        
        outLou = open(figDIR + "sum.txt2_" + colorClass,"w")
        outLou.write("Clusters" + "\t" + "Group" + "\t" + "No.Peptides" + "\t" + "No.Total" + "\t" + "Percent" + "\t" + "Meta" + "\t" + "Frequency" + "\n" )
        
        refGroups = {}
        
        header = refPep.get("Peptides")
        if (colorClass == "groups"):
            index = 0
            for head in header:
                if ("_" in head):
                    refGroups[head.split("_")[0] + "_" + head.split("_")[1]] = 1
        elif (colorClass == "rounds"):
            index = 1
            for head in header:
                if ("_" in head):
                    refGroups[head.split("_")[0]] = 1
        else:
            index = 2
            for head in header:
                if ("_" in head):
                    refGroups[head.split("_")[1]] = 1

        arrGroups = list(refGroups)
        
        for l in L:
            ##write the peptides 
            pepout = open(louvainClusterDIR + "Cluster" +  l + ".txt", 'w')
            pepout.write("Cluster" + l + "_Peptides"+ "\t" + "\t".join(refPep.get("Peptides")) + "\n")

            ldata = adata[adata.obs['louvain'] == l] 
            
            peps = list(ldata.obs['louvain'].index)
            for pep in peps:
                pepout.write(pep + "\t" + "\t".join(refPep.get(pep)) + "\n")
            pepout.close()

            nrows = 1
            height = 2.5

            fig = plt.figure(figsize=(10,height * nrows))
            logomat  = lm.alignment_to_matrix(sequences=peps)        
            logomatI = lm.transform_matrix(df = logomat, from_type='counts', to_type='information')
            newTitle = "Cluster " +  l
            
            axis = fig.add_subplot(1,1,1)
            logo = showlogo(logomatI, newTitle, "logo",axis)
            plt.savefig(louvainClusterDIR + "Cluster" + l + ".pdf")
            plt.close('all')
            
            myGroups = ldata.obs[colorClass].unique()
            total=len(ldata)
            refGroups = {}
            for myGroup in myGroups:
                myGroupTmps = myGroup.split("+")
                for myGroupTmp in myGroupTmps:
                    refGroups[myGroupTmp] = 1


            myAllGroups = list(refGroups)
            
            for myGroupTmp in myAllGroups:
                refPepsTmp = {}
                for pep in peps:
                    if (myGroupTmp in refPep.get(pep)[index].split("+")):
                        refPepsTmp[pep] = 1
                
                myPercent = "{:0.2f}".format((len(refPepsTmp)/(total + 0.0))*100)
                outLou.write("Cluster" + l + "\t" + myGroupTmp + "\t" + str(len(refPepsTmp)) + "\t" + str(total) + "\t" + str(myPercent) + "\n"  )

        outLou.close()


        df = pd.read_csv(figDIR + "sum.txt2_" + colorClass,index_col=None, header=0,sep="\t",dtype=str)  

        df["Clusters"] = df["Clusters"].str.replace("Cluster","")
        df["Clusters"] = df["Clusters"].astype(int)
        df["No.Peptides"] = df["No.Peptides"].astype(int)
        df["Percent"] = df["Percent"].astype(float)
        df["No.Total"] = df["No.Total"].astype(int)


        df = pd.pivot_table(df,index=["Clusters","No.Total"], columns=["Group"], values=["No.Peptides","Percent"], aggfunc='first')
        new_header = ['%s%s' % (b, '_%s' % a if b else '') for a, b in df.columns]
        df.columns = new_header
        groups_suffix =["No.Peptides","Percent"]
        
        
        new_comb = []
        for arrGroup in arrGroups:
            for group_suffix in groups_suffix:
                comb = arrGroup + "_" + group_suffix
                if (comb in list(df.columns)):
                    new_comb.append(comb)
        
        df = df[new_comb]
               
        df.fillna("NA").to_csv(sumDIR + "sum_stat_" + colorClass + ".txt", sep="\t", index = True, header=True)
        #df2=pd.DataFrame(df.to_records())
        #df2.fillna("NA").to_json("sum_" + colorClass + ".json")
        #######END######    


def main():
    
    args = parse_arguments()
    refArgs = vars(args)
    
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") 
    
    colorClasses = refArgs.get("colorClass") 
    blosum_matrix = refArgs.get("b") ## "BLOSUM62.txt", "blosum62_freq.txt"
    
    tSNE = refArgs.get("tsne") ## whether or not to plot tsne
    wOption = refArgs.get("wOption")## Do you want to set weights for each position?
    weight = ""
    if (wOption == "yes"):
        weight = refArgs.get("w") 
   
    
    delimiters = ["X",",", ";", " "]
    regexPattern = '|'.join(map(re.escape, delimiters))
    weightsTmp= re.split(regexPattern, weight)
    weights = [x for x in weightsTmp if x] #weight for each position
    

    blosum_option=""
    blosum_option_message=""
    if (refArgs.get("o")=="1"):
        blosum_option = "blosum"
        blosum_option_message = "No weight, original Blosum matrix"
    elif (refArgs.get("o")=="2"):
        blosum_option = "blosumMultiplyWeight" 
        blosum_option_message = "Peptide frequency used as weight"
    elif (refArgs.get("o")=="3"):
        blosum_option = "blosumFreqAddfreq" 
        blosum_option_message = "Peptide frequency used as an additional feature"
    else:
        sys.exit("option is invalid. It should be 1,2 or 3.")
    print ("Peptides will be colored by: ", colorClasses) 
    print ("Your Blosum matrix: ", blosum_matrix.split("/")[-1].replace(".txt",""))
    print ("Peptide matrix option: ", blosum_option_message)
    
    refBlo, AAs = read_blossum_matrix(blosum_matrix)
    
    figDIR =  "figures_PCA_UMap_Louvain/"
    louvainClusterDIR =  "Louvain_Clusters/"
    sumDIR = "sum_stats/"
    
    os.system("mkdir " + figDIR)
    os.system("mkdir " + louvainClusterDIR )
    os.system("mkdir " + sumDIR )


    refScore,refPep = run_analysis(myFile, wOption, refBlo, AAs,colorClasses,weights,blosum_option)
    
    adata = figurePlots(refScore, figDIR,tSNE,colorClasses)
    write_data_out(adata, figDIR, louvainClusterDIR,sumDIR,refPep,colorClasses)

if __name__ == "__main__":
    main()






