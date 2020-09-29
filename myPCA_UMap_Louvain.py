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
    my_parser.add_argument('-colorClass', nargs='?', default="groups",
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


def showlogo(df,myTitle):
    mylogo=lm.Logo(df,font_name="DejaVu Sans")
    mylogo.style_spines(visible=False)
    mylogo.style_spines(spines=["left","bottom"],visible=True)
    mylogo.ax.set_ylabel("Information (bits)")
    mylogo.ax.set_title(myTitle)
    mylogo.ax.xaxis.set_ticks_position("none")
    mylogo.ax.xaxis.set_tick_params(pad=-1)
    mylogo.ax.grid(False)
    mylogo.fig.tight_layout()


def figurePlots(fileTmp,wOption, refBlo, AAs,colorClass,weights,blosum_option):  
                
    refScore = {}    
    refScore["Peps"] = ["Classes"] 
    refPep = {}
    sorted_groups = [] 
    
    with open(fileTmp, 'r') as file0:                
        for line0 in file0:            
            if (len(line0.rstrip()) == 0):
                break
        
            arr = line0.rstrip().split("\t")
            
            if (len(arr[0]) == 0):
            
                for i in range (1,len(arr)):
                    if ("_" in arr[i]):
                        group_tmp = ""
                        if (colorClass == "groups"):
                            group_tmp=arr[i].split("_")[0] + "_" + arr[i].split("_")[1]
                        elif (colorClass == "rounds"):
                            group_tmp = arr[i].split("_")[0]
                        else:
                            group_tmp = arr[i].split("_")[1]
                        if (group_tmp not in sorted_groups):
                            sorted_groups.append(group_tmp)

                refPep["Peptides"]= arr[1:] 

            else:
                
                pep=arr[0]

                if (colorClass == "groups"):
                    classes = arr[1]
                elif (colorClass == "rounds"):
                    classes = arr[2]
                else:
                    classes = arr[3]

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
                        refScore[pep] = [classes] + score          ##option1: no frequency ##option2: add freq as weight into matrix
                    elif (blosum_option == "blosumFreqAddfreq"):
                        refScore[pep] = [classes] + score + [freq] ##option3: add freq as a feature
                    else:
                        print ("Blosum option is not correct.")

                else:
                    print ("Your weights: ", weights, ", length is " + str(len(weights)) + ". Length of peptides: " + str(len(peps)) )
                    sys.exit("The number of weights does not match the number of positions.")
    arrTitle = []        
    for tmp in AAs:
        for leng in range(1,len(weights) + 1):
            arrTitle.append("P" + tmp + str(leng))
   
    refScore["Peps"] = ["Classes"] + arrTitle
    print ("Your weights on each AA position:", ", ".join (weights))
    return (refScore, refPep, sorted_groups)


def run_analysis(refScore, figDIR, tSNE):
   
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
        #######START######    
        rawdf=df.loc[:, df.columns != 'Classes']
        #rawdf.to_csv("test_xli", sep="\t", index = True, header=True)

        adata = sc.AnnData(rawdf)    
        adata.obs['Classes'] = df['Classes']

        myColors = ["#490092","#24FF24","#FF6DB6", "#009292", "#396AB1","#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#e1f5fe","#1A4756","#33918D","#F6D82D","#E39517","#D4342D"]
  
        sc.tl.pca(adata, svd_solver='arpack', n_comps=49)
        sc.pl.pca(adata,color=['Classes'],title="PCA",palette = myColors,components=['1,2'])
        plt.savefig(figDIR + "PCA_1vs2.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca(adata,color=['Classes'],title="PCA",palette = myColors,components=['1,3'])
        plt.savefig(figDIR + "PCA_1vs3.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca(adata,color=['Classes'],title="PCA",palette = myColors,components=['2,3'])
        plt.savefig(figDIR + "PCA_2vs3.pdf",bbox_inches='tight')
        plt.close("all")
        sc.pl.pca_variance_ratio(adata,n_pcs=10)
        plt.savefig(figDIR + "PCA_variance_ratio.pdf",bbox_inches='tight')
        plt.close("all")

        if (tSNE == "1"):
            
            sc.tl.tsne(adata, n_pcs = 49)
            sc.pl.tsne(adata, color='Classes',title="tSNE")
            plt.savefig(figDIR + "tSNE.pdf",bbox_inches='tight')
            plt.close("all")

        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=49)
        sc.tl.umap(adata)        
        sc.pl.umap(adata, color=['Classes'],title="UMAP",palette = myColors)
        plt.savefig(figDIR + "UMAP.pdf",bbox_inches='tight')
        plt.close("all")

        sc.tl.louvain(adata)
        sc.pl.umap(adata, color=['louvain'],title="Louvain")
        plt.savefig(figDIR + "Louvain.pdf",bbox_inches='tight')
        plt.close("all")
        
        return adata

    #########Output the clusters:
def write_data_out(adata, figDIR, louvainClusterDIR, refPep, colorClass, sorted_groups):
        
    with open(figDIR + "sum.txt2","w") as outLou:
        
        outLou.write("Clusters" + "\t" + "Group" + "\t" + "No.Peptides" + "\t" + "No.Total" + "\t" + "Percent" + "\t" + "Meta" + "\t" + "Frequency" + "\n" )
        L = list(set(adata.obs['louvain'].values))
        L.sort()

        for l in L:
            ldata = adata[adata.obs['louvain'] == l]      
            myGroups = ldata.obs['Classes'].unique()
            total=len(ldata)

            with open(louvainClusterDIR + "Cluster" +  l + ".txt", 'w') as pepout:
                pepout.write("Cluster" + l + "_Peptides"+ "\t" + "\t".join(refPep.get("Peptides")) + "\n")
                peps = list(adata[adata.obs['louvain'] == str(l)].obs['louvain'].index) # list of peptides
                for pep in peps:
                    pepout.write(pep + "\t" + "\t".join(refPep.get(pep)) + "\n")
            

            refCountList = {}
            refPepList = {}
            refRankList = {}
            refFreqList= {}
            for myGroup in myGroups:
                
                myCount = len (ldata[ldata.obs['Classes'] == myGroup])      
                myPercent = "{:0.2f}".format((myCount/(total + 0.0))*100)
                myGroupTmps = myGroup.split("+") ##To handle the groups with "+"
                 
                for myGroupTmp in myGroupTmps:
                    ######check count
                    myCount = len (ldata[ldata.obs['Classes'] == myGroup])           
                    myPercent = "{:0.2f}".format((myCount/(total + 0.0))*100)
                    myPeps = list( ldata[(ldata.obs['Classes'] == myGroup )].obs['louvain'].index)
                    ######check count
                    index_rank = 1000
                    index_freq = 1000
                    index_count= 1000
                    header = refPep.get("Peptides")
                       
                    for i in range(0,len(header)):
                        if ("_" in header[i] ):                       
                            suffix = header[i].split("_")[-1]
                            rou_exp = ""
                            if (colorClass == "groups"):
                                rou_exp = header[i].split("_")[0] + "_" + header[i].split("_")[1]
                            elif (colorClass == "rounds"):
                                rou_exp = header[i].split("_")[0]
                            else:
                                rou_exp = header[i].split("_")[1]
                            

                            if (rou_exp==myGroupTmp):
                                if (suffix == "Rank"):
                                    index_rank = i
                                if (suffix == "Frequency"):
                                    index_freq = i
                                if (suffix == "TotalCount"):
                                    index_count = i
                    arr_rank = []
                    arr_freq = []
                    arr_count= []

                    for myPep in myPeps:
                        arr_rank.append(refPep.get(myPep)[index_rank])
                        arr_freq.append(refPep.get(myPep)[index_freq])
                        arr_count.append(refPep.get(myPep)[index_count])

                    if (myGroupTmp not in refPepList):
                        refPepList[myGroupTmp] = myPeps
                    else:
                        refPepList[myGroupTmp] = refPepList.get(myGroupTmp) + myPeps

                    if (myGroupTmp not in refRankList):
                        refRankList[myGroupTmp] = arr_rank
                    else:
                        refRankList[myGroupTmp] = refRankList.get(myGroupTmp) + arr_rank
                    if (myGroupTmp not in refCountList):
                        refCountList[myGroupTmp] = arr_count
                    else:
                        refCountList[myGroupTmp] = refCountList.get(myGroupTmp) + arr_count
                    if (myGroupTmp not in refFreqList):
                        refFreqList[myGroupTmp] = arr_freq
                    else:
                        refFreqList[myGroupTmp] = refFreqList.get(myGroupTmp) + arr_freq
                    
                    
            for keys in refPepList:
                myPercent = "{:0.2f}".format((len(refPepList.get(keys))/(total + 0.0))*100)
                outLou.write("Cluster" + l + "\t" + keys + "\t" + str(len(refPepList.get(keys))) + "\t" + str(total) + "\t" + str(myPercent) + "\t" + "[[" + ";".join(refPepList.get(keys)) + "]" + ",[" + ";".join(refRankList.get(keys))  + "]" + ",[" + ";".join(refCountList.get(keys)) + "]]" + "\t" + "[" + ";".join(refFreqList.get(keys)) + "]" + "\n")

            fig = plt.figure()
            logomat  = lm.alignment_to_matrix(sequences=peps)        
            logomatI = lm.transform_matrix(df = logomat, from_type='counts', to_type='information')
            newTitle = "Cluster " +  l
            logo = showlogo(logomatI, newTitle)

            plt.savefig(louvainClusterDIR + "Cluster" + l + ".pdf")
            plt.close('all')


    df = pd.read_csv(figDIR + "sum.txt2",index_col=None, header=0,sep="\t",dtype=str)  
        
    df["Clusters"] = df["Clusters"].str.replace("Cluster","")
    df["Clusters"] = df["Clusters"].astype(int)
    df["No.Peptides"] = df["No.Peptides"].astype(int)
    df["Percent"] = df["Percent"].astype(float)
    df["No.Total"] = df["No.Total"].astype(int)
        
        
    df = pd.pivot_table(df,index=["Clusters","No.Total"], columns=["Group"], values=["Meta","Frequency","No.Peptides","Percent"], aggfunc='first')
    new_header = ['%s%s' % (b, '_%s' % a if b else '') for a, b in df.columns]
    df.columns = new_header
    groups_prefix=["No.Peptides","Percent","Meta","Frequency"]
    new_comb=[]
    for group in sorted_groups:
        for group_prefix in groups_prefix:
            new_comb_tmp = group + "_" + group_prefix 
            if (new_comb_tmp in df.columns):
                new_comb.append(new_comb_tmp)

    df = df[new_comb]
    all_columns = list(df.columns)
    sim_columns = []
    for col in all_columns:
        if ("_Frequency" not in col and "_Meta" not in col):
            sim_columns.append(col)

    df[sim_columns].fillna("NA").to_csv("sum.txt", sep="\t", index = True, header=True)
    df2=pd.DataFrame(df.to_records())
    df2.fillna("NA").to_json("sum.json")
        #######END######    


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

def main():
    
    args = parse_arguments()
    refArgs = vars(args)
    
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") ##peptide file with frequency
    
    colorClass = refArgs.get("colorClass") ##PCA UMAP color by classes, rounds or cellTypes
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
    print ("Peptides will be colored by: ", colorClass) 
    print ("Your Blosum matrix: ", blosum_matrix.split("/")[-1].replace(".txt",""))
    print ("Peptide matrix option: ", blosum_option_message)
    
    refBlo, AAs = read_blossum_matrix(blosum_matrix)
    
    figDIR =  "figures_PCA_UMap_Louvain/"
    louvainClusterDIR =  "Louvain_Clusters/"

    os.system("mkdir " + figDIR)
    os.system("mkdir " + louvainClusterDIR )

    refScore,refPep,sorted_groups = figurePlots(myFile, wOption, refBlo, AAs,colorClass,weights,blosum_option)
    
    adata = run_analysis(refScore, figDIR,tSNE)
    write_data_out(adata, figDIR, louvainClusterDIR,refPep,colorClass,sorted_groups)

if __name__ == "__main__":
    main()





