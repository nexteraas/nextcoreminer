#!/bin/python
import argparse
import os
import pandas as pd

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='?',help='The input should be a list of peptides')
    my_parser.add_argument('-g',  nargs='?', default="1-15",help='Gibbs cluster')
    
    return my_parser.parse_args()
    
def main():

    args = parse_arguments()
    refArgs = vars(args)

    myFile = refArgs.get("f") 

    gibbs = "/home/galaxy/galaxy/tools/nextera/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl"
    gibs_cluster = refArgs.get("g").strip()

    cmd0="cut -f 1 " + myFile + " |awk '$0 ~ /^[A-Z]+$/' > myTmp0"
    cmd1=gibbs + " -f myTmp0 -g " + gibs_cluster +"  -P myGibbs -T -j 2 -S 5 -k 5 > myGibbs.log.txt"
    os.system(cmd0)
    os.system(cmd1)

    cmd2 = '''ls myGibbs_100/cores/ |grep -v "cores" > myTmp1'''
    cmd3 = "mkdir myGibbs_100/cores_real/"
    cmd4 = "mkdir myGibbs_100/res_real/"
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    
    df1 =  pd.read_csv("myTmp0",header=None,index_col=None)
    peps_all = list(df1.iloc[:,0])

    file0=open("myTmp1")
    while True:
        line0=file0.readline().rstrip()
        if (len(line0)==0):
            break
        txt = line0.split(".")[1]
        cmd5 = '''sed 1i\gibbs.''' + txt + "_Cores myGibbs_100/cores/" + line0 + " > myGibbs_100/cores_real/" + line0 + ".txt"
        os.system(cmd5)
        df2 =  pd.read_csv("myGibbs_100/cores_real/" + line0 + ".txt",header=0,index_col=None)  
        myCores = list(df2.iloc[:,0])
        
        refPep = {}
        for myCore in myCores:
            for s in filter(lambda x: myCore in x, peps_all) :
                refPep[s] = 1

        out = open("myGibbs_100/res_real/" + line0.replace("core","pep") + ".txt","w")
        out.write("gibbs." + txt + "_Peptides" + "\n")
        pepList = list(refPep)
        for pep in pepList:
            out.write(pep + "\n")
        out.close()

    file0.close()

    cmd6 = "rm myTmp0  myTmp1"
    os.system(cmd6)

if __name__ == "__main__":
    main()

