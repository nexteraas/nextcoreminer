#!/bin/python
import argparse
import os

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='?',help='The input should be a list of peptides')
    my_parser.add_argument('-a',  nargs='?', default="HLA-DQA10501-DQB10201",help='HLA allele')
    return my_parser.parse_args()
    
def main():

    args = parse_arguments()
    refArgs = vars(args)

    myFile = refArgs.get("f") 

    netMHC2pan = "/home/galaxy/galaxy/tools/testing/netMHCIIpan-4.0/netMHCIIpan"
    allele = refArgs.get("a")

    cmd0="cut -f 1 " + myFile + " |awk '$0 ~ /^[A-Z]+$/' > myTmp0"
    cmd1=netMHC2pan + " -f myTmp0 -inptype 1 -a " + allele + " > myTmpOut"
    os.system(cmd0)
    os.system(cmd1)
    cmd2='''gawk 'length($0)>1' myTmpOut |gawk '$0!~"---"' |grep -v "#"  > myNetMHC2pan.txt'''
    os.system(cmd2)

    file0=open("myNetMHC2pan.txt")
    outCount=open("myNetMHC2pan.BinderCount.txt" ,"w")
    count_sb=0
    count_wb_sb=0
    count_all=0
    outMotif=open("myNetMHC2pan.9mer.txt" ,"w")
    while True:
        line0=file0.readline().rstrip()
        if (len(line0)==0):
            break
        arr=line0.split()
    
        if (arr[0]!="Seq" and arr[0]!="Number"):
            count_all = count_all + 1
            if ("SB" in arr[-1]):
                count_sb=count_sb + 1
                count_wb_sb=count_wb_sb + 1
                outMotif.write(arr[4] + "\n")
            if ("WB" in arr[-1]):
                count_wb_sb=count_wb_sb + 1
                outMotif.write(arr[4] + "\n")
    file0.close()
    outCount.write("Total:" + str(count_wb_sb) + "\t" + "SB:" + str(count_sb) + "\t" + "WB:" + str(count_wb_sb - count_sb) + "\n")
    outCount.close()
    outMotif.close()

    cmd3 = "rm myTmp0 myTmpOut "
    os.system(cmd3)

if __name__ == "__main__":
    main()

