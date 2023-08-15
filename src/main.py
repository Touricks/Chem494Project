from Round import Round
from Sequence import Sequence
from Cluster import Cluster
import subprocess
import pandas as pd
import os
import sys

def clear_temp_files():
    cmd="rm temp_ct.txt"
    subprocess.call([cmd], shell=True)
    cmd="rm temp_seqs.txt"
    subprocess.call([cmd], shell=True)

def get_info(seq:Sequence):
    return [seq.get_aptamer_seq(),seq.get_aptamer_count()]

def display_comparison_result():
    Finallist=list(map(get_info,finaldata.get_sequence_list()))
    Finallist=Finallist[:300]
    for round in range(len(roundData)):
        data=[]
        for seq in roundData[round].get_sequence_list():
            for x in Finallist:
                if seq.get_aptamer_seq()==x[0]:
                    clunum=seq.get_cluster_num()
                    if clunum==-1:
                        kmer1='None'
                        kmer2='None'
                        size=0
                    else:
                        cluster=roundData[round].cludata
                        theHead=cluster.get_cluster_head(clunum)
                        kmer1=theHead.eigenkmer1
                        kmer2=theHead.eigenkmer2
                        size=theHead.get_size()
                    data.append([seq.get_aptamer_seq(),seq.get_aptamer_count(),seq.get_aptamer_ct(),clunum,kmer1,kmer2,size,x[1]])
                    break
        df=pd.DataFrame(data,columns=['Sequence','Count in curr rnd','ct','ClusterNum','Kmer1','Kmer2','SizeofClu','Count in Final rnd'])
        path= str(os.getcwd()) + "/testexcel"
        df.to_csv(path+f'/Round {rnd[round]} sequences which appears in Round {final[0]}.csv', index=False)

if __name__ == '__main__': 
    parameter=[True,None,None]
    syspara=sys.argv[1:]
    if len(syspara)==2:
        if syspara[0]=='0':
            parameter[1]='with_primer'
        else:
            parameter[1]='without_primer'
        parameter[2]=int(syspara[1])
    elif len(syspara)==1:
        if syspara[0]=='0':
            parameter[1]='with_primer'
        else:
            parameter[1]='without_primer'
        parameter[2]=15
    elif len(syspara)==0:
        parameter[1]='with_primer'
        parameter[2]=15
    else:
        print('ERROR, Number of argument is wrong')
        sys.exit()

    if parameter[1]=='with_primer':
        binding_target='theophylline'
        rnd=['15','16','20']
        final=['22']
        roundData=[]
        for x in final:
            clear_temp_files()
            finaldata=Round(binding_target,x)
            finaldata.set_para(parameter[0],parameter[1],parameter[2])
            finaldata.set_data()
        for x in rnd:
            clear_temp_files()
            data=Round(binding_target,x)
            data.set_para(parameter[0],parameter[1],parameter[2])
            data.set_data()
            roundData.append(data)
        display_comparison_result()
    else:
        clear_temp_files()
        data=Round(None,None)
        data.set_para(parameter[0],parameter[1],parameter[2])
        data.set_data()