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
   
def evaluation(predData,evaData,rnd):
    for feature in range(len(predData)):
        candidates=predData[feature].get_recommendation()
        winner=[]
        for candidate in candidates:
            abundance=[]
            for eva_rnd in evaData:
                number=eva_rnd.get_abundance(candidate)
                abundance.append(number)
            result=True
            for x in range(len(abundance)-1):
                if abundance[x]>abundance[x+1]:
                    result=False
            if result==True:
                winner.append(candidate) #Test
                #winner.append(candidate)
        print(f"The final recommendation for {rnd[feature]} is:")
        if len(winner)>0:
            print(f"There are {len(winner)} aptamer for this target. They are")    
            for x in range(len(winner)):
                color=hash(winner[x][0])%8+(3*hash(winner[x][0]))%2*60+30
                print(f"\033[{color}m{winner[x]}\033[0m")
                # print(f"3[{color}m{winner[x]}[0m")
                # print("3[30mSuixinBlog: https://suixinblog.cn3[0m")
        else:
            print("We didn't find any aptamer based on this evaluation")
        df=pd.DataFrame(winner, columns=['recommendation']) #Test
        #df=pd.DataFramce(winner)
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Recommendation after evaluation for target {rnd[feature]}.csv', index=False)
    
if __name__ == '__main__': 
    # parameter= [lenofseq,goal,num1,num2] 
    # 起始值与文件中的primer有关
    # lenofseq: 序列长度
    # goal: 选取的序列数量
    # num1: kmerlen<8时的统计数量，推荐为40
    # num2: kmerlen>8时的统计数量，推荐为1000
    parameter=[46,10,50,1000]
    #parameter=['with_primer',20,4,10]  # Test
    rnd=['atpr10']
    final=[]
    predData=[]
    evaData=[]
    for x in rnd:
        clear_temp_files()
        data=Round(x)
        data.set_para(parameter[0],parameter[1],parameter[2],parameter[3])
        data.set_data()
        predData.append(data)
        
    for x in final:
        clear_temp_files()
        finaldata=Round(x)
        finaldata.set_para(parameter[0],parameter[1],parameter[2],parameter[3])
        finaldata.set_seq_abundance()
        evaData.append(finaldata)
    evaluation(predData,evaData,rnd)