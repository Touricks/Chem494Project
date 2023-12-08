from Round import Round
from Sequence import Sequence
from Cluster import Cluster
import subprocess
import pandas as pd
import os

def clear_temp_files():
    cmd="rm temp_ct.txt"
    subprocess.call([cmd], shell=True)
    cmd="rm temp_seqs.txt"
    subprocess.call([cmd], shell=True)

def similarity(str1:str,str2:str,length):
        dis=0
        for x in range(len(str1)):
            if str1[x]!=str2[x]:
                dis+=1
        if dis/length<=0.3:
            return True
        else:
            return False
        
def evaluation_ed2(predData,rnd,seqlen):
    rnd1=predData[0].get_cludata()
    rnd2=predData[1].get_cludata()
    data=[]
    filterdata=[]
    for feature in rnd2:
        if (feature in rnd1)and(rnd1[feature]['relcount']<rnd2[feature]['relcount']):
            data.append([feature,rnd2[feature]['relcount'],rnd1[feature]['relcount'],
                        round(rnd2[feature]['relcount']/rnd1[feature]['relcount'],5),
                         rnd2[feature]['represent'],rnd2[feature]['mfestructure']])
            filterdata.append([feature,rnd1[feature]['represent'],rnd1[feature]['count'],rnd2[feature]['count'],0])
    data.sort(key=lambda x:x[3],reverse=True)
    #print(data)
    print(predData[0].cludata.count_for_test)
    print(predData[1].cludata.count_for_test)
    df=pd.DataFrame(data, columns=['feature',rnd[1],rnd[0],'ratio','seq','ct'])
    path= str(os.getcwd()) + "/dataset"
    df.to_csv(path+f'/Evaluation_from_{rnd[0]}_to_{rnd[1]}.csv', index=False)

    for x in range(len(filterdata)-2):
        if filterdata[x][4]==0:
            for y in range(x+1,len(filterdata)):
                if similarity(filterdata[x][1],filterdata[y][1],seqlen):
                    filterdata[x][2]+=filterdata[y][2]
                    filterdata[x][3]+=filterdata[y][3]
                    filterdata[y][4]=1
    filterdata=list(filter(lambda x:x[4]==0,filterdata))
    for x in filterdata:
        x[4]=round(x[3]/x[2],5)
    filterdata.sort(key=lambda x:x[4],reverse=True)
    df=pd.DataFrame(filterdata, columns=['feature','seq','count1','count2','ratio'])
    path= str(os.getcwd()) + "/dataset"
    df.to_csv(path+f'/Simplified_evaluation_from_{rnd[0]}_to_{rnd[1]}.csv', index=False)

def evaluation_ed3(predData,rnd,seqlen):
    rnd1=predData[0].get_cludata()
    rnd2=predData[1].get_cludata()
    rnd3=predData[2].get_cludata()
    data=[]
    filterdata=[]
    for feature in rnd3:
        if (feature in rnd1)and(feature in rnd2)and(rnd1[feature]['relcount']<rnd2[feature]['relcount'])and(rnd2[feature]['relcount']<rnd3[feature]['relcount']):
            data.append([feature,rnd3[feature]['count'],rnd2[feature]['count'], rnd1[feature]['count'],
                        rnd3[feature]['relcount']/rnd2[feature]['relcount'], rnd2[feature]['relcount']/rnd1[feature]['relcount'],
                        rnd3[feature]['relcount']/rnd1[feature]['relcount'],rnd2[feature]['represent'],rnd2[feature]['mfestructure']])
            filterdata.append([feature,rnd1[feature]['represent'],rnd1[feature]['count'],rnd2[feature]['count'],rnd3[feature]['count'],0])
    data.sort(key=lambda x:x[4]*x[5],reverse=True)
    #print(data)
    print(predData[0].cludata.count_for_test)
    print(predData[1].cludata.count_for_test)
    print(predData[2].cludata.count_for_test)
    df=pd.DataFrame(data, columns=['feature',rnd[2],rnd[1],rnd[0],'ratio32','ratio21','ratio31','seq','ct'])
    path= str(os.getcwd()) + "/dataset"
    df.to_csv(path+f'/Evaluation_from_{rnd[0]}_through_{rnd[1]}_to_{rnd[2]}.csv', index=False)

    for x in range(len(filterdata)-2):
        if filterdata[x][5]==0:
            for y in range(x+1,len(filterdata)):
                if similarity(filterdata[x][1],filterdata[y][1],seqlen):
                    filterdata[x][2]+=filterdata[y][2]
                    filterdata[x][3]+=filterdata[y][3]
                    filterdata[x][4]+=filterdata[y][4]
                    filterdata[y][5]=1
    filterdata=list(filter(lambda x:x[5]==0,filterdata))
    for x in filterdata:
        x[5]=round(x[4]/x[2],5)
    filterdata.sort(key=lambda x:x[5],reverse=True)
    df=pd.DataFrame(filterdata, columns=['feature','seq','count1','count2','count3','ratio'])
    path= str(os.getcwd()) + "/dataset"
    df.to_csv(path+f'/Simplified_evaluation_from_{rnd[0]}_through_{rnd[1]}_to_{rnd[2]}.csv', index=False)

if __name__ == '__main__': 
    # parameter= [lenofseq,goal，num1,num2] 
    # lenofseq: 序列长度有多少
    # goal: 要选取的序列数量
    # num1: kmerlen<8时的统计数量，推荐为40
    # num2: kmerlen>8时的统计数量，推荐为1000
    parameter=[46,10,50,1000]
    rnd=['atpr6','atpr9','atpr10']
    predData=[]
    for x in rnd:
        clear_temp_files()
        data=Round(x)
        data.set_para(parameter[0],parameter[1],parameter[2],parameter[3])
        data.set_data_ed2()
        predData.append(data)
    #evaluation_ed2(predData,rnd,parameter[0])
    evaluation_ed3(predData,rnd,parameter[0])