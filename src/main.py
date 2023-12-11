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
   
def evaluation(predData,evaData,rnd,eval):
    if len(evaData)==0:
        return
    for feature in range(len(predData)):
        candidates=predData[feature].get_recommendation()
        winner=[]
        times=[]
        for candidate in candidates:
            abundance=[]
            for eva_rnd in evaData:
                number=eva_rnd.get_abundance(candidate)
                abundance.append(number)
            result=True
            '''
            for x in range(len(abundance)-1):
                if abundance[x]>abundance[x+1]:           
            '''
            if abundance[0]>abundance[1]:
                result=False
            if result==True:
                winner.append(candidate) #Test
                if len(abundance)==3:
                    times.append([abundance[0],abundance[1],abundance[2]])
                else:
                    times.append([abundance[0],abundance[1]])
                #winner.append(candidate)
        print(f"The final recommendation for {rnd[feature]} is:")
        if len(winner)>0:
            print(f"There are {len(winner)} aptamer for this target. They are")    
            for x in range(len(winner)):
                color=hash(winner[x])%8+(3*hash(winner[x]))%2*60+30
                if len(times[x])==2:
                    print(f"\033[{color}m{winner[x]}\033[0m,which appears {times[x][0]} times in {eval[0]}, {times[x][1]} times in {eval[1]}")
                else:
                    print(f"\033[{color}m{winner[x]}\033[0m,which appears {times[x][0]} times in {eval[0]}, {times[x][1]} times in {eval[1]}, {times[x][2]} times in {eval[2]}")
                # print(f"3[{color}m{winner[x]}[0m")
                # print("3[30mSuixinBlog: https://suixinblog.cn3[0m")
        else:
            print("We didn't find any aptamer based on this evaluation")
        df=pd.DataFrame(winner, columns=['recommendation']) #Test
        #df=pd.DataFramce(winner)
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Recommendation after evaluation for target {rnd[feature]}.csv', index=False)
    
if __name__ == '__main__': 
    # parameter= [lenofseq,goal,num1,num2,op] 
    parameter=[36,20,50,1000,False]
    pred=['caff20']
    eval=['caff12','caff20']
    predData=[]
    evaData=[]
    for x in pred:
        clear_temp_files()
        data=Round(x)
        data.set_para(parameter[0],parameter[1],parameter[2],parameter[3],parameter[4])
        data.set_data()
        predData.append(data)
        
    for x in eval:
        clear_temp_files()
        finaldata=Round(x)
        finaldata.set_para(parameter[0],parameter[1],parameter[2],parameter[3],parameter[4])
        finaldata.set_seq_abundance()
        evaData.append(finaldata)
    evaluation(predData,evaData,pred,eval)
