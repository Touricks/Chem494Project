from Aptamer import Aptamer
from Kmer import Kmer
from Sequence import Sequence
from Cluster import Cluster
import os,subprocess,multiprocessing
from pyfastaq import sequences
import pandas as pd
import math
import hashlib
import shutil

class Round():
    def __init__(self,round=None):
        #Step1:initialize Round instance
        #Step2:Run set_para()
        #Step3:Run set_data()
        self.round=round
        self.seqdata=[]
        self.weightdata={}
        self.kmerdata = None
        self.cludata = None
        self.lenofseq= 0
        self.result=[]
        self.recommend=None
        self.amount_seq=0
        
    def set_para(self,op,goal,num1,num2,switch):
        self.lenofseq = op
        self.goal= goal
        self.num1=num1
        self.num2=num2
        self.switch=switch

    def set_data(self):
        seq_list=self.__get_sequence_with_primer(self.round,'prediction') #Build seq_list from fasta file
        count_seq=self.__get_count_single_rnd(seq_list) #Build count_seq from seq_list to gather count data
        self.__add_RNAfold_data(count_seq) #Add RNAfold data to count_seq, gather dg/ct data, making self.seqdata
        self.kmerdata=Kmer(self.seqdata,self.round,self.num1,self.num2)# Build kmerdata from seqdata, initialize with instance of Kmer
        self.kmerdata.set_candidate()# Set candidate Eigenkmers for the execution of method get_kmer_list()
        kmerlist=self.kmerdata.get_kmer_list()# Get kmerlist from kmerdata
        self.cludata=Cluster(self.seqdata,kmerlist,self.round,self.goal,self.amount_seq,self.num1**2*2+self.num2*5) #Build cludata from seqdata, initialize with instance of Cluster
        self.cludata.set_cluster_data() #Make cluster based on seqdata, then gather size/kmer data for each cluster
        self.result=self.cludata.get_recommend_aptamer()
        self.__result_display()
        self.__record_sequence_data()
    
    def set_data_ed2(self):
        '''
        This time, we use entire feature as indication
        '''
        seq_list=self.__get_sequence_with_primer(self.round,'prediction') #Build seq_list from fasta file
        print('finished reading')
        count_seq=self.__get_count_single_rnd(seq_list) #Build count_seq from seq_list to gather count data
        print('finished counting')
        self.__add_RNAfold_data(count_seq) #Add RNAfold data to count_seq, gather dg/ct data
        print('finished adding RNAfold')
        self.cludata=Cluster(self.seqdata,None,self.round,self.goal,self.amount_seq,self.num1**2*2+self.num2*5)
        self.cludata.set_cluster_data2()
        print('finished clustering')

    def set_seq_abundance(self):
        seq_list=self.__get_sequence_with_primer(self.round,'evaluation') #Build seq_list from fasta file
        count_seq=self.__get_count_single_rnd(seq_list) #Build count_seq from seq_list to gather count data
        seq_list=count_seq['seq'].tolist() #将序列转化为列表
        count_list=count_seq['count'].tolist() #将出现次数转化为列表
        for x in range(len(seq_list)):
            seqid = hashlib.sha1(str(seq_list[x]).encode('utf8')).hexdigest()[:10]
            if seqid not in self.weightdata.keys():
                self.weightdata[seqid]=count_list[x]
        data=[]
        for x in range(len(seq_list)):
            data.append([seq_list[x],count_list[x]])
        df=pd.DataFrame(data, columns=['sequence','count'])
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Abundance_Data_for_{self.round}.csv', index=False)

        
    def get_sequence_list(self):
        return self.seqdata
    
    def get_abundance(self,seq):
        seqid = hashlib.sha1(str(seq).encode('utf8')).hexdigest()[:10]
        if seqid in self.weightdata.keys():
            return self.weightdata[seqid]
        else:
            return 0
    
    def get_recommendation(self):
        return self.recommend
    
    def __get_sequence_with_primer(self,rnd,foldername):
        current_file_path = os.path.abspath(__file__)
        print(current_file_path)
        data_file=current_file_path.split('src')[0]+'data/'
        path_seq = data_file+f'{foldername}/'
        path_primer = data_file+f"primer.txt"
        with open(path_primer,"r") as f:
            primer = f.readlines()
        f.close()
        df = []
        for file in os.listdir(path_seq):
            if file.replace(".fastq", "") == rnd:
                seq_reader = sequences.file_reader(path_seq + "/" + file)
                for sequence in seq_reader:
                    seq, succ = self.__getN30(sequence.seq,primer)
                    if succ == False:
                        continue
                    df.append(seq)
                break
        return df

    def __getN30(self,seq,primer):
        pos=seq.find(primer[0].strip())
        if pos==-1:
            success=False
            return None, success    
        start=pos+len(primer[0].strip())-6
        N30=seq[start:start+self.lenofseq+12]
        success = True
        return N30, success
    
    def __get_count_single_rnd(self,seq_list):
        df=pd.DataFrame(seq_list,columns=['seq'])
        count_df=df.groupby('seq',as_index=False).value_counts(sort=True,ascending=False)
        seq_sum=count_df['count'].sum()
        self.amount_seq=seq_sum
        count_df['freq'] = count_df['count'] / seq_sum
        count_df = count_df.sort_values(by='freq', ascending=False)
        freq=count_df['freq'].tolist()
        cum_freq = [freq[0]]
        for i in range(1, len(freq)):
            cum_freq.append(cum_freq[i-1]+freq[i])
        count_df['cum_freq'] = cum_freq
        length=len(count_df)
        count_df.index=list(range(1,length+1))
        print(f"The amount of sequence is {seq_sum}")
        print(count_df)
        return count_df
    
    def __add_RNAfold_data(self,count_df):

        print(f"COUNT_DF LENGTH=",len(count_df))
        seq_list=count_df['seq'].tolist() 
        print(f"SEQ_LIST LENGTH=",len(seq_list))
        count_list=count_df['count'].tolist() 

        flag=0
        new_ct_file=None

        path_ct = str(os.getcwd())+f"/ctfiledatabase/"
        for file in os.listdir(path_ct):
            if file.replace("_ct.txt", "") == self.round:
                new_ct_file=path_ct+file
                flag=1
                print('ct file already exists')
                break
        if (self.switch==False)or(flag==0):
            seq_file = str(os.getcwd()) + "/temp_seqs.txt"
            with open(seq_file, 'w') as f:
                for x in seq_list:
                    f.write(x + "\n")
            f.close()
            ct_file="temp_ct.txt"
            new_ct_file=str(os.getcwd()) + f"/ctfiledatabase/{self.round}_ct.txt"
            cmd = f"RNAfold --jobs={multiprocessing.cpu_count()} --infile={seq_file} --outfile={ct_file.split('/')[-1]} --noPS -T {37.0} --noconv"
            subprocess.call([cmd], shell=True)
            shutil.move(ct_file, new_ct_file)
        ct_list=[]
        dg_list=[]
        with open(new_ct_file, 'r') as f:
            for i, line in enumerate(f):
                if i % 2 == 0:
                    sequence = line.strip()
                else: 
                    ct = line.strip().split(" ")[0]
                    dg = float(line.strip().split(" (")[-1].replace("(", "").replace(")", ""))
                    ct_list.append(ct)
                    dg_list.append(dg)
            f.close()
        for x in range(len(seq_list)):
            apt=Aptamer(seq_list[x],ct_list[x],dg_list[x])
            sequence=Sequence(apt,count_list[x])
            self.seqdata.append(sequence)
    
    def __result_display(self):
        #Result display part
        print("Cluster display:")
        if self.cludata is not None:
            self.cludata.cluster_display()
        else:
            print("No cluster data available.")
        print("Recommend sequence:")
        cnt=0
        data=[]
        data1=[]
        for seq in self.result:
            content=f"No.{cnt+1} sequence: {seq.get_aptamer_seq()} | ct: {seq.get_aptamer_ct()} | Count: {seq.get_aptamer_count()} | dG: {seq.get_aptamer_dG()}"
            print(content)
            data.append(seq.get_aptamer_seq())
            data1.append(content)
            cnt+=1
        df=pd.DataFrame(data1, columns=['recommendation'])
        self.recommend=data
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Recommendation for target {self.round}.csv', index=False)
                  
    def __record_sequence_data(self):
        #Display the most abundance sequence with their cluster
        print("Most abundance sequence:")
        cnt=0
        for seq in self.seqdata:
            if cnt==self.goal:
                break
            Numofclu=seq.get_cluster_num()
            print(f"No.{cnt+1} sequence: {seq.get_aptamer_seq()} | ct: {seq.get_aptamer_ct()} | Count:{seq.get_aptamer_count()} | affiliated: {Numofclu}")
            cnt+=1
        #Test
        data=[]
        for seq in self.seqdata:
            data.append([seq.get_aptamer_seq(),seq.get_aptamer_count(),seq.get_aptamer_ct(),seq.get_cluster_num()])
        df=pd.DataFrame(data, columns=['sequence','counts','ct','cluster'])
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Sequence_Data_for_{self.round}_round.csv', index=False)
