from Aptamer import Aptamer
from Kmer import Kmer
from Sequence import Sequence
from Cluster import Cluster
import os,subprocess,multiprocessing
from pyfastaq import sequences
import pandas as pd
import math

class Round():
    def __init__(self,binding_target=None,round=None):
        #Step1:initialize Round instance
        #Step2:Run set_para()
        #Step3:Run set_data()
        self.binding_target=binding_target
        self.round=round
        self.seqdata=[]
        self.kmerdata = None
        self.cludata = None
        self.weighted = True
        self.mode= None
        self.result=[]
        
    def set_para(self,weighted,op,goal):
        self.weighted = weighted
        self.mode = op
        self.goal= goal

    def set_data(self):
        if self.mode=='with_primer':
            seq_list=self.__get_sequence_with_primer(self.binding_target,self.round) #Build seq_list from fasta file
        else:
            seq_list=self.__get_sequence_without_primer() 
        count_seq=self.__get_count_single_rnd(seq_list) #Build count_seq from seq_list to gather count data
        self.__add_RNAfold_data(count_seq) #Add RNAfold data to count_seq, gather dg/ct data
        self.kmerdata=Kmer(self.seqdata,self.weighted,self.round)# Build kmerdata from seqdata, initialize with instance of Kmer
        self.kmerdata.set_candidate()# Set candidate Eigenkmers for the execution of method get_kmer_list()
        kmerlist=self.kmerdata.get_kmer_list()# Get kmerlist from kmerdata
        self.cludata=Cluster(self.seqdata,self.weighted,kmerlist,self.round,self.goal) #Build cludata from seqdata, initialize with instance of Cluster
        self.cludata.set_cluster_data() #Make cluster based on seqdata, then gather size/kmer data for each cluster
        self.result=self.cludata.get_recommend_aptamer()
        self.__result_display()
        self.__record_sequence_data()
    
    def get_sequence_list(self):
        return self.seqdata
    
    def __get_sequence_with_primer(self,binding_target,rnd):
        # get_sequence: 从fastq文件中提取序列
        # binding_target: 本次实验的靶标
        # rnd: 本次实验的轮次
        # return: 对应Selex轮次中所有序列组成的列表
        current_file_path = os.path.abspath(__file__)
        data_file=current_file_path.split('src')[0]+'data'
        path_seq = data_file+f"/{binding_target}_fastq_r2"
        path_primer = data_file+f"/{binding_target[:4]}_primers.txt"
        with open(path_primer,"r") as f:
            primer = f.readlines()
        f.close()
        df = []
        for file in os.listdir(path_seq):
            if file.replace(binding_target[:4], "").replace(".fastq", "") == rnd:
                seq_reader = sequences.file_reader(path_seq + "/" + file)
                for sequence in seq_reader:
                    seq, succ = self.__getN30(sequence.seq,primer)
                    if succ == False:
                        continue
                    #seq=sequence.seq[18:48]
                    #if seq=='GATTGTGGTCTATTCATAGGCGTCCGCTGA':
                    #    print(f"find this in rnd {rnd}")
                    df.append(seq)
                break
        return df
    
    def __get_sequence_without_primer(self):
        current_file_path = os.path.abspath(__file__)
        data_file=current_file_path.split('src')[0]+'data'
        #data_file=current_file_path.split('src')[0]+'src'+current_file_path.split('src')[1]+'data'
        path_seq=data_file+'/'+str(os.listdir(data_file)[0])
        df = []
        print(f"HEREIS {path_seq}")
        seq_reader = sequences.file_reader(path_seq)
        for sequence in seq_reader:
            seq=sequence.seq[18:48]
            df.append(seq)
        return df

    def __getN30(self,seq,primer):
        # getN30: 通过primer信息从序列seq中提取N30
        # seq: 本次实验的序列
        # primer_path: 本次实验的引物文件路径
        # return: N30序列
        # 通过实验验证第二段底物在序列数据中完全不存在(一个都没有)，
        # 第一段底物出现概率为88% （且各轮数据均相同），因此策略为
        # 不统计第一段底物不出现的序列，忽略第二段底物，直接从第一段底物后面取30bp
        pos=seq.find(primer[0].strip())
        if pos==-1:
            success=False
            return None, success
        start=pos+len(primer[0].strip())
        N30=seq[start:start+30]
        success = True
        return N30, success
    
    def __get_count_single_rnd(self,seq_list):
        # get_count_single_rnd: 计算单轮数据中每个序列的出现次数
        # seq_list: 单轮数据的序列列表
        # return: 单轮数据中每个序列的出现次数
        df=pd.DataFrame(seq_list,columns=['seq'])
        count_df=df.groupby('seq',as_index=False).value_counts(sort=True,ascending=False)
        seq_sum=count_df['count'].sum()
        count_df['freq'] = count_df['count'] / seq_sum
        count_df = count_df.sort_values(by='freq', ascending=False)
        freq=count_df['freq'].tolist()
        cum_freq = [freq[0]]
        for i in range(1, len(freq)):
            cum_freq.append(cum_freq[i-1]+freq[i])
        count_df['cum_freq'] = cum_freq
        length=len(count_df)
        count_df.index=list(range(1,length+1))
        print(count_df)
        return count_df
    
    def __add_RNAfold_data(self,count_df):
        # add_RNAfold_data: 将RNAfold计算得到的二级结构和自由能添加到库中
        # binding_target: 本次实验的靶标
        # rnd: 本次实验的轮次
        # count_df: 本次实验的序列计数数据
        # return: 无
        seq_list=count_df['seq'].tolist() #将序列转化为列表
        count_list=count_df['count'].tolist() #将出现次数转化为列表
        #造ct，dg数据
        seq_file = str(os.getcwd()) + "/temp_seqs.txt"
        with open(seq_file, 'w') as f:
            for x in seq_list:
                f.write(x + "\n")
        f.close()
        ct_file = str(os.getcwd()) + "/temp_ct.txt"
        cmd = f"RNAfold --jobs={multiprocessing.cpu_count()} --infile={seq_file} --outfile={ct_file.split('/')[-1]} --noPS -T {37.0} --noconv"
        subprocess.call([cmd], shell=True)
        ct_list=[]
        dg_list=[]
        with open(ct_file, 'r') as f:
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
            sequence=Sequence(apt,count_list[x],self.weighted)
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
        for seq in self.result:
            print(f"No.{cnt+1} sequence: {seq.get_aptamer_seq()} | ct: {seq.get_aptamer_ct()} | Count: {seq.get_aptamer_count()} | dG: {seq.get_aptamer_dG()}")
            cnt+=1
    
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
        path= str(os.getcwd()) + "/testexcel"
        df.to_csv(path+f'/Sequence_Data_for_{self.round}_round.csv', index=False)