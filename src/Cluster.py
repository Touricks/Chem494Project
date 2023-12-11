from Sequence import Sequence
import math
import pandas as pd
import os
import sys
from scipy.stats import poisson
import numpy as np

class Member():
    def __init__(self,seq:Sequence,kmer1:str,kmer2:str,N:int):
        '''
        self.size: Number of aptamers in the cluster
        '''
        self.eigenkmer1=kmer1
        self.eigenkmer2=kmer2
        self.score=0
        self.rank=0
        self.merged=False
        self.family=[seq]
        self.N=N
        self.size=seq.get_aptamer_count()
    
    def calculate_score(self):
        s=self.get_size()
        x=5
        if self.eigenkmer2=='':
            x=len(self.eigenkmer1)-8
        f1=[1,1,1,1,1,4] #N=20
        f2=[1,4,16,64,256,128]
        kmerscore=list(map(lambda x,y:x/y,f1,f2))
        self.score=s/kmerscore[x]

    def equal(self,kmer1,kmer2):
        if (self.eigenkmer1==kmer1)and(self.eigenkmer2==kmer2):
            return True
        else:
            return False
        
    def add_family(self,seq:Sequence):
        self.family.append(seq)
        self.size+=seq.get_aptamer_count()

    def get_size(self):
        return self.size
    
    def get_score(self):
        return self.score
    
    def set_rank(self,rank):
        self.rank=rank
        for seq in self.family:
            if seq.get_cluster_num()==-1:
                seq.set_cluster_num(rank)

class Member2():
    def __init__(self,feature,count,seq):
        self.feature=feature
        self.count=count
        self.family=[seq]
        self.relcount=0

    def addcount(self,count):
        self.count+=count
    
    def setrelcount(self,N):
        self.relcount=self.count/N


class Cluster():
    def __init__(self,seqdata:list[Sequence],kmerlist,round,goal,N,kmeramount):
        '''
        N:序列的总数量
        '''
        self.head=[]
        self.head2=[]
        self.head2feature=[]
        self.data=seqdata
        self.kmers=kmerlist
        self.size=0 #number of clusters
        self.testdata=[]
        self.round=round
        self.goal=goal
        self.N=N
        self.kmercnt=kmeramount

    def set_cluster_data(self):
        self.__get_cluster()
        print("Clustering finished, sorting...")
        self.__sort_cluster()
        print("Combining clusters...")
        self.__combine_cluster()
        print("Combining finished, predicting...")
        self.__sort_family()
        self.__set_rank()    
        self.__record_cluster_data()

    def issimilar(self,feature:list[str],featurelist:list[list[str]]):
        for x in featurelist:
            #print(feature,x)
            str1=''.join(feature)
            str2=''.join(x)
            #print(str1,str2)
            if self.__get_distance2(str1,str2)<0.3:
                return True,x
        return False,0
    
    def set_cluster_data2(self):
        Clu_cnt=0 #number of clusters
        Clu_seq_cnt=0 #number of sequences in clusters
        Seq_cnt=0 #number of testing sequences 
        for seq in self.data:
            Seq_cnt+=1
            feature=seq.get_kmer_list2()
            flag,clufeature=self.issimilar(feature,self.head2feature)
            if flag:
                Clu_seq_cnt+=1
                index=self.head2feature.index(clufeature)
                self.head2[index].addcount(seq.get_aptamer_count())
            elif Clu_cnt<100:
                self.head2.append(Member2(feature,seq.get_aptamer_count(),seq))
                self.head2feature.append(feature)
                Clu_cnt+=1
                Clu_seq_cnt+=1
            if Clu_cnt%10==0:
                print(Clu_cnt,Clu_seq_cnt,Seq_cnt)
        for x in self.head2:
            x.setrelcount(self.N)

    def get_recommend_aptamer(self):
        cnt=0
        recommend=[]
        for head in self.head:
            if cnt==self.goal:
                break
            candidate=head.family
            if not(candidate[0] in recommend):
                recommend.append(candidate[0])
                cnt+=1
        return recommend

    def get_cluster_head(self,rank):
        return self.head[rank]

    def get_cluster_head2(self):
        return self.head2
        
    def cluster_display(self):
        """
        Prints information about each cluster, including the cluster number, the number of aptamers in the cluster,
        and the head aptamer sequence. For each cluster, it also prints information about each member aptamer,
        including the aptamer sequence and its count.
        """
        cnt=1
        for head in self.head:
            print(f'Cluster {cnt}:, Eigenkmer: {head.eigenkmer1} {head.eigenkmer2}, Cluster size: {head.get_size()}, Cluster score: {head.get_score()}')
            cnt+=1
            tot=1
            print('The first 10 members with maximum abundance aptamer in this cluster:')
            for member in head.family:
                print(f'    Seq{tot}: {member.get_aptamer_seq()} | count: {member.get_aptamer_count()} | ct={member.get_aptamer_ct()} | dG= {member.get_aptamer_dG()}')
                tot+=1
                if tot==10:
                    break
            if cnt==10:
                break

    def __get_cluster(self):
        flag=1 # Handle eigenkmer with k>=8
        Kmer_cnt=0
        for x in range(len(self.kmers)):
            Kmer_cnt+=1
            if Kmer_cnt%100==0:
                print(f"{Kmer_cnt}/{self.kmercnt}")
            if len(self.kmers[x])==1:
                self.__set_longkmer(self.kmers[x])
            elif flag: # Reset the state of all sequences (If they have been clustered, set state to False)
                self.__renew_seqrep() 
                self.__set_shortkmer(self.kmers[x])
                flag=0
            else:
                self.__set_shortkmer(self.kmers[x])

    def __set_longkmer(self,eigenkmer):
        eigenkmer=eigenkmer[0]
        for seq in self.data:
            if seq.get_state()==True:
                continue
            kmer,pos,kmer_len=seq.apt.get_kmer()
            if eigenkmer in kmer:
                self.__add_cluster_longkmer(seq,eigenkmer)
                seq.set_state(True)

    def __set_shortkmer(self,eigenkmer):
        kmer1=eigenkmer[0]
        kmer2=eigenkmer[1]
        for seq in self.data:
            if seq.get_state()==True:
                continue
            kmer,pos,kmer_len=seq.apt.get_kmer() 
            if (kmer1 in kmer)and(kmer2 in kmer):
                if self.__is_valid(kmer,pos,kmer_len,kmer1,kmer2):
                    self.__add_cluster_shortkmer(seq,kmer1,kmer2)
                    seq.set_state(True)

    def __add_cluster_longkmer(self,seq:Sequence,eigenkmer):
        if (self.head==[]):
            newhead=Member(seq,eigenkmer,'',self.N)
            self.head.append(newhead)
            self.size+=1
        else:
            member=self.head[-1]
            if member.equal(eigenkmer,''):
                member.add_family(seq)
            else:
                newhead=Member(seq,eigenkmer,'',self.N)
                self.head.append(newhead)
                self.size+=1
                
    def __add_cluster_shortkmer(self,seq:Sequence,x1:str,x2:str):
        if (self.head==[]):
            newhead=Member(seq,x1,x2,self.N)
            self.head.append(newhead)
            self.size+=1
        else:
            member=self.head[-1]
            if member.equal(x1,x2):
                member.add_family(seq)
            else:
                newhead=Member(seq,x1,x2,self.N)
                self.head.append(newhead)
                self.size+=1
                
    def __renew_seqrep(self):
        for seq in self.data:
            seq.set_state(False)

    def __sort_cluster(self):
        for x in self.head:
            x.calculate_score()
        self.head.sort(key = lambda x: x.score,reverse=True)   
    
    def __combine_cluster(self):
        print(len(self.head))
        for x in range(len(self.head)):
            if self.head[x].merged==False:
                for y in range(x+1,len(self.head)):
                    if self.head[y].merged==False:
                        if self.__similarity(self.head[x],self.head[y]):
                            self.testdata.append([self.head[x].eigenkmer1,self.head[x].eigenkmer2,self.head[x].score,self.head[x].size,
                                                  self.head[y].eigenkmer1,self.head[y].eigenkmer2,self.head[y].score,self.head[y].size])
                            for seq in self.head[y].family:
                                self.head[x].add_family(seq)
                            self.head[x].score+=self.head[y].score
                            self.head[y].merged=True      
        newhead=[]
        for x in range(len(self.head)):
            if self.head[x].merged==False:
                newhead.append(self.head[x])
        self.head=newhead

    def __similarity(self,mem1,mem2):
        eigenkmer1_1=mem1.eigenkmer1
        eigenkmer1_2=mem1.eigenkmer2
        eigenkmer2_1=mem2.eigenkmer1
        eigenkmer2_2=mem2.eigenkmer2
        if abs(len(eigenkmer1_2)-len(eigenkmer2_2))>2:
            return False
        distance1=self.__get_distance(eigenkmer1_1,eigenkmer2_1)+self.__get_distance(eigenkmer1_2,eigenkmer2_2)
        distance2=self.__get_distance(eigenkmer1_1,eigenkmer2_2)+self.__get_distance(eigenkmer1_2,eigenkmer2_1)
        distance=min(distance1,distance2)
        if distance<=2:
            return True
        else:
            return False
        
    def __sort_family(self):
        self.head.sort(key = lambda x: x.score,reverse=True)  
        for head in self.head:
            head.family.sort(key=lambda x: x.get_aptamer_count(),reverse=True) 
    
    def __set_rank(self):
        tot=0
        for head in self.head:
            head.set_rank(tot)
            tot+=1
    
    def __get_distance(self, s1, s2):
        '''
        只允许开头结尾有不相同之处
        求两字符串编辑距离
        '''
        m = len(s1)
        n = len(s2)
        dp = [[0]*(n+1) for _ in range(2)]
        length = 0
        endi=0
        endj=0
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                try:
                    if s1[i - 1] == s2[j - 1]:
                        dp[i%2][j] = dp[(i - 1)%2][j - 1] + 1
                        if dp[i%2][j] > length:
                            length = dp[i%2][j]
                            endi = i
                            endj = j
                    else:
                        dp[i%2][j] = 0
                except IndexError:
                    print(m,n,i,j)
                    print(s1,s2)
                    sys.exit()
        # return the cost
        costFirst=max(endi,endj)-length
        costLast=max(m-endi,n-endj)
        return costFirst+costLast
    
    def __get_distance2(self, s1, s2):
        '''
        允许中间有不相同之处
        求两字符串编辑距离
        '''
        m = len(s1)
        n = len(s2)
        dp = [[0]*(n+1) for i in range(m+1)]
        for i in range(m+1):
            dp[i][0] = i
        for j in range(n+1):
            dp[0][j] = j
        for i in range(1,m+1):
            for j in range(1,n+1):
                if s1[i-1] == s2[j-1]:
                    dp[i][j] = min(dp[i-1][j-1],dp[i-1][j]+1,dp[i][j-1]+1)
                else:
                    dp[i][j] = min(dp[i-1][j],dp[i][j-1],dp[i-1][j-1])+1
        return dp[m][n]/max(m,n)

    def __is_valid(self,kmer,pos,kmer_len,kmer1,kmer2):
        position1 = [i for i, x in enumerate(kmer) if x == kmer1]
        position2 = [i for i, x in enumerate(kmer) if x == kmer2]
        for x in position1:
            for y in position2:
                if pos[x]<pos[y]:
                    if pos[y]-pos[x]>=kmer_len[x]:
                        return True
                else:
                    if pos[x]-pos[y]>=kmer_len[y]:
                        return True
        return False
    
    def __record_cluster_data(self):
        data=[]
        print(f'Number of Final cluster for round {self.round} is {len(self.head)}')
        cnt=0
        for x in range(len(self.head)):
            firstelement=self.head[x].family[0]
            aptamer=firstelement.get_aptamer_seq()
            count=firstelement.get_aptamer_count()
            data.append([cnt,self.head[x].eigenkmer1, self.head[x].eigenkmer2,self.head[x].size,self.head[x].score,aptamer,count])
            cnt+=1
            if cnt==100:
                break
        df=pd.DataFrame(data, columns=['index','kmer1','kmer2','size','score','representative','countforrepresentative'])
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Final_cluster_for_round{self.round}.csv', index=False)
