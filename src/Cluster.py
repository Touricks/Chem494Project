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
        f1=[13,12,11,10,9,45] #N=20
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
            seq.set_cluster_num(rank)

class Cluster():
    def __init__(self,seqdata:list[Sequence],kmerlist,round,goal,N):
        self.head=[]
        self.head2={}
        self.data=seqdata
        self.kmers=kmerlist
        self.size=0 #number of clusters
        self.testdata=[]
        self.round=round
        self.goal=goal
        self.N=N
        self.count_for_test=0
    
    def set_cluster_data(self):
        self.__get_cluster()
        self.__sort_cluster()
        self.__combine_cluster()
        self.__sort_family()
        self.__set_rank()    
        self.__record_cluster_data()

    def set_cluster_data2(self):
        for seq in self.data:
            feature=tuple(seq.get_kmer_list2())
            if feature in self.head2.keys():
                self.head2[feature]['count']+=seq.get_aptamer_count()
                self.count_for_test+=1
            else:
                self.head2[feature]={}
                self.head2[feature]['represent']=seq.get_aptamer_seq()
                self.head2[feature]['mfestructure']=seq.get_aptamer_ct()
                self.head2[feature]['count']=seq.get_aptamer_count()
                self.head2[feature]['flag']=0
        '''
        cnt=0
        for x in list(self.head2.keys()):
            cnt+=1
            print(cnt)
            if self.head2[x].get('flag')!=0:
                continue
            else:
                self.head2[x]['flag']=1
                for y in list(self.head2.keys()):
                    if self.head2[y].get('flag')!=0:
                        continue
                    else:
                        if self.__similarity2(self.head2[x]['represent'],self.head2[y]['represent']):
                            self.head2[x]['count']+=self.head2[y]['count']
                            self.head2[y]['flag']=2
        self.head2 = {key: value for key, value in self.head2.items() if value['flag'] < 2}
        '''
        for x in self.head2:
            self.head2[x]['relcount']=self.head2[x]['count']/self.N

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
        data=[]
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
        '''
        cnt=1
        for head in self.head:
            data.append([cnt,head.eigenkmer1,head.eigenkmer2,head.get_size(),head.get_score()])
            cnt+=1
            if cnt==100:
                break
        df=pd.DataFrame(data, columns=['cluster','kmer1','kmer2','size','score'])
        path= str(os.getcwd()) + "/dataset"
        df.to_csv(path+f'/Cluster_data_for_round{self.round}.csv', index=False)
        '''

    def __levenshtein_distance(self,str1,str2):
        if len(str1) < len(str2):
            return self.__levenshtein_distance(str2,str1)

        if len(str2) == 0:
            return len(str1)

        previous_row = range(len(str2) + 1)
        for i, c1 in enumerate(str1):
            current_row = [i + 1]
            for j, c2 in enumerate(str2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        return previous_row[-1]
    
    def __similarity2(self,str1:str,str2:str):
        dis=0
        for x in range(len(str1)):
            if str1[x]!=str2[x]:
                dis+=1
        if dis/max(len(str1),len(str2))<=0.06:
            return True
        else:
            return False
        
    def __get_cluster(self):
        flag=1
        for x in range(len(self.kmers)):
            if len(self.kmers[x])==1:
                self.__set_longkmer(self.kmers[x])
            elif flag:
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
        length=max(len(eigenkmer1_1)+len(eigenkmer1_2),len(eigenkmer2_1)+len(eigenkmer2_2))
        if distance/length<=0.1:
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