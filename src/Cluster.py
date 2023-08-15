from Sequence import Sequence
import math
import pandas as pd
import os
import sys

class Member():
    def __init__(self,seq:Sequence,kmer1:str,kmer2:str,weighted:bool):
        '''
        self.size: Number of aptamers in the cluster
        '''
        self.eigenkmer1=kmer1
        self.eigenkmer2=kmer2
        self.weighted=weighted
        self.score=0
        self.rank=0
        self.merged=False
        self.family=[seq]
        if self.weighted:
            self.size=seq.get_aptamer_count()
        else:
            self.size=1

    def calculate_score(self):
        size=self.get_size()
        if self.eigenkmer2=='':
            kmer_score=4**(len(self.eigenkmer1)-8)
        else:
            kmer_score=16
        self.score=(size-1)*kmer_score

    def equal(self,kmer1,kmer2):
        if (self.eigenkmer1==kmer1)and(self.eigenkmer2==kmer2):
            return True
        else:
            return False
        
    def add_family(self,seq:Sequence):
        self.family.append(seq)
        if self.weighted:
            self.size+=seq.get_aptamer_count()
        else:
            self.size+=1

    def get_size(self):
        return self.size
    
    def get_score(self):
        return self.score
    
    def set_rank(self,rank):
        self.rank=rank
        for seq in self.family:
            seq.set_cluster_num(rank)

class Cluster():
    def __init__(self,seqdata:list[Sequence],weighted:bool,kmerlist:list[str],round,goal):
        self.head=[]
        self.data=seqdata
        self.weighted=weighted
        self.kmers=kmerlist
        self.size=0 #number of clusters
        self.testdata=[]
        self.round=round
        self.goal=goal
    
    def set_cluster_data(self):
        self.__get_cluster()
        self.__sort_cluster()
        self.__combine_cluster()
        self.__sort_family()
        self.__set_rank()    
        self.__record_cluster_data()

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
            newhead=Member(seq,eigenkmer,'',self.weighted)
            self.head.append(newhead)
            self.size+=1
        else:
            member=self.head[-1]
            if member.equal(eigenkmer,''):
                member.add_family(seq)
            else:
                newhead=Member(seq,eigenkmer,'',self.weighted)
                self.head.append(newhead)
                self.size+=1
                
    def __add_cluster_shortkmer(self,seq:Sequence,x1:str,x2:str):
        if (self.head==[]):
            newhead=Member(seq,x1,x2,self.weighted)
            self.head.append(newhead)
            self.size+=1
        else:
            member=self.head[-1]
            if member.equal(x1,x2):
                member.add_family(seq)
            else:
                newhead=Member(seq,x1,x2,self.weighted)
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
        if distance/length<=0.2:
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
        for x in range(len(self.head)):
            data.append([self.head[x].eigenkmer1, self.head[x].eigenkmer2,self.head[x].size,self.head[x].score])
        df=pd.DataFrame(data, columns=['kmer1','kmer2','size','score'])
        path= str(os.getcwd()) + "/testexcel"
        df.to_csv(path+f'/Final_cluster_for_round{self.round}.csv', index=False)