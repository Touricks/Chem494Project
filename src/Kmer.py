from Sequence import Sequence
from Aptamer import Aptamer
import pandas as pd
import os
class Kmer():
    def __init__(self,seqdata:list[Sequence],round,num1,num2):
        self.candidate={'5':{},'6':{},'7':{},'8':{},'9':{},'10':{},'11':{},'12':{}}
        self.cand_pair={'5':[],'6':[],'7':[]}
        self.data=seqdata
        self.round=round
        self.num1=num1
        self.num2=num2

    def __display_candidate(self):
        for kmerlen in range(5,8):
            print(f'kmerlen={kmerlen}')
            for seq in self.candidate[str(kmerlen)].keys():
                print(f'{seq} : {self.candidate[str(kmerlen)][seq]}')
            print('')

    def set_candidate(self):
        for data in self.data:
            x=data.apt
            y=data.count
            if y>1:
                for id in range(len(x.kmer)):
                    seq=x.kmer[id]
                    kmerlen=x.kmerlen[id]
                    if seq not in self.candidate[str(kmerlen)]:
                        self.candidate[str(kmerlen)][seq]=y
                    else:
                        self.candidate[str(kmerlen)][seq]+=y
        #elf.__display_candidate()
        #Display the chosen kmer to excel file
        data=[]
        for kmerlen in range(5,13):
            for seq in self.candidate[str(kmerlen)].keys():
                data.append([seq,kmerlen,self.candidate[str(kmerlen)][seq]])
        df=pd.DataFrame(data, columns=['kmer','length','count'])
        path= str(os.getcwd()) + "/testexcel"
        df.to_csv(path+f'/Featured Kmer for Round {self.round}.csv', index=False)

    def get_kmer_list(self):
        kmerlist=[]
        for kmerlen in range(12,4,-1):
            keys=self.candidate[str(kmerlen)].keys()
            candidate_list=[]
            for x in keys:
                candidate_list.append([x,self.candidate[str(kmerlen)][x]])
            candidate_list.sort(key=lambda x: x[1],reverse=True)
            
            if kmerlen<8:
                maxlen=min(self.num1,len(candidate_list))
                self.cand_pair[str(kmerlen)]=list(map(lambda x:x[0],candidate_list[:maxlen]))
            else:
                maxlen=min(self.num2,len(candidate_list))
                for x in candidate_list[:maxlen]:
                    kmerlist.append([x[0]])

        for kmer1 in self.cand_pair['6']:
            for kmer2 in self.cand_pair['6']:
                kmerlist.append([kmer1,kmer2])
        for kmer1 in self.cand_pair['5']:
            for kmer2 in self.cand_pair['7']:
                kmerlist.append([kmer1,kmer2])
        return kmerlist