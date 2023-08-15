from Sequence import Sequence
from Aptamer import Aptamer
import pandas as pd
import os
class Kmer():
    def __init__(self,seqdata:list[Sequence],weighted,round):
        self.candidate={'5':{},'6':{},'7':{},'8':{},'9':{},'10':{},'11':{},'12':{}}
        self.cand_pair={'5':[],'6':[],'7':[]}
        self.data=seqdata
        self.weighted=weighted
        self.round=round

    def set_candidate(self):
        for data in self.data:
            x=data.apt
            y=data.count
            if y>1:
                for id in range(len(x.kmer)):
                    seq=x.kmer[id]
                    kmerlen=x.kmerlen[id]
                    if self.weighted:
                        if seq not in self.candidate[str(kmerlen)]:
                            self.candidate[str(kmerlen)][seq]=y
                        else:
                            self.candidate[str(kmerlen)][seq]+=y
                    else:
                        if seq not in self.candidate[str(kmerlen)]:
                            self.candidate[str(kmerlen)][seq]=1
                        else:
                            self.candidate[str(kmerlen)][seq]+=1

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
            key=self.candidate[str(kmerlen)].keys()
            candidate_list=[]
            for x in key:
                candidate_list.append([x,self.candidate[str(kmerlen)][x]])
            candidate_list.sort(key=lambda x: x[1],reverse=True)
            if kmerlen<8:
                maxlen=min(30,len(candidate_list))
            else:
                maxlen=min(500,len(candidate_list))
            if kmerlen<8:
                self.cand_pair[str(kmerlen)]=list(map(lambda x:x[0],candidate_list[:maxlen]))
            else:
                for x in candidate_list[:maxlen]:
                    kmerlist.append([x[0]])
        for kmer1 in self.cand_pair['6']:
            for kmer2 in self.cand_pair['6']:
                kmerlist.append([kmer1,kmer2])
        for kmer1 in self.cand_pair['5']:
            for kmer2 in self.cand_pair['7']:
                kmerlist.append([kmer1,kmer2])
        return kmerlist