class Aptamer():
    def __init__(self,seq,ct,dg):
        self.seq=seq
        self.ct=ct
        self.dg=dg
        self.origin_kmer=[]
        self.kmer=[]
        self.kmerpos=[]
        self.kmerlen=[]
        self.hairpinnum=0
        self.hairpinpos=[]
        self._add_hairpin()
        if self.hairpinnum!=0:
            self._add_kmer()

    def _add_hairpin(self):
        s=self.ct
        start=0
        endpos=0
        i=0
        startstack=[]
        #print(s)
        while i<len(s):
            if (s[i]=='(')and(s[i+1]!='('):
                startstack.append(i)
            elif (s[i]==')'):
                cnt=1
                endpos=i
                #print(startstack)
                start=startstack.pop()
                while (start>0)and(endpos<len(s)-1)and(s[start-1]=='(')and(s[endpos+1]==')'):
                    start-=1
                    endpos+=1 
                    cnt+=1
                if cnt>=4:
                    self.hairpinnum+=2
                    self.hairpinpos.append([start,start+cnt-1])
                    self.hairpinpos.append([endpos-cnt+1,endpos])
                if (start>0)and(s[start-1]=='('):
                    startstack.append(start-1)
                i=endpos
            i+=1

    def _add_kmer(self):
        for x in range(self.hairpinnum+1):
            startpos=0
            if x==0:
                seq=self.seq[0:self.hairpinpos[x][0]]
            elif x==self.hairpinnum:
                seq=self.seq[self.hairpinpos[x-1][1]+1:]
                startpos=self.hairpinpos[x-1][1]+1
            else:
                seq=self.seq[self.hairpinpos[x-1][1]+1:self.hairpinpos[x][0]]
                startpos=self.hairpinpos[x-1][1]+1
            if seq!='':
                self.origin_kmer.append(seq)
            for kmerlen in range(12,4,-1):
                for x in range(len(seq)-kmerlen+1):
                    self.kmer.append(seq[x:x+kmerlen])
                    self.kmerpos.append(x+startpos)
                    self.kmerlen.append(kmerlen)

    def get_kmer(self):
        return self.kmer,self.kmerpos,self.kmerlen
    
    def get_hairpin_number(self):
        return self.hairpinnum
    
    def get_sequence(self):
        return self.seq
    
    def get_dg(self):
        return self.dg
    
    def get_ct(self):
        return self.ct

    def get_kmer2(self):
        return self.origin_kmer
