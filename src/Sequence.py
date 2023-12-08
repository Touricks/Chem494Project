from Aptamer import Aptamer
class Sequence():
    def __init__(self,apt:Aptamer,count):
        self.apt = apt
        self.count = count
        self.used=False
        self.clunum=-1

    def set_cluster_num(self,num):
        self.clunum=num

    def get_cluster_num(self):
        return self.clunum
    
    def set_state(self,state):
        self.used=state

    def get_state(self):
        return self.used

    def get_aptamer_seq(self):
        x=self.apt
        return x.get_sequence()
    
    def get_aptamer_count(self):
        return self.count
    
    def get_aptamer_dG(self):
        x=self.apt
        return x.get_dg()
    
    def get_aptamer_ct(self):
        x=self.apt
        return x.get_ct()
    
    def get_kmer_list(self):
        x=self.apt
        kmers,pos,length=x.get_kmer()
        return kmers
    
    def get_kmer_list2(self):
        x=self.apt
        kmers=x.get_kmer2()
        return kmers