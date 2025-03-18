
import sys
sys.path.insert(1, 'CODE AND EXPERIMENTS/CGR-pcmer/')
import VIRUSES

from VIRUSES.module import *
from VIRUSES.functions_CGR_PCMER import count_kmers

class PCmer:
  def __init__(self, n):
      self.seqKmers_PPAKWS = [0] * n

    
  def __set_seqKmers_PPAKWS(self, seq_kmers, weakStrong, aminoKeto, purPyr, index):
        self.seqKmers_PPAKWS.pop(index)
        self.seqKmers_PPAKWS.insert(index, [seq_kmers, [purPyr, aminoKeto, weakStrong]])
     
   
  def get_seqKmers_PPAKWS(self, i=None):
      if i == None: return self.seqKmers_PPAKWS # return all the list kmers_PPAKWS
      return self.seqKmers_PPAKWS[i]
  
  # CGR_i = 0.5(CGR_i-1 + CGR_si) with i = 1,...,ls and CGR_0 = (0.5, 0.5)
  def build_pcmer(self, seq, index):

    weakStrong = list(); aminoKeto = list(); purPyr = list()

    for i in range(len(seq)):

      #weak or strong
      if seq[i] == "C" or seq[i] == "G": # S
        weakStrong.append('S')

      elif seq[i] == "A" or seq[i] == "T": # W
          weakStrong.append('W')

      #amino or keto       
      if seq[i] == "C" or seq[i] == "A": # M
          aminoKeto.append('M')

      elif seq[i] == "G" or seq[i] == "T": # K
          aminoKeto.append('K')
      
      #purine or pyrimidine
      if seq[i] == "C" or seq[i] == "T": # Y
          purPyr.append('Y')
    
      elif seq[i] == "A" or seq[i] == "G": # R
          purPyr.append('R')

    self.__set_seqKmers_PPAKWS(seq, weakStrong, aminoKeto, purPyr, index)
   
    return self
  
  

