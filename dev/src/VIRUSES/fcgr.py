import sys
sys.path.insert(1, 'dev/src/')
from module import *
import VIRUSES

from VIRUSES.cgr import CGR
from VIRUSES.pcmer import PCmer
from VIRUSES.functions_CGR_PCMER import count_kmers

## fcgr algorithms

class FCGR:
  def __init__(self, seq = "", kmer = 0, kmer_freq = {}):
    self.seq = seq
    self.kmer = kmer
    self.kmer_freq = kmer_freq

  def get_kmer(self):
    return self.kmer
    
  def get_seq(self):
    return self.seq
    
  def get_kmer_freq(self):
    return self.kmer_freq

  def __set_kmer(self, k):
    self.kmer = k

  def __set_seq(self, seq):
    self.seq = seq


  def __set_kmer_freq(self, index, seq_kmer, freq):
    if len(self.get_kmer_freq()) >= index + 1:
      self.kmer_freq.pop(index)
    self.kmer_freq.insert(index, [seq_kmer, freq])


  def init_seq_kmer(self, array_size):
    for i in range(0, array_size * array_size):
      self.__set_kmer_freq(i, 0, 0)
      

  #def __set_kmer_freq(self, matrix_fcgr, seq_fcgr):
  #  list_freq = list()
  #  for j in matrix_fcgr:
  #    for j_i in j:
  #      freq = j_i
  #      list_freq.append(freq)
    
  #  list_seq = list()
  #  for j in seq_fcgr:
  #   for j_i in j:
  #      seq = j_i
  #     list_seq.append(seq)
  #  for i in range(0, len(list_freq)):
  #    self.kmer_freq.append([list_seq[i], list_freq[i]])

  def tot_freq(self):
    totFreq = 0
    for kmer_freq in self.kmer_freq:
      totFreq += kmer_freq[1]

    return totFreq
  
  def max_freq(self):
    max = 0
    for kmer_freq in self.kmer_freq:
      if max < kmer_freq[1]:
        max = kmer_freq[1]
        
    return max
    

  # Design the FCGR in file .png
  def __design_fcgr(self, title, directory, matrix_fcgr, seq_fcgr):
    fig, (ax) = plt.subplots(1, 1, figsize=(20, 18))
    k = self.kmer
    im = ax.imshow(matrix_fcgr, interpolation='nearest', cmap=cm.gray_r)
    # Loop over data dimensions and create text annotations.
    #for i in range(int(math.sqrt(2**(2*k)))):
      #for j in range(int(math.sqrt(2**(2*k)))):
          #if k <= 2:
            #text = ax.text(j, i, str(seq_fcgr[i][j]) + '(' + str(matrix_fcgr[i][j]) + ')', fontsize = 25,
                           #ha="center", va="center", color="c", fontweight = 'bold')
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    #ax.set_title(title + ', k = ' + str(k), fontsize=35)
    plt.savefig(directory + '/' + 'FCGR' + title + "(k = "+ str(k) + ")")
    plt.clf(); plt.close()

    # save in GS format
    img = Image.open(directory + '/' + 'FCGR' + title + "(k = "+ str(k) + ").png").convert('L')
    img.save(directory + '/' + 'FCGR' + title + "(k = "+ str(k) + ").png")

  def build_fcgr(self, seq, k, title = None, directory = None, flag_design = False):
    self.__set_seq(seq)
    self.__set_kmer(k)
    kmers, frequencies = count_kmers(seq, k)
    array_size = int(math.sqrt(2**(2*k))) # cells for row/colomn present in the fcgr bi-dimensional matrix
    self.init_seq_kmer(array_size)

    matrix_fcgr = []; seq_fcgr = []
    for i in range(array_size):
      matrix_fcgr.append([0]*array_size)
      seq_fcgr.append([0]*array_size)
    for n in range(0, len(kmers)):
      seq_kmer = kmers[n]; freq = frequencies[n]
      CGR_tmp = CGR('N', list()).build_cgr(seq_kmer)
      cgr_seq = CGR_tmp.get_cgr_seq()
      x = cgr_seq[k][0]; y = cgr_seq[k][1] # coordinates (x,y) in fcgr for the k-mers
      
      x_supp = -1; y_supp = -1
      i = -1; j = array_size
      ratio_cell = 2 / array_size

      # calculate the exact x position in fcgr matrix for the k-mers (sequence_data)
      while x_supp <= x:
        x_supp += ratio_cell
        i += 1 # 0..array_size-1

      # calculate the exact y postion in fcgr matrix for the k-mers (sequence_data)
      while y_supp <= y:
        y_supp += ratio_cell
        j -= 1 # array_size-1..0 (i.e., 'AA': y_coord = 3 and y = 1; 0.5+0.5+0.5+0.5, y=4-1-1-1=3)


      matrix_fcgr[j][i] = freq # [j][i] coord y, x = riga j, colonna i
      seq_fcgr[j][i] = seq_kmer

      index = j * array_size + i
      self.__set_kmer_freq(index, seq_kmer, freq)

    #self.__set_kmer_freq(matrix_fcgr, seq_fcgr)

    if flag_design == True:
      self.__design_fcgr(title, directory, matrix_fcgr, seq_fcgr)

    return self
  

  def set_fcgr_intofcgr_pcmer_rgb(self, seq, k, index, seq_kmer, freq):
    self.__set_seq(seq)
    self.__set_kmer(k)
    self.__set_kmer_freq(index, seq_kmer, freq)
  



# REWRITE CODE

class FCGR_PCMER_RGB2(FCGR, PCmer):
  def __init__(self, seq = "", kmer = 0, kmer_freq = list(), seqKmers_PPAKWS_rgb = list()):
    FCGR.__init__(self, seq = "", kmer = 0, kmer_freq = list())
    PCmer.__init__(self, seq_kmers = list(), kmer = 0, seqKmers_PPAKWS_rgb = list())

  def __design_fcgr(self, title, directory):
    len_square = int(math.sqrt(2**(2*self.kmer)))
    new_list = []; list_rgb = []
    count = 0; pc = 0
    freq = self.kmer_freq
    pcmer = self.seqKmers_PPAKWS_rgb
    for i in freq:
        count += 1
        if i == [0, 0]: list_rgb.append([0, 0, 0])
        else: list_rgb.append(pcmer[pc][2]); pc += 1

        if count == len_square: new_list.append(list_rgb); count = 0; list_rgb = []

    grid = np.array(new_list, dtype=np.uint8)

    # save in RGB format
    im = Image.fromarray(grid).resize((1600,1600), resample=Image.NEAREST).convert('RGB')

    # displaying the title
    #plt.title(title + ', k = ' + str(self.kmer), fontsize=10)
    plt.axis('off')
    plt.imshow(im)
    plt.savefig(directory + '/' + 'FCGR-RGB' + title + "(k = "+ str(self.kmer) + ")")
    plt.clf(); plt.close()

    im = Image.open(directory + '/' + 'FCGR-RGB' + title + "(k = "+ str(self.kmer) + ")" + ".png")

    background = Image.new("RGB", im.size, (255, 255, 255))
    background.paste(im, mask = im.split()[3])
    background.save(directory + '/' + 'FCGR-RGB' + title + "(k = "+ str(self.kmer) + ").png", "PNG", quality=100)
  
  
    return self


  def build_fcgr(self, seq, k, method_rgb = 'binarization', title = None, directory = None, flag_design_gs = False):
    super().build_fcgr(seq, k, title, directory, flag_design_gs)
    super().build_pcmer(seq, k, method_rgb)
    self.__design_fcgr(title, directory)
    return self
  




## REWRITE CODE
class FCGR_RGB(FCGR):
  def __init__(self, seq = "", kmer = 0, kmer_freq = list()):
    FCGR.__init__(self, seq = "", kmer = 0, kmer_freq = list())
    self.rgb = {}

  def __set_rgb(self, dict_rgb):
    self.rgb = dict_rgb

  def __set_kmer(self, k):
    self.kmer = k

  def __set_seq(self, seq):
    self.seq = seq

  def __set_kmer_freq(self, matrix_fcgr, seq_fcgr):
    list_freq = list()
    for j in matrix_fcgr:
      for j_i in j:
        freq = j_i
        list_freq.append(freq)
    
    list_seq = list()
    for j in seq_fcgr:
      for j_i in j:
        seq = j_i
        list_seq.append(seq)
    for i in range(0, len(list_freq)):
      self.kmer_freq.append([list_seq[i], list_freq[i]])
    

  def __design_fcgr_classical(self, title, directory):
    len_square = int(math.sqrt(2**(2*self.kmer)))
    new_list = []; list_rgb = []
    count = 0
    for color in self.rgb.values():
      count += 1
      list_rgb.append(color)
      if count == len_square: new_list.append(list_rgb); count = 0; list_rgb = []
    grid = np.array(new_list, dtype=np.uint8)

    # save in RGB format
    im = Image.fromarray(grid).resize((1600,1600), resample=Image.NEAREST).convert('RGB')

    # displaying the title
    #plt.title(title + ', k = ' + str(self.kmer), fontsize=10)
    plt.axis('off')
    plt.imshow(im)
    plt.savefig(directory + '/' + 'FCGR-RGBClassical' + title + "(k = "+ str(self.kmer) + ")")
    plt.clf(); plt.close()

    im = Image.open(directory + '/' + 'FCGR-RGBClassical' + title + "(k = "+ str(self.kmer) + ")" + ".png")
    background = Image.new("RGB", im.size, (255, 255, 255))
    background.paste(im, mask = im.split()[3])
    background.save(directory + '/' + 'FCGR-RGBClassical' + title + "(k = "+ str(self.kmer) + ").png", "PNG", quality=100)
  
    return self
  
  # brg map: blue is assigned to the kmers less frequencies, green to the more frequencies and red to the middle frequencies (scaling)
  def __design_fcgr_mapBRG(self, title, directory, matrix_fcgr, seq_fcgr):
    fig, (ax) = plt.subplots(1, 1, figsize=(20, 18))
    k = self.kmer
    im = ax.imshow(matrix_fcgr, interpolation='nearest', cmap=cm.brg)
    # Loop over data dimensions and create text annotations.
    #for i in range(int(math.sqrt(2**(2*k)))):
      #for j in range(int(math.sqrt(2**(2*k)))):
          #if k <= 2:
            #text = ax.text(j, i, str(seq_fcgr[i][j]) + '(' + str(matrix_fcgr[i][j]) + ')', fontsize = 25,
                           #ha="center", va="center", color="c", fontweight = 'bold')
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    #ax.set_title(title + ', k = ' + str(k), fontsize=35)
    plt.savefig(directory + '/' + 'FCGR-BRG' + title + "(k = "+ str(k) + ")")
    plt.clf(); plt.close()

    # save in RGBA format
    img = Image.open(directory + '/' + 'FCGR-BRG' + title + "(k = "+ str(k) + ").png").convert('RGB')

    im = Image.open(directory + '/' + 'FCGR-BRG' + title + "(k = "+ str(k) + ").png")
  

    background = Image.new("RGB", im.size, (255, 255, 255))
    background.paste(im, mask = im.split()[3])
    background.save(directory + '/' + 'FCGR-BRG' + title + "(k = "+ str(k) + ").png", "PNG", quality=100)

 # In the building FCGR is count the points in CGR with the information frequency
  def build_fcgr(self, seq, k, title = None, directory = None, flag_design = False):
    self.__set_seq(seq)
    self.__set_kmer(k)
    frequency = count_kmers(seq, k)
    array_size = int(math.sqrt(2**(2*k))) # cells for row/colomn present in the fcgr bi-dimensional matrix
    matrix_fcgr = []; seq_fcgr = []
    for i in range(array_size):
      matrix_fcgr.append([0]*array_size)
      seq_fcgr.append([0]*array_size)
    for seq_kmer, freq in frequency.items():
      CGR_tmp = CGR('N', list()).build_cgr(seq_kmer)
      cgr_seq = CGR_tmp.get_cgr_seq()
      x = cgr_seq[k][0]; y = cgr_seq[k][1] # coordinates (x,y) in fcgr for the k-mers
      x_supp = -1; y_supp = -1
      i = -1; j = array_size
      ratio_cell = 2 / array_size
      # calculate the exact x position in fcgr matrix for the k-mers (sequence_data)
      while x_supp <= x:
        x_supp += ratio_cell
        i += 1 # 0..array_size-1
      # calculate the exact y postion in fcgr matrix for the k-mers (sequence_data)
      while y_supp <= y:
        y_supp += ratio_cell
        j -= 1 # array_size-1..0 (i.e., 'AA': y_coord = 3 and y = 1; 0.5+0.5+0.5+0.5, y=4-1-1-1=3)
      matrix_fcgr[j][i] = freq # [j][i] coord y, x = riga j, colonna i
      seq_fcgr[j][i] = seq_kmer


      index = j * array_size + i
      self.__set_kmer_freq(index)
    
    #self.__set_kmer_freq(matrix_fcgr, seq_fcgr)

    if flag_design == True:
      self.__design_fcgr_mapBRG(title, directory, matrix_fcgr, seq_fcgr)

    return self
  
  # This method using Schema Colours and no doesn't distinguish between kmer 'AT' and 'TA'
  def __classicalNOFrequcencyRGB(kmer_freq):
    sumMaximumFrequency = 0
    dict_rgb = {}
    for i in kmer_freq:
      freq = i[1]
      sumMaximumFrequency += freq
  
    for i in kmer_freq:
      kmer = i[0]; freq = i[1]; rgb_Red = 0; rgb_Green = 0; rgb_Blue = 0
      if kmer == 0:
        rgb_Red = 0; rgb_Green = 0; rgb_Blue = 0
        rgb_value = [rgb_Red, rgb_Green, rgb_Blue]
        dict_rgb.update({kmer:(rgb_value)})
      else:
        for base in kmer: ## DNA color scheme (A azure, T tweety bird, G green, C carmine)
          if base == 'A': rgb_Red += 0; rgb_Green += 127; rgb_Blue += 255
          elif base == 'T': rgb_Red += 255; rgb_Green += 255; rgb_Blue += 0
          elif base == 'C': rgb_Red += 0; rgb_Green += 128; rgb_Blue += 0
          elif base == 'G': rgb_Red += 150; rgb_Green += 0; rgb_Blue += 24
    
        # the maximum value of the channel is 255
        if rgb_Red > 255: rgb_Red = 255
        if rgb_Green > 255: rgb_Green = 255
        if rgb_Blue > 255: rgb_Blue = 255
        rgb_value = [rgb_Red, rgb_Green, rgb_Blue]
        dict_rgb.update({kmer:(rgb_value)})
  
    return dict_rgb



  def build_fcgr_NOFREQUENCY(self, seq, k, title = None, directory = None, flag_design_gs = False):
    super().build_fcgr(seq, k, title, directory, flag_design_gs)
    self.__set_rgb(self.__classicalNOFrequcencyRGB(self.kmer_freq))
    self.__design_fcgr_classical(title, directory)
    return self


  

  
    
