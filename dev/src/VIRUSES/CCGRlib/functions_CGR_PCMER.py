import sys
sys.path.insert(1, 'dev/src/')
import VIRUSES

from VIRUSES.CCGRlib.module import *
from VIRUSES.CCGRlib.cgr import CGR
import collections
from collections import OrderedDict

# Remove characters to the genome sequences improper (e.g., -, N, n, c, g, a, t)
# G C T/U A
def preprocessing_line(processing_line):
    processing_line = processing_line.replace('-', "").replace('N', "").replace('n', "").replace('c', "C").replace('g', "G").replace('a', "A").replace('t', "T")
    processing_line = processing_line.replace('U', "T").replace('u', "T")
    processing_line = processing_line.replace('B', "").replace('D', "").replace('E', "").replace('F', "").replace('H', "").replace('I', "").replace('L', "").replace('M', "").replace('N', "").replace('O', "").replace('P', "").replace('Q', "").replace('R', "").replace('S', "").replace('V', "").replace('Z', "").replace('J', "").replace('K', "").replace('W', "").replace('X', "").replace('Y', "")
    processing_line = processing_line.replace("'","").replace(' ', "").replace("[","").replace(']', "")
    return processing_line.replace('b', "").replace('d', "").replace('e', "").replace('f', "").replace('h', "").replace('i', "").replace('l', "").replace('m', "").replace('n', "").replace('o', "").replace('p', "").replace('q', "").replace('r', "").replace('s', "").replace('v', "").replace('z', "").replace('j', "").replace('k', "").replace('w', "").replace('x', "").replace('y', "")


# The fasta file is trasformed in genome sequence
def parse_sequence(fastafile):
  with open(fastafile) as linefile:
    next(linefile) # skip the header line in the fasta file
    lines_genome = list(csv.reader(linefile))
  seq = ''
  for line in lines_genome[:-1]:
    processed_line = preprocessing_line(line[0])
    seq = seq + processed_line
  return seq

# Generate a collection composed by {"k-mers names : number of occurance of the k-mers"} e.g., 'GAT' : 114}
def count_kmers(sequence, k):
  sequence = preprocessing_line(sequence)
  kmer_count = collections.defaultdict(int)
  for i in range(len(sequence)-(k-1)):
    key = sequence[i:i+k]
    kmer_count[key] += 1
  return list(kmer_count.keys()), list(kmer_count.values())

# count_kmers with jellyfish
def count_kmers_jellyfish(fasta, k):
  my_dir = 'dev/lib/jellyfish/'
  cmd = my_dir + 'jellyfish-binary'

  # Ensure the binary is executable
  if not os.access(cmd, os.X_OK):
      print(f"Error: {cmd} not executable")
      return

    # Run the jellyfish count command
  try:
      subprocess.run([cmd, 'count', '-m', str(k), '-s', '100M', '-t', '10', '-o', my_dir + 'mer_counts.jf', fasta], check=True)
      #subprocess.run([cmd + ' dump ' + my_dir + 'mer_counts.jf -c', '> ' + '], check=True) 
      with open(my_dir + 'mer_counts_dumps.csv', 'w') as output_file:
        # Run the jellyfish dump command and write the output to the file
        subprocess.run([cmd, 'dump', my_dir + 'mer_counts.jf', '-c'], stdout=output_file, check=True)
  except subprocess.CalledProcessError as e:
      print(f"Error during command execution: {e}")
  except PermissionError as e:
      print(f"Permission error: {e}")


  count_kmers = pd.read_csv(my_dir + 'mer_counts_dumps.csv', delimiter = ' ', skiprows=1, header=None, low_memory=False)
  kmers = count_kmers.iloc[:, 0].values
  frequency = count_kmers.iloc[:, 1].values
  os.remove(my_dir + 'mer_counts.jf')
  os.remove(my_dir + 'mer_counts_dumps.csv')

  return kmers, frequency

# ADDITIONAL EXPERIMENTS
# additional method for codifyRGB_MaxEncodingFrequency: return for the list 'arr' the maximum configuration for 'arr'
def sorted_inMaxRepresentation(arr, len_vector):
  for idx in range(0, len_vector):
    tmp = arr[idx]
    max_value = arr[idx]; max_idx = idx
    for idx_n in range(idx+1, len_vector):
      if max_value < arr[idx_n]: max_value = arr[idx_n]; max_idx = idx_n
    if max_value != arr[idx]:
      arr.pop(idx); arr.insert(idx, max_value); arr.pop(max_idx); arr.insert(max_idx, tmp)
  return arr

# The functions realize the map {pcmer rapresentation:color} with a hash table H
def encodingColour(k, items, threshold): 
  dict_encoding_rgb = {}
  codeRGB = 0; range = (256 * threshold)/2**k
  for item in product(items, repeat=k):
    dict_encoding_rgb.update({item : codeRGB})
    codeRGB += range


  return dict_encoding_rgb

# ADDITIONAL EXPERIMENTS
# return for the kmer codify in chemical/physical structure l the presence and the char count e.g., RRRYY -> 3,2,R,Y
def countPresence(l):
  char_one = ""; char_two = ""
  for i in l:
    if i != char_one and char_one == "": char_one = i
    elif i != char_one and char_one != "": char_two = i
  count_one = 0; count_two = 0    
  for j in l:
    if j == char_one: count_one += 1
    else: count_two += 1     
  return count_one, count_two, char_one, char_two 

# ADDITIONAL EXPERIMENTS
# support method of encodingColourNORIPETITIONS; return if a item is already present in the presents' dictionary
def checkNORIPETIONS(item, dict_encoding_rgb):
  count_oneItem, count_twoItem, char_oneItem, char_twoItem = countPresence(item) # e.g., 3,2,R,Y
  char_count_item = [str(count_oneItem) + char_oneItem, str(count_twoItem) + char_twoItem] # e.g., 3R, 2Y
  for encoding_rgb in dict_encoding_rgb:

    # check if 3R, 2Y already exist in encoding_rgb
    if char_count_item[0] == encoding_rgb[0] and char_count_item[1] == encoding_rgb[1]: return True, None
    elif char_count_item[0] == encoding_rgb[1] and char_count_item[1] == encoding_rgb[0]: return True, None
  
  # not exist a ripetitions
  return False, char_count_item


# ADDITIONAL EXPERIMENTS
# The functions realize the map {pcmer rapresentation:color} with a dictionary without ripetitions (it count only the presence)
def encodingColourNORIPETITIONS(k, Npresences, items, threshold): 
  dict_encoding_rgb = {}
  codeRGB = 0; range = (256 * threshold)/Npresences # a range between [0..255] that contains already Npresences
  for item in product(items, repeat=k):
    ripetitions, char_count_item = checkNORIPETIONS(item, dict_encoding_rgb) # check if item is already presents in dict_encoding_rgb
    if ripetitions == False: # if it not present then is update the dictionary with the correspondent codeRGB
      dict_encoding_rgb.update({tuple(char_count_item) : codeRGB}) # e.g., [2Y, 1M]
      codeRGB += range

  return dict_encoding_rgb

# Calculate the ratio between max range for frequency and the max frequency of the kmers in fcgr
def ratioFreq(maxFreq, channel_freq):
  return channel_freq / maxFreq

# Return the value of kmer pcmer that occurance of more
def max_freqPCMER(dict_d):
  max = 0
  for pcmer_count in list(dict_d.values()):
    if max < pcmer_count:
      max = pcmer_count

  return max

# Return the frequency in the fcgr for the kmer
def findFreqKmer(kmer, kmer_freq):
  for item in kmer_freq:
    if item[0] == kmer:
      return item[1]
    
  else: return -1







