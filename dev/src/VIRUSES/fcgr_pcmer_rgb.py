import sys
sys.path.insert(1, 'dev/src/')
import VIRUSES

from VIRUSES.module import *
from VIRUSES.pcmer import PCmer
from VIRUSES.fcgr import FCGR
from VIRUSES.cgr import CGR
from VIRUSES.functions_CGR_PCMER import encodingColour, ratioFreq, count_kmers, count_kmers_jellyfish

class FCGR_PCMER_RGB(FCGR, PCmer):
    def __init__(self, n, seq = "", kmer = 0, kmer_freq = list(), kmers_rgb = {}):
        FCGR.__init__(self, seq = "", kmer = 0, kmer_freq = list())
        PCmer.__init__(self, n)
        self.kmers_rgb = kmers_rgb

    def __set_kmers_rgb(self, threshold, list_kmer_rgb):
        self.kmers_rgb.update({threshold : list_kmer_rgb})


    def get_kmers_rgb(self):
        return self.kmers_rgb
    
#### kCCGR ##################################################################################################################################
       
## maxFreq is calculate in the maximum frequency existing in the FCGR

    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR(self, maxFreq, thresholds):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency


        for threshold in thresholds:
            list_kmer_rgb = list()
            if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
            elif threshold == 0: channel_freq = 256
            else: channel_freq = 0
            ratio_freq = ratioFreq(maxFreq, channel_freq) 
            index = 0



            if threshold > 0: # structural encoding
                threshold_k = int(math.log2(256 * threshold))
                items_pp = ['R', 'Y']; items_ak = ['M', 'K']; items_ws = ['W', 'S']
                struct_codeRGB_pp_threshold_k = encodingColour(threshold_k, items_pp, threshold); struct_codeRGB_pp_mod = encodingColour(self.kmer%threshold_k, items_pp, threshold)
                struct_codeRGB_ak_threshold_k = encodingColour(threshold_k, items_ak, threshold); struct_codeRGB_ak_mod = encodingColour(self.kmer%threshold_k, items_ak, threshold)
                struct_codeRGB_ws_threshold_k = encodingColour(threshold_k, items_ws, threshold); struct_codeRGB_ws_mod = encodingColour(self.kmer%threshold_k, items_ws, threshold)

                count_div = math.ceil(self.kmer/threshold_k)

            else: threshold_k = 0 # only frequency


            for entry in self.get_seqKmers_PPAKWS():

                if entry == 0: index += 1; continue
                kmer = entry[0]
                pp = entry[1][0]
                ak = entry[1][1]
                ws = entry[1][2]
                freq = (self.get_kmer_freq())[index][1]

                inc_freq = freq * ratio_freq

                # 0 in frequency is reserved for empy k-mer, 255 is reserved for the two k-mer more frequency
                # in structural encoding are distribuited [0,255] value

                if inc_freq > 255: inc_freq = 255

    
                k = self.kmer
                pp_rest = pp; ak_rest = ak; ws_rest = ws
                start = 0; rgb_pp = 0; rgb_ak = 0; rgb_ws = 0
                if threshold_k > 0:
                    while k > 0: # until of k > 0 devide kmer in number of pieces
                        if k >= threshold_k:
                            struct_codeRGB_pp = struct_codeRGB_pp_threshold_k
                            struct_codeRGB_ak = struct_codeRGB_ak_threshold_k
                            struct_codeRGB_ws = struct_codeRGB_ws_threshold_k
                        else:
                            struct_codeRGB_pp = struct_codeRGB_pp_mod
                            struct_codeRGB_ak = struct_codeRGB_ak_mod
                            struct_codeRGB_ws = struct_codeRGB_ws_mod



                        if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                        pp_aux = pp_rest[start:start+threshold_k]; ak_aux = ak_rest[start:start+threshold_k]; ws_aux = ws_rest[start:start+threshold_k]
                        pp_rest = pp_rest[start+threshold_k:]; ak_rest = ak_rest[start+threshold_k:]; ws_rest = ws_rest[start+threshold_k:]
                        k -= threshold_k

                        rgb_pp += struct_codeRGB_pp[tuple(pp_aux)] # code rgb for pp_aux
                        rgb_ak += struct_codeRGB_ak[tuple(ak_aux)]
                        rgb_ws += struct_codeRGB_ws[tuple(ws_aux)]
   
                    rgb_pp /= count_div # devide rgb in number of pieces
                    rgb_ak /= count_div
                    rgb_ws /= count_div
                
    
                    rgb_pp += inc_freq # sum between frequency and pp structure
                    rgb_ak += inc_freq
                    rgb_ws += inc_freq



                else: ## threshold_k == 0 and theshold == 0 (only frequence)
                    rgb_pp += inc_freq; rgb_ak += inc_freq; rgb_ws += inc_freq

                

                list_kmer_rgb.append([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])
                index += 1

                # reset threshold_k for next k-mer in PPAKWS
                if threshold > 0: threshold_k = int(math.log2(256 * threshold))
                else: threshold_k = 0


            self.__set_kmers_rgb(threshold, list_kmer_rgb)

#### pcCCGR ###################################################################################################################################

# Counter PCMER is the frequency and PCMER is the codify in hashmap
# Count the kmers in the structural encoding          
    def counter_pcmer(self, pp, ak, ws, maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws):
        if tuple(pp) in dicFreq_pp.keys(): dicFreq_pp[tuple(pp)] += 1
        else: dicFreq_pp.update({tuple(pp):1})
        if maxFreq_pp < dicFreq_pp[tuple(pp)]: maxFreq_pp = dicFreq_pp[tuple(pp)]

        if tuple(ak) in dicFreq_ak.keys(): dicFreq_ak[tuple(ak)] += 1
        else: dicFreq_ak.update({tuple(ak):1})
        if maxFreq_ak < dicFreq_ak[tuple(ak)]: maxFreq_ak = dicFreq_ak[tuple(ak)]

        if tuple(ws) in dicFreq_ws.keys(): dicFreq_ws[tuple(ws)] += 1
        else: dicFreq_ws.update({tuple(ws):1})
        if maxFreq_ws < dicFreq_ws[tuple(ws)]: maxFreq_ws = dicFreq_ws[tuple(ws)]

        return maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws
    
## maxFreq is calculate in the maximum frequency existing in the PCMER
    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYPCMER(self, maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws, thresholds):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency
        for threshold in thresholds:
            list_kmer_rgb = list()
            if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
            elif threshold == 0: channel_freq = 256
            else: channel_freq = 0 

            ratio_freq_pp = ratioFreq(maxFreq_pp, channel_freq)
            ratio_freq_ak = ratioFreq(maxFreq_ak, channel_freq)
            ratio_freq_ws = ratioFreq(maxFreq_ws, channel_freq)

            index = 0

            if threshold > 0: # structural encoding
                threshold_k = int(math.log2(256 * threshold))
                items_pp = ['R', 'Y']; items_ak = ['M', 'K']; items_ws = ['W', 'S']
                struct_codeRGB_pp_threshold_k = encodingColour(threshold_k, items_pp, threshold); struct_codeRGB_pp_mod = encodingColour(self.kmer%threshold_k, items_pp, threshold)
                struct_codeRGB_ak_threshold_k = encodingColour(threshold_k, items_ak, threshold); struct_codeRGB_ak_mod = encodingColour(self.kmer%threshold_k, items_ak, threshold)
                struct_codeRGB_ws_threshold_k = encodingColour(threshold_k, items_ws, threshold); struct_codeRGB_ws_mod = encodingColour(self.kmer%threshold_k, items_ws, threshold)

                count_div = math.ceil(self.kmer/threshold_k)

            else: threshold_k = 0 # only frequency

            for entry in self.get_seqKmers_PPAKWS():
        
                if entry == 0 : index += 1; continue
                kmer = entry[0]
                pp = entry[1][0]; freq_pp = int(dicFreq_pp[tuple(pp)]); inc_freq_pp = freq_pp * ratio_freq_pp
                ak = entry[1][1]; freq_ak = int(dicFreq_ak[tuple(ak)]); inc_freq_ak = freq_ak * ratio_freq_ak
                ws = entry[1][2]; freq_ws = int(dicFreq_ws[tuple(ws)]); inc_freq_ws = freq_ws * ratio_freq_ws


                # 0 in frequency is reserved for empy k-mer, 255 is reserved for the two k-mer more frequency
                # in structural encoding are distribuited [0,255] value

                if inc_freq_pp > 255: inc_freq_pp = 255
                if inc_freq_ak > 255: inc_freq_ak = 255
                if inc_freq_ws > 255: inc_freq_ws = 255


   
                k = self.kmer
                pp_rest = pp; ak_rest = ak; ws_rest = ws
                start = 0; rgb_pp = 0; rgb_ak = 0; rgb_ws = 0

                if threshold_k > 0:
                    while k > 0: # until of k > 0 devide kmer in number of pieces
                        if k >= threshold_k: 
                            struct_codeRGB_pp = struct_codeRGB_pp_threshold_k
                            struct_codeRGB_ak = struct_codeRGB_ak_threshold_k
                            struct_codeRGB_ws = struct_codeRGB_ws_threshold_k
                        else:
                            struct_codeRGB_pp = struct_codeRGB_pp_mod
                            struct_codeRGB_ak = struct_codeRGB_ak_mod
                            struct_codeRGB_ws = struct_codeRGB_ws_mod
                        

                        if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                        pp_aux = pp_rest[start:start+threshold_k]; ak_aux = ak_rest[start:start+threshold_k]; ws_aux = ws_rest[start:start+threshold_k]
                        pp_rest = pp_rest[start+threshold_k:]; ak_rest = ak_rest[start+threshold_k:]; ws_rest = ws_rest[start+threshold_k:]
                        k -= threshold_k

                        rgb_pp += struct_codeRGB_pp[tuple(pp_aux)] # code rgb for pp_aux
                        rgb_ak += struct_codeRGB_ak[tuple(ak_aux)]
                        rgb_ws += struct_codeRGB_ws[tuple(ws_aux)]

                    rgb_pp /= count_div # devide rgb in number of pieces
                    rgb_pp += inc_freq_pp # sum between frequency and pp structure
                    rgb_ak /= count_div; rgb_ak += inc_freq_ak
                    rgb_ws /= count_div; rgb_ws += inc_freq_ws
            
                else: ## threshold_k == 0 and theshold == 0 (only frequence)
                    rgb_pp += inc_freq_pp; rgb_ak += inc_freq_ak; rgb_ws += inc_freq_ws

            
                list_kmer_rgb.append([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])
                index += 1

                # reset threshold_k for next k-mer in PPAKWS
                if threshold > 0: threshold_k = int(math.log2(256 * threshold))
                else: threshold_k = 0


            self.__set_kmers_rgb(threshold, list_kmer_rgb)

##################################################################################################################################################################

    def __design_fcgr(self, title, directory, type_encodingColour, dir_classes):
        len_square = int(math.sqrt(2**(2*self.kmer)))

        for threshold, kmer_rgb in self.kmers_rgb.items():
            new_list = []; list_rgb = []; count = 0; freq = self.kmer_freq

            directory_th = directory + '/CCGR (k='+ str(self.kmer) +' T=' + str(threshold) + " " + type_encodingColour + ')/' + dir_classes
            if not os.path.exists(directory_th):
                os.makedirs(directory_th)
            print(directory_th)
            
            pcmer_kmer_rgb = list(kmer_rgb)
            
        
            for i in freq:
                count += 1
                kmer = i[0] ## fcgr
                if kmer == 0: 
                    if (threshold == 0 and "kCCGR" in type_encodingColour):# or  (threshold == 0 and "pcCCGR" in type_encodingColour) : 
                        list_rgb.append(0)
                    else: list_rgb.append([0, 0, 0]) 
                else:
                    for j in pcmer_kmer_rgb: ## rgb
                        if kmer == j[0]: 
                            if threshold == 0 and "kCCGR" in type_encodingColour: list_rgb.append(j[1][0]); break # freq fcgr gs
                            #elif threshold == 0 and "pcCCGR" in type_encodingColour: list_rgb.append(j[1][1]); break # pc-mer gs
                            else: list_rgb.append(j[1]); break # rgb
                    


                if count == len_square: new_list.append(list_rgb); count = 0; list_rgb = []
            grid = np.array(new_list, dtype=np.uint8)

            # save in GS/RGB format
            if (threshold == 0 and "kCCGR" in type_encodingColour):# or (threshold == 0 and "pcCCGR" in type_encodingColour):
                im = Image.fromarray(grid).resize((1024,1024), resample=Image.NEAREST)

                plt.subplots(1, 1, figsize=(256/100, 256/100))
                plt.axis('off')
                plt.imshow(im)
                title = title.replace(' ',"")
                plt.margins(x=0, y=0) # for eliminate the image box (vs 2.0)
                plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0, rect=(0,0,0,0)) # vs 2.0
                plt.savefig(directory_th+ '/' + 'CCGR-GS' + title + "(k = "+ str(self.kmer) + ")", pad_inches = 0)
                plt.clf(); plt.close()

                # save in GS format
                img = Image.open(directory_th + '/' + 'CCGR-GS' + title + "(k = "+ str(self.kmer) + ").png").convert('L')
                img.save(directory_th + '/' + 'CCGR-GS' + title + "(k = "+ str(self.kmer) + ").png")

            else:
                im = Image.fromarray(grid).resize((1024,1024), resample=Image.NEAREST).convert('RGB')

                plt.subplots(1, 1, figsize=(256/100, 256/100))
                plt.axis('off')
                plt.imshow(im)
                title = title.replace(' ',"")
                plt.margins(x=0, y=0) # for eliminate the image box (vs 2.0)
                plt.tight_layout(pad=0.0, h_pad=0.0, w_pad=0.0, rect=(0,0,0,0)) # vs 2.0
                plt.savefig(directory_th + '/CCGR-RGB' + title + "(k = "+ str(self.kmer) + ")", pad_inches = 0)
                plt.clf(); plt.close()

                im = Image.open(directory_th + '/CCGR-RGB' + title + "(k = "+ str(self.kmer) + ")" + ".png")

                background = Image.new("RGB", im.size, (255, 255, 255))
                background.paste(im, mask = im.split()[3])
                background.save(directory_th + '/CCGR-RGB' + title + "(k = "+ str(self.kmer) + ").png", "PNG", quality=100)


        return self
    

    def build_fcgr_pcmerMAXFREQUENCYFCGR(self, seq, k, thresholds, jellyfish, fileFASTA):

        if jellyfish == False: kmers, frequencies = count_kmers(seq, k)
        else: kmers, frequencies = count_kmers_jellyfish(fileFASTA, k)

        maxFreq = 0; array_size = int(math.sqrt(2**(2*k))) # cells for row/colomn present in the fcgr bi-dimensional matrix
        super().init_seq_kmer(array_size)

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
            
            
            j -= 1 # array_size-1..0 
          matrix_fcgr[j][i] = freq # [j][i] coord y, x = riga j, colonna i
          seq_fcgr[j][i] = seq_kmer

          index = j * array_size + i
          # PCMER codify
          super().build_pcmer(seq_kmer, index)
          self.set_fcgr_intofcgr_pcmer_rgb(seq, k, index, seq_kmer, freq)

          if freq > maxFreq: maxFreq = freq

        self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR(maxFreq, thresholds)


        return self
    

    def build_fcgr_pcmerMAXFREQUENCYPCMER(self, seq, k, thresholds, jellyfish, fileFASTA):
 
        maxFreq_pp = 0; maxFreq_ak = 0; maxFreq_ws = 0
        dicFreq_pp = {}; dicFreq_ak = {}; dicFreq_ws = {}

        if jellyfish == False: kmers, frequencies = count_kmers(seq, k)
        else: kmers, frequencies = count_kmers_jellyfish(fileFASTA, k)

        array_size = int(math.sqrt(2**(2*k))) # cells for row/colomn present in the fcgr bi-dimensional matrix
        super().init_seq_kmer(array_size)

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
          # PCMER codify
          super().build_pcmer(seq_kmer, index)
          self.set_fcgr_intofcgr_pcmer_rgb(seq, k, index, seq_kmer, freq)

          pp = self.seqKmers_PPAKWS[index][1][0]
          ak = self.seqKmers_PPAKWS[index][1][1]
          ws = self.seqKmers_PPAKWS[index][1][2]

          maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws = self.counter_pcmer(pp, ak, ws, maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws)
        self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYPCMER(maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws, thresholds)
        return self
    

    def build_fcgr(self, seq, k, type_encodingColour, thresholds, jellyfish, fileFASTA, title = None, directory = None, dir_classes = None, flag_design_gs = False):
        if "kCCGR" in type_encodingColour:
            self.build_fcgr_pcmerMAXFREQUENCYFCGR(seq, k, thresholds, jellyfish, fileFASTA)
        elif "pcCCGR" in type_encodingColour:
            self.build_fcgr_pcmerMAXFREQUENCYPCMER(seq, k, thresholds, jellyfish, fileFASTA)
        self.__design_fcgr(title, directory, type_encodingColour, dir_classes)
        return self  
    
    
