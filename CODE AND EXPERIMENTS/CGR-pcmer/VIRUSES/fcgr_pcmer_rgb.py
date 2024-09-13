import sys
sys.path.insert(1, 'CODE AND EXPERIMENTS/CGR-pcmer/')
import VIRUSES

from VIRUSES.module import *
from VIRUSES.pcmer import PCmer
from VIRUSES.fcgr import FCGR
from VIRUSES.cgr import CGR
from VIRUSES.functions_CGR_PCMER import encodingColour, encodingColourNORIPETITIONS, ratioFreq, findFreqKmer, countPresence, max_freqPCMER, count_kmers, count_kmers_jellyfish

class FCGR_PCMER_RGB(FCGR, PCmer):
    def __init__(self, n, seq = "", kmer = 0, kmer_freq = list(), kmers_rgb = {}):
        FCGR.__init__(self, seq = "", kmer = 0, kmer_freq = list())
        PCmer.__init__(self, n)
        self.kmers_rgb = kmers_rgb

    def __set_kmers_rgb(self, threshold, list_kmer_rgb):
        self.kmers_rgb.update({threshold : list_kmer_rgb})


    def get_kmers_rgb(self):
        return self.kmers_rgb
    
#### Frequency FCGR ##################################################################################################################################
       
## maxFreq is calculate in the maximum frequency existing in the FCGR (1 version)

    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR(self, maxFreq, thresholds):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency


        for threshold in thresholds:
            #print('threshold', threshold)
            if threshold == 0: start0_time = time.time()
            if threshold == 0.5: start05_time = time.time()
            if threshold == 1: start1_time = time.time()
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

            if threshold == 0: end0_time = time.time()
            if threshold == 0.5: end05_time = time.time()
            if threshold == 1: end1_time = time.time()
            if threshold == 0: print('Each kmers rgb 0 codify', end0_time - start0_time)
            if threshold == 0.5: print('Each kmers rgb 0.5 codify', end05_time - start05_time)
            if threshold == 1: print('Each kmers rgb 1 codify', end1_time - start1_time)


            self.__set_kmers_rgb(threshold, list_kmer_rgb)

    ## maxFreq is calculate in the maximum frequency existing of the kmers (2 version)
    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYKMERS(self, threshold_range = []):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency

        maxFreq = len(self.seq) - self.kmer + 1 # n - k + 1
        for threshold in threshold_range:
            list_kmer_rgb = list()
            if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
            elif threshold == 0: channel_freq = 256
            else: channel_freq = 0
            ratio_freq = ratioFreq(maxFreq, channel_freq) 
            for entry in self.get_seqKmers_PPAKWS():
            
                kmer = entry[0]
                pp = entry[1][0]
                ak = entry[1][1]
                ws = entry[1][2]
                freq = findFreqKmer(kmer, self.kmer_freq)

                inc_freq = freq * ratio_freq - 1 # freq * (channel_freq / n - k + 1)

                if threshold > 0: threshold_k = int(math.log2(256 * threshold))
                else: threshold_k = 0

                k = self.kmer
                pp_rest = pp; ak_rest = ak; ws_rest = ws
                start = 0; count_div = 0; rgb_pp = 0; rgb_ak = 0; rgb_ws = 0
                items_pp = ['R', 'Y']; items_ak = ['M', 'K']; items_ws = ['W', 'S']
                if threshold_k > 0:
                    while k > 0: # until of k > 0 devide kmer in number of pieces
                        if k >= threshold_k: 
                            struct_codeRGB_pp = encodingColour(threshold_k, items_pp, threshold)
                            struct_codeRGB_ak = encodingColour(threshold_k, items_ak, threshold)
                            struct_codeRGB_ws = encodingColour(threshold_k, items_ws, threshold)
                        else: 
                            struct_codeRGB_pp = encodingColour(k, items_pp, threshold)
                            struct_codeRGB_ak = encodingColour(k, items_ak, threshold)
                            struct_codeRGB_ws = encodingColour(k, items_ws, threshold)

                        count_div += 1

                        if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                        pp_aux = pp_rest[start:start+threshold_k]; ak_aux = ak_rest[start:start+threshold_k]; ws_aux = ws_rest[start:start+threshold_k]
                        pp_rest = pp_rest[start+threshold_k:]; ak_rest = pp_rest[start+threshold_k:]; pp_rest = ws_rest[start+threshold_k:]
                        k -= threshold_k

                        rgb_pp += struct_codeRGB_pp[tuple(pp_aux)] # code rgb for pp_aux
                        rgb_ak += struct_codeRGB_ak[tuple(pp_aux)]
                        rgb_ws += struct_codeRGB_ws[tuple(pp_aux)]

                    rgb_pp /= count_div # devide rgb in number of pieces
                    rgb_ak /= count_div
                    rgb_ws /= count_div

                    rgb_pp += inc_freq # sum between frequency and pp structure
                    rgb_ak += inc_freq 
                    rgb_ws += inc_freq 

                else: ## threshold_k == 0 and theshold == 0 (only frequence)
                    rgb_pp += inc_freq; rgb_ak += inc_freq; rgb_ws += inc_freq




                list_kmer_rgb.append([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])

            self.__set_kmers_rgb(threshold, list_kmer_rgb)
###
    
    # 1 version
    def __seqKmers_PPAKWS_rgbAllEncodingWeightedAverage_MAXFREQUENCYFCGR(self, threshold):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency

        maxFreq = self.max_freq()
        if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
        elif threshold == 0: channel_freq = 256
        else: channel_freq = 0
        ratio_freq = ratioFreq(maxFreq, channel_freq) 

        for entry in self.get_seqKmers_PPAKWS():
        
            kmer = entry[0]
            pp = entry[1][0]
            ak = entry[1][1]
            ws = entry[1][2]
            freq = findFreqKmer(kmer, self.kmer_freq)

            inc_freq = freq * ratio_freq

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0
   
            ## pp
            k = self.kmer
            pp_rest = pp
            start = 0; n_pieces = 0; rgb_pp = 0; items = ['R', 'Y']
            if threshold_k > 0:
                while k > 0: # until of k > 0 devide kmer in number of pieces
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k
  

                    if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                    pp_aux = pp_rest[start:start+threshold_k]
                    pp_rest = pp_rest[start+threshold_k:]
                    k -= threshold_k

                    rgb_pp += struct_codeRGB[tuple(pp_aux)] * n_pieces  # code rgb for pp_aux

                rgb_pp /= self.kmer # devide rgb in number of pieces (k)
                rgb_pp += inc_freq # sum between frequency and pp structure
            
            else: ## threshold_k == 0 and theshold == 0 (only frequence)
                rgb_pp += inc_freq

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ak
            k = self.kmer
            ak_rest = ak
            start = 0; n_pieces = 0; rgb_ak = 0; items = ['M', 'K']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k

                    if k-threshold_k <= 0: threshold_k = k
                    ak_aux = ak_rest[start:start+threshold_k]
                    ak_rest = ak_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ak += struct_codeRGB[tuple(ak_aux)] * n_pieces

                rgb_ak /= self.kmer
                rgb_ak += inc_freq # sum between frequency and ak structure

            else: ## threshold_k == 0 and theshold == 0
                rgb_ak += inc_freq


            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ws
            k = self.kmer
            ws_rest = ws
            start = 0; n_pieces = 0; rgb_ws = 0; items = ['W', 'S']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k

                    if k-threshold_k <= 0: threshold_k = k
                    ws_aux = ws_rest[start:start+threshold_k]
                    ws_rest = ws_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ws += struct_codeRGB[tuple(ws_aux)] * n_pieces


                rgb_ws /= self.kmer
                rgb_ws += inc_freq # sum between frequency and ws structure (only if threshold != 1)

            else: ## threshold_k == 0 and theshold == 0
                rgb_ws += inc_freq

            self.__set_kmers_rgb([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])


    # 2 version
    def __seqKmers_PPAKWS_rgbAllEncodingWeightedAverage_MAXFREQUENCYKMERS(self, threshold):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency

        maxFreq = len(self.seq) - self.kmer + 1 # n - k + 1
        if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
        elif threshold == 0: channel_freq = 256
        else: channel_freq = 0
        ratio_freq = ratioFreq(maxFreq, channel_freq) 
        for entry in self.get_seqKmers_PPAKWS():
        
            kmer = entry[0]
            pp = entry[1][0]
            ak = entry[1][1]
            ws = entry[1][2]
            freq = findFreqKmer(kmer, self.kmer_freq)

            inc_freq = freq * ratio_freq # freq * (channel_freq / n - k + 1)

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0
   
            ## pp
            k = self.kmer
            pp_rest = pp
            start = 0; n_pieces = 0; rgb_pp = 0; items = ['R', 'Y']
            if threshold_k > 0:
                while k > 0: # until of k > 0 devide kmer in number of pieces
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k
  

                    if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                    pp_aux = pp_rest[start:start+threshold_k]
                    pp_rest = pp_rest[start+threshold_k:]
                    k -= threshold_k

                    rgb_pp += struct_codeRGB[tuple(pp_aux)] * n_pieces  # code rgb for pp_aux

                rgb_pp /= self.kmer # devide rgb in number of pieces (k)
                rgb_pp += inc_freq # sum between frequency and pp structure
            
            else: ## threshold_k == 0 and theshold == 0 (only frequence)
                rgb_pp += inc_freq

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ak
            k = self.kmer
            ak_rest = ak
            start = 0; n_pieces = 0; rgb_ak = 0; items = ['M', 'K']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k

                    if k-threshold_k <= 0: threshold_k = k
                    ak_aux = ak_rest[start:start+threshold_k]
                    ak_rest = ak_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ak += struct_codeRGB[tuple(ak_aux)] * n_pieces

                rgb_ak /= self.kmer
                rgb_ak += inc_freq # sum between frequency and ak structure

            else: ## threshold_k == 0 and theshold == 0
                rgb_ak += inc_freq


            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ws
            k = self.kmer
            ws_rest = ws
            start = 0; n_pieces = 0; rgb_ws = 0; items = ['W', 'S']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold); n_pieces = threshold_k
                    else: struct_codeRGB = encodingColour(k, items, threshold); n_pieces = k

                    if k-threshold_k <= 0: threshold_k = k
                    ws_aux = ws_rest[start:start+threshold_k]
                    ws_rest = ws_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ws += struct_codeRGB[tuple(ws_aux)] * n_pieces


                rgb_ws /= self.kmer
                rgb_ws += inc_freq # sum between frequency and ws structure (only if threshold != 1)

            else: ## threshold_k == 0 and theshold == 0
                rgb_ws += inc_freq

            self.__set_kmers_rgb([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])




###
            
    # 1 version
    def __seqKmers_PPAKWS_rgbOnlyPresences_MAXFREQUENCYFCGR(self, threshold_range = []):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency
        ## if k>50 memory error risk
        for threshold in threshold_range:
            list_kmer_rgb = list()
            max_k = 256 * threshold - 1
            if threshold > 0:
                    presences_k = self.kmer+1 # the presences in k is k+1
            else:
                presences_k = 0


            maxFreq = self.max_freq()
            if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
            elif threshold == 0: channel_freq = 256
            else: channel_freq = 0
            ratio_freq = ratioFreq(maxFreq, channel_freq) 

            for entry in self.get_seqKmers_PPAKWS():
                kmer = entry[0]
                pp = entry[1][0]
                ak = entry[1][1]
                ws = entry[1][2]
                freq = findFreqKmer(kmer, self.kmer_freq)

                inc_freq = freq * ratio_freq
    
                ## pp
                k = self.kmer
                pp_rest = pp, ak_rest = ak; ws_rest = ws
                start = 0; count_div = 0; gb_pp = 0; rgb_ak = 0; rgb_ws = 0
                items_pp = ['R', 'Y']; items_ak = ['M', 'K']; items_ws = ['W', 'S']
                if threshold > 0:
                    while k > 0: # until of k > 0 devide kmer in number of pieces
                        if k < max_k: 
                            struct_codeRGB_pp = encodingColourNORIPETITIONS(k, presences_k, items_pp, threshold) # k < 255 -> presences_k = k+1 <= [0..255]
                            struct_codeRGB_ak = encodingColourNORIPETITIONS(k, presences_k, items_ak, threshold)
                            struct_codeRGB_ws = encodingColourNORIPETITIONS(k, presences_k, items_ws, threshold)
                            max_k = k
                        else: 
                            struct_codeRGB_pp = encodingColourNORIPETITIONS(k, max_k+1, items_pp, threshold)
                            struct_codeRGB_ak = encodingColourNORIPETITIONS(k, max_k+1, items_ak, threshold)
                            struct_codeRGB_ws = encodingColourNORIPETITIONS(k, max_k+1, items_ws, threshold)

                        count_div += 1

                        pp_aux = pp_rest[start:start+max_k]; ak_aux = ak_rest[start:start+max_k]; ws_aux = ws_rest[start:start+max_k]
                        pp_rest = pp_rest[start+max_k:]; ak_rest = ak_rest[start+max_k:]; ws_rest = ws_rest[start+max_k:]
                        k -= max_k

                        count_oneBase_pp, count_twoBase_pp, char_oneBase_pp, char_twoBase_pp = countPresence(pp_aux)
                        count_oneBase_ak, count_twoBase_ak, char_oneBase_ak, char_twoBase_ak = countPresence(ak_aux)
                        count_oneBase_ws, count_twoBase_ws, char_oneBase_ws, char_twoBase_ws = countPresence(ws_aux)

                        pp_idx = tuple([str(count_oneBase_pp) + char_oneBase_pp, str(count_twoBase_pp) + char_twoBase_pp])
                        if pp_idx not in struct_codeRGB_pp.keys(): pp_idx = tuple([str(count_twoBase_pp) + char_twoBase_pp, str(count_oneBase_pp) + char_oneBase_pp])
                        rgb_pp += struct_codeRGB_pp[pp_idx] # code rgb for pp_aux

                        ak_idx = tuple([str(count_oneBase_ak) + char_oneBase_ak, str(count_twoBase_ak) + char_twoBase_ak])
                        if ak_idx not in struct_codeRGB_ak.keys(): ak_idx = tuple([str(count_twoBase_ak) + char_twoBase_ak, str(count_oneBase_ak) + char_oneBase_ak])
                        rgb_ak += struct_codeRGB_ak[ak_idx]

                        ws_idx = tuple([str(count_oneBase_ws) + char_oneBase_ws, str(count_twoBase_ws) + char_twoBase_ws])
                        if ws_idx not in struct_codeRGB_ws.keys(): ws_idx = tuple([str(count_twoBase_ws) + char_twoBase_ws, str(count_oneBase_ws) + char_oneBase_ws])
                        rgb_ws += struct_codeRGB_ws[ws_idx]

                    rgb_pp /= count_div # devide rgb in number of pieces
                    rgb_ak /= count_div
                    rgb_ws /= count_div

                    rgb_pp += inc_freq # sum between frequency and pp structure
                    rgb_ak += inc_freq
                    rgb_ws += inc_freq

                else: ## theshold == 0 (only frequence)
                    rgb_pp += inc_freq; rgb_ak += inc_freq; rgb_ws += inc_freq


                
                list_kmer_rgb.append([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])

            end_time = time.time()
            print('Each kmers rgb codify', end_time - start_time)


            self.__set_kmers_rgb(threshold, list_kmer_rgb)


# 2 version
    def __seqKmers_PPAKWS_rgbOnlyPresences_MAXFREQUENCYKMERS(self, threshold = 0):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency
        ## if k>50 memory error risk
        
        max_k = 256 * threshold - 1
        if threshold > 0:
                presences_k = self.kmer+1 # the presences in k is k+1
        else:
            presences_k = 0


        maxFreq = len(self.seq) - self.kmer + 1 # n - k + 1
        if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
        elif threshold == 0: channel_freq = 256
        else: channel_freq = 0
        ratio_freq = ratioFreq(maxFreq, channel_freq) 
        for entry in self.get_seqKmers_PPAKWS():

            kmer = entry[0]
            pp = entry[1][0]
            ak = entry[1][1]
            ws = entry[1][2]
            freq = findFreqKmer(kmer, self.kmer_freq)

            inc_freq = freq * ratio_freq
   
            ## pp
            k = self.kmer
            pp_rest = pp
            start = 0; count_div = 0; rgb_pp = 0; items = ['R', 'Y']
            if threshold > 0:
                while k > 0: # until of k > 0 devide kmer in number of pieces
                    if k < max_k: struct_codeRGB = encodingColourNORIPETITIONS(k, presences_k, items, threshold); max_k = k # k < 255 -> presences_k = k+1 <= [0..255]
                    else: struct_codeRGB = encodingColourNORIPETITIONS(k, max_k+1, items, threshold)
                    count_div += 1

                    pp_aux = pp_rest[start:start+max_k]
                    pp_rest = pp_rest[start+max_k:]
                    k -= max_k

                    count_oneBase, count_twoBase, char_oneBase, char_twoBase = countPresence(pp_aux)
                    pp_idx = tuple([str(count_oneBase) + char_oneBase, str(count_twoBase) + char_twoBase])
                    if pp_idx not in struct_codeRGB.keys(): pp_idx = tuple([str(count_twoBase) + char_twoBase, str(count_oneBase) + char_oneBase])
                    rgb_pp += struct_codeRGB[pp_idx] # code rgb for pp_aux

                rgb_pp /= count_div # devide rgb in number of pieces
                rgb_pp += inc_freq # sum between frequency and pp structure
            
            else: ## theshold == 0 (only frequence)
                rgb_pp += inc_freq

            ## ak
            k = self.kmer
            ak_rest = ak
            start = 0; count_div = 0; rgb_ak = 0; items = ['M', 'K']
            if threshold > 0:
                while k > 0:
                    if k < max_k: struct_codeRGB = encodingColourNORIPETITIONS(k, presences_k, items, threshold); max_k = k
                    else: struct_codeRGB = encodingColourNORIPETITIONS(k, max_k+1, items, threshold)
                    count_div += 1

                    ak_aux = ak_rest[start:start+max_k]
                    ak_rest = ak_rest[start+max_k:]
                    k -= max_k

                    count_oneBase, count_twoBase, char_oneBase, char_twoBase = countPresence(ak_aux)
                    ak_idx = tuple([str(count_oneBase) + char_oneBase, str(count_twoBase) + char_twoBase])
                    if ak_idx not in struct_codeRGB.keys(): ak_idx = tuple([str(count_twoBase) + char_twoBase, str(count_oneBase) + char_oneBase])
                    rgb_ak += struct_codeRGB[ak_idx]

                rgb_ak /= count_div
                rgb_ak += inc_freq # sum between frequency and ak structure

            else: ## threshold_k == 0 and theshold == 0
                rgb_ak += inc_freq

            ## ws
            k = self.kmer
            ws_rest = ws
            start = 0; count_div = 0; rgb_ws = 0; items = ['W', 'S']
            if threshold > 0:
                while k > 0:
                    if k < max_k: struct_codeRGB = encodingColourNORIPETITIONS(k, presences_k, items, threshold); max_k = k
                    else: struct_codeRGB = encodingColourNORIPETITIONS(k, max_k+1, items, threshold)
                    count_div += 1

                    ws_aux = ws_rest[start:start+max_k]
                    ws_rest = ws_rest[start+max_k:]
                    k -= max_k

                    count_oneBase, count_twoBase, char_oneBase, char_twoBase = countPresence(ws_aux)
                    ws_idx = tuple([str(count_oneBase) + char_oneBase, str(count_twoBase) + char_twoBase])
                    if ws_idx not in struct_codeRGB.keys(): ws_idx = tuple([str(count_twoBase) + char_twoBase, str(count_oneBase) + char_oneBase])
                    rgb_ws += struct_codeRGB[tuple(ws_idx)]

                rgb_ws /= count_div
                rgb_ws += inc_freq # sum between frequency and ws structure

            else: ## threshold_k == 0 and theshold == 0
                rgb_ws += inc_freq

            self.__set_kmers_rgb([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])

#### Frequency PCMER ###################################################################################################################################
            
# Counter PCMER is the frequency and PCMER is the codify in hashmap
    
## maxFreq is calculate in the maximum frequency existing in the PCMER (1 version)
    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYPCMER(self, maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws, thresholds):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency
        for threshold in thresholds:
            if threshold == 0: start0_time = time.time()
            if threshold == 0.5: start05_time = time.time()
            if threshold == 1: start1_time = time.time()
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
                
                
            if threshold == 0: end0_time = time.time()
            if threshold == 0.5: end05_time = time.time()
            if threshold == 1: end1_time = time.time()
            if threshold == 0: print('Each kmers rgb 0 codify', end0_time - start0_time)
            if threshold == 0.5: print('Each kmers rgb 0.5 codify', end05_time - start05_time)
            if threshold == 1: print('Each kmers rgb 1 codify', end1_time - start1_time)


            self.__set_kmers_rgb(threshold, list_kmer_rgb)



     ## maxFreq is calculate in the maximum frequency existing of the pcmer (2 version)
    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYKMERSPCMER(self, threshold = 0):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency
        dict_pp, dict_ak, dict_ws = self.counter_pcmer()

        maxFreq = len(self.seq) - self.kmer + 1 # n - k + 1
        if threshold > 0 and threshold < 1: channel_freq = 256-2**int(math.log2(256 * threshold))
        elif threshold == 0: channel_freq = 256
        else: channel_freq = 0

        ratio_freq = ratioFreq(maxFreq, channel_freq)

        for entry in self.get_seqKmers_PPAKWS():
        
            kmer = entry[0]
            pp = entry[1][0]; freq_pp = int(dict_pp[tuple(pp)])
            ratio_freq = ratioFreq(maxFreq, channel_freq); inc_freq_pp = freq_pp * ratio_freq
            ak = entry[1][1]; freq_ak = int(dict_ak[tuple(ak)])
            ratio_freq = ratioFreq(maxFreq, channel_freq); inc_freq_ak = freq_ak * ratio_freq
            ws = entry[1][2]; freq_ws = int(dict_ws[tuple(ws)])
            ratio_freq = ratioFreq(maxFreq, channel_freq); inc_freq_ws = freq_ws * ratio_freq

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0
   
            ## pp
            k = self.kmer
            pp_rest = pp
            start = 0; count_div = 0; rgb_pp = 0; items = ['R', 'Y']
            if threshold_k > 0:
                while k > 0: # until of k > 0 devide kmer in number of pieces
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold)
                    else: struct_codeRGB = encodingColour(k, items, threshold)
                    count_div += 1

                    if k-threshold_k <= 0: threshold_k = k # index set of pp_aux where k < threshold_k
                    pp_aux = pp_rest[start:start+threshold_k]
                    pp_rest = pp_rest[start+threshold_k:]
                    k -= threshold_k

                    rgb_pp += struct_codeRGB[tuple(pp_aux)] # code rgb for pp_aux

                rgb_pp /= count_div # devide rgb in number of pieces
                rgb_pp += inc_freq_pp # sum between frequency and pp structure
            
            else: ## threshold_k == 0 and theshold == 0 (only frequence)
                rgb_pp += inc_freq_pp

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ak
            k = self.kmer
            ak_rest = ak
            start = 0; count_div = 0; rgb_ak = 0; items = ['M', 'K']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold)
                    else: struct_codeRGB = encodingColour(k, items, threshold)
                    count_div += 1

                    if k-threshold_k <= 0: threshold_k = k
                    ak_aux = ak_rest[start:start+threshold_k]
                    ak_rest = ak_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ak += struct_codeRGB[tuple(ak_aux)]

                rgb_ak /= count_div
                rgb_ak += inc_freq_ak # sum between frequency and ak structure

            else: ## threshold_k == 0 and theshold == 0
                rgb_ak += inc_freq_ak

            if threshold > 0: threshold_k = int(math.log2(256 * threshold))
            else: threshold_k = 0

            ## ws
            k = self.kmer
            ws_rest = ws
            start = 0; count_div = 0; rgb_ws = 0; items = ['W', 'S']
            if threshold_k > 0:
                while k > 0:
                    if k >= threshold_k: struct_codeRGB = encodingColour(threshold_k, items, threshold)
                    else: struct_codeRGB = encodingColour(k, items, threshold)
                    count_div += 1

                    if k-threshold_k <= 0: threshold_k = k
                    ws_aux = ws_rest[start:start+threshold_k]
                    ws_rest = ws_rest[start+threshold_k:]
                    k -= threshold_k
                    rgb_ws += struct_codeRGB[tuple(ws_aux)]


                rgb_ws /= count_div
                rgb_ws += inc_freq_ws # sum between frequency and ws structure (only if threshold != 1)

            else: ## threshold_k == 0 and theshold == 0
                rgb_ws += inc_freq_ws

            self.__set_kmers_rgb([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])


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
    
##Frequency FCGR;  Codify PCMER Frequency FCGR (with max frequency FCGR - 1 version) #################################################################################################################################################
    

    ## maxFreq is calculate in the maximum frequency existing in the FCGR (1 version)
    def __seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR_MAXFREQUENCYPCMER(self, threshold = 0):
        # threshold = 0, no pcmer
        # threshold = 1, no frequency

        maxFreq = self.max_freq()
        if threshold > 0 and threshold < 1: channel_freq_FCGR = 256-2**int(math.log2(256 * threshold)) 
        elif threshold == 0: channel_freq_FCGR = 256
        else: channel_freq_FCGR = 0 
        channel_freq_PCMER = 256 - channel_freq_FCGR

        ratio_freq_FCGR = ratioFreq(maxFreq, channel_freq_FCGR)
        dict_pp, dict_ak, dict_ws = self.counter_pcmer() # counts the kmers in the chemical-physical structure for pcmer's frequency

        for entry in self.get_seqKmers_PPAKWS():
        
            if entry == 0 : continue
            kmer = entry[0]
            pp = entry[1][0]
            ak = entry[1][1]
            ws = entry[1][2]

            freq = findFreqKmer(kmer, self.kmer_freq)
            inc_freq_FCGR = freq * ratio_freq_FCGR

            maxFreqPCMER_pp = max_freqPCMER(dict_pp); freqPCMER_pp = int(dict_pp[tuple(pp)])
            ratio_freqPCMER_pp = ratioFreq(maxFreqPCMER_pp, channel_freq_PCMER); inc_freqPCMER_pp = freqPCMER_pp * ratio_freqPCMER_pp
            maxFreqPCMER_ak = max_freqPCMER(dict_ak); freqPCMER_ak = int(dict_ak[tuple(ak)])
            ratio_freqPCMER_ak = ratioFreq(maxFreqPCMER_ak, channel_freq_PCMER); inc_freqPCMER_ak = freqPCMER_ak * ratio_freqPCMER_ak
            maxFreqPCMER_ws = max_freqPCMER(dict_ws); freqPCMER_ws = int(dict_ws[tuple(ws)])
            ratio_freqPCMER_ws = ratioFreq(maxFreqPCMER_ws, channel_freq_PCMER); inc_freqPCMER_ws = freqPCMER_ws * ratio_freqPCMER_ws

            
            ## pp
            rgb_pp = 0
            rgb_pp += inc_freqPCMER_pp # pcmer frequency
            rgb_pp += inc_freq_FCGR # fcgr frequency

            ## ak
            rgb_ak = 0
            rgb_ak += inc_freqPCMER_ak # pcmer frequency
            rgb_ak += inc_freq_FCGR # fcgr frequency
            
            ## ws
            rgb_ws = 0
            rgb_ws += inc_freqPCMER_ws # pcmer frequency
            rgb_ws += inc_freq_FCGR # fcgr frequency

            self.__set_kmers_rgb([kmer, [round(rgb_pp), round(rgb_ak), round(rgb_ws)]])

 

    def seqKmers_PPAKWS_rgb(self, maxFreq, type_encodingColour, thresholds): 
        # fcgr frequency (structural)
        if type_encodingColour == "all encoding MAXFREQUENCYFCGR": self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR(maxFreq, thresholds)
        elif type_encodingColour == "all encoding MAXFREQUENCYKMERS": self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYKMERS(thresholds)
        elif type_encodingColour == "only presences MAXFREQUENCYFCGR": self.__seqKmers_PPAKWS_rgbOnlyPresences_MAXFREQUENCYFCGR(thresholds)
        elif type_encodingColour == "only presences MAXFREQUENCYKMERS": self.__seqKmers_PPAKWS_rgbOnlyPresences_MAXFREQUENCYKMERS(thresholds)
        elif type_encodingColour == "all encoding WA MAXFREQUENCYFCGR": self.__seqKmers_PPAKWS_rgbAllEncodingWeightedAverage_MAXFREQUENCYFCGR(thresholds)
        elif type_encodingColour == "all encoding WA MAXFREQUENCYKMERS": self.__seqKmers_PPAKWS_rgbAllEncodingWeightedAverage_MAXFREQUENCYKMERS(thresholds)
        # pcmer frequency (structural)
        elif type_encodingColour == "all encoding MAXFREQUENCYPCMER": self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYPCMER(thresholds)
        elif type_encodingColour == "all encoding MAXFREQUENCYKMERSPCMER": self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYKMERSPCMER(thresholds)
        ## fcgr frequency and pcmer frequency
        elif type_encodingColour == "all encoding MAXFREQUENCYFCGR_MAXFREQUENCYPCMER": self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR_MAXFREQUENCYPCMER(thresholds)
        

    def __design_fcgr(self, title, directory, type_encodingColour, dir_classes):
        start = time.time()
        #print(self.kmers_rgb)
        len_square = int(math.sqrt(2**(2*self.kmer)))

        for threshold, kmer_rgb in self.kmers_rgb.items():
            new_list = []; list_rgb = []; count = 0; freq = self.kmer_freq

            directory_th = directory + '/FCGR RGB-PCMER-FREQ OPTIMIZATION (k='+ str(self.kmer) +' threshold=' + str(threshold) + " " + type_encodingColour + ')/' + dir_classes
            if not os.path.exists(directory_th):
                os.makedirs(directory_th)
            print(directory_th)
            
            pcmer_kmer_rgb = list(kmer_rgb)
            
        
            for i in freq:
                count += 1
                kmer = i[0] ## fcgr
                if kmer == 0: 
                    if (threshold == 0 and "MAXFREQUENCYFCGR" in type_encodingColour):# or  (threshold == 0 and "MAXFREQUENCYPCMER" in type_encodingColour) : 
                        list_rgb.append(0)
                    else: list_rgb.append([0, 0, 0]) 
                else:
                    for j in pcmer_kmer_rgb: ## rgb
                        if kmer == j[0]: 
                            if threshold == 0 and "MAXFREQUENCYFCGR" in type_encodingColour: list_rgb.append(j[1][0]); break # freq fcgr gs
                            #elif threshold == 0 and "MAXFREQUENCYPCMER" in type_encodingColour: list_rgb.append(j[1][1]); break # pc-mer gs
                            else: list_rgb.append(j[1]); break # rgb
                    


                if count == len_square: new_list.append(list_rgb); count = 0; list_rgb = []
            grid = np.array(new_list, dtype=np.uint8)
            #print('grid', grid)

            # save in GS/RGB format
            if (threshold == 0 and "MAXFREQUENCYFCGR" in type_encodingColour):# or (threshold == 0 and "MAXFREQUENCYPCMER" in type_encodingColour):
                im = Image.fromarray(grid).resize((1024,1024), resample=Image.NEAREST)
                #print('grid GS', grid)
                plt.axis('off')
                plt.imshow(im)
                title = title.replace(' ',"")
                save_fig_start = time.time()
                plt.savefig(directory_th+ '/' + 'FCGR-GS-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ")")
                save_fig_end = time.time()
                plt.clf(); plt.close()
                print('save fig', save_fig_end - save_fig_start)

                # save in GS format
                start_convert = time.time()
                img = Image.open(directory_th + '/' + 'FCGR-GS-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ").png").convert('L')
                img.save(directory_th + '/' + 'FCGR-GS-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ").png")
                end_convert = end.time()
                print('Convert image in GS', end_convert - start_convert)

            else:
                im = Image.fromarray(grid).resize((1024,1024), resample=Image.NEAREST).convert('RGB')
                #print('grid RGB', grid)

                # displaying the title
                #plt.title(title + ', k = ' + str(self.kmer), fontsize=10)
                plt.axis('off')
                plt.imshow(im)
                title = title.replace(' ',"")
                save_fig_start = time.time()
                plt.savefig(directory_th + '/FCGR-RGB-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ")")
                save_fig_end = time.time()
                print('save fig RGB', save_fig_end - save_fig_start)
                plt.clf(); plt.close()

                start_convert = time.time()
                im = Image.open(directory_th + '/FCGR-RGB-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ")" + ".png")

                background = Image.new("RGB", im.size, (255, 255, 255))
                background.paste(im, mask = im.split()[3])
                background.save(directory_th + '/FCGR-RGB-PCMERFREQ' + title + "(k = "+ str(self.kmer) + ").png", "PNG", quality=100)
                end_convert = time.time()
                print('Convert image in RGB', end_convert - start_convert)


        end = time.time()
        print('Time design FCGR ', end - start)
        return self
    

    def build_fcgr_pcmerMAXFREQUENCYFCGR(self, seq, k, thresholds, fileFASTA):
        start0_time = time.time() 
        maxFreq = 0
        kmers, frequencies = count_kmers(seq, k)
        #kmers, frequencies = count_kmers_jellyfish(fileFASTA, k)
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

          if freq > maxFreq: maxFreq = freq

        start_time = time.time()
        self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYFCGR(maxFreq, thresholds)
        end_time = time.time()
        end0_time = time.time()
        print('Time applied colour', end_time - start_time)
        print('Time FCGR PCMER RGB', end0_time - start0_time)

        return self
    

    def build_fcgr_pcmerMAXFREQUENCYPCMER(self, seq, k, thresholds, fileFASTA):
        start0_time = time.time()    
        maxFreq_pp = 0; maxFreq_ak = 0; maxFreq_ws = 0
        dicFreq_pp = {}; dicFreq_ak = {}; dicFreq_ws = {}
        kmers, frequencies = count_kmers(seq, k)
        #kmers, frequencies = count_kmers_jellyfish(fileFASTA, k)
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
        start_time = time.time()
        self.__seqKmers_PPAKWS_rgbAllEncoding_MAXFREQUENCYPCMER(maxFreq_pp, maxFreq_ak, maxFreq_ws, dicFreq_pp, dicFreq_ak, dicFreq_ws, thresholds)
        end_time = time.time()
        end0_time = time.time()
        print('Time applied colour', end_time - start_time)
        print('Time FCGR PCMER RGB', end0_time - start0_time)
        return self
    

    def build_fcgr(self, seq, k, type_encodingColour, thresholds, fileFASTA, title = None, directory = None, dir_classes = None, flag_design_gs = False):
        start_time = time.time()
        if "MAXFREQUENCYFCGR" in type_encodingColour:
            self.build_fcgr_pcmerMAXFREQUENCYFCGR(seq, k, thresholds, fileFASTA)
            end_time = time.time()
        elif "MAXFREQUENCYPCMER" in type_encodingColour:
            self.build_fcgr_pcmerMAXFREQUENCYPCMER(seq, k, thresholds, fileFASTA)
            end_time = time.time()
        #print('build fcgr pcmer rgb ', end_time - start_time)
        self.__design_fcgr(title, directory, type_encodingColour, dir_classes)
        end_time = time.time()
        #print('design fcgr', end_time - start_time)
        return self  
    
    
