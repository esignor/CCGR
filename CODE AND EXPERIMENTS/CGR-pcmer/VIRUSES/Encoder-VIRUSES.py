
## MODULE
import sys
sys.path.insert(1, 'CODE AND EXPERIMENTS/CGR-pcmer/') # my kmer-cout
import VIRUSES

from VIRUSES.module import * 
from VIRUSES.cgr import CGR
from VIRUSES.fcgr import FCGR, FCGR_RGB
from VIRUSES.fcgr_pcmer_rgb import FCGR_PCMER_RGB
from VIRUSES.functions_CGR_PCMER import parse_sequence

from os import listdir
from os.path import isfile, join
import os

if __name__ == '__main__':
    start_time = time.time()
# DATASET HIV-1 (1) n_classes = 12, HIV-1 (2) n_classes = 37, DENGUE n_classes=4, CORONAVIRUSES n_classes=7, HEPATITISC n_classes=6, HEPATITIS B (1) n_classes=8, HEPATITIS B (2) n_classes=13, INFLUENZA A n_classes = 56
    
    jellyfish = False # for used jellyfish set the flag to True
    if jellyfish == True: print('jellyfish multi-threader k-mers counter active')
    else: print('jellyfish multi-threader k-mers counter inactive')

    MAIN_FOLDER = 'CODE AND EXPERIMENTS'
    out_directory = MAIN_FOLDER + '/CGR-pcmer/VIRUSES/ColoredChaosGameRepresentation_ENCODER/Viruses'
    print(out_directory)
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    #folder_Dataset = MAIN_FOLDER + '/DATASET/12classes_hiv dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/37classes_hiv dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/4classes_dengue dataset'
    folder_Dataset = MAIN_FOLDER + '/DATASET/7classes_coronaviruses dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/6classes_hepatitisC dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/8classes_hepatitisB dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/13classes_hepatitisB dataset'
    #folder_Dataset = MAIN_FOLDER + '/DATASET/56classes_influenzaA_reduced dataset'
    

    list_folder_VIRUSES = [folder_Dataset]

    for folder_viruses in list_folder_VIRUSES:

        mypath = folder_viruses
        print(mypath)

        if 'hiv' and '12' in mypath:
            viruses = 'HIV1'
        elif 'hiv' and '37' in mypath:
            viruses = 'HIV2'
        elif 'dengue' in mypath:
            viruses = 'Dengue'
        elif 'coronaviruses' in mypath:
            viruses = 'Coronaviruses'
        elif 'hepatitisC' in mypath:
            viruses = 'HepatitisC'
        elif 'hepatitisB' and '8' in mypath:
            viruses = 'HepatitisB1'
        elif 'hepatitisB' and '13' in mypath:
            viruses = 'HepatitisB2'
        elif 'influenzaA' in mypath:
            viruses = 'InfluenzaA'


        w, h = 1024, 1024
        my_list = os.listdir(mypath)

        for i in range(len(my_list)):
            Subpath=mypath+'/'+my_list[i]
            imgpath=my_list[i]
            onlyfiles = [f for f in listdir(Subpath) if isfile(join(Subpath, f))]
            for s in range(0,len(onlyfiles)):

                data = np.zeros((h, w), dtype=np.uint8)

                fileFASTA = Subpath + '/' + onlyfiles[s]
                genome_viruses = parse_sequence(fileFASTA)
                title = fileFASTA.split('/')[4].split('.fasta')[0].replace('.', '-')
                kmer = 10 # k-mers size: 4, 6, 8, and 10
                type_encodingColour = "kCCGR" # Colour Chaos Game Rapresentation (CCGR), kCCGR e pcCCGR
                threshold = [0, 0.5, 1]
                
                directory_png = out_directory + '/' +  viruses 
                
                if not os.path.exists(directory_png):
                    os.makedirs(directory_png)
                print(directory_png)
                
                FCGR_VIRUSES = FCGR_PCMER_RGB(2**(2*kmer), "", 0, list(), {}).build_fcgr(genome_viruses, kmer, type_encodingColour, threshold, jellyfish, fileFASTA, title, directory_png, my_list[i]) ## --RGB pcmer
                end_time = time.time()
                print('TIME:', end_time - start_time)

                    


