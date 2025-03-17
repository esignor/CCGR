import os
from collections import Counter
import operator

def classes_HepatitisC(y_class):
    if '_' in y_class: return True, "_"
    elif '1' in y_class: return True, "1"
    elif '2' in y_class: return True, "2"
    elif '3' in y_class: return True, "3"
    elif '4' in y_class: return True, "4"
    elif '5' in y_class: return True, "5"
    elif '6' in y_class: return True, "6"
    elif '7' in y_class: return True, "7"
    elif '8' in y_class: return True, "8"

    return False, 0


def prerocessingDataset_InfluenzaA(avg, masterpath, directory, ClassesN_dict, N):
    lines_influenzaA = open(masterpath, "r")
    content = lines_influenzaA.readlines()
    lines_influenzaA.close()
    count_dict = {}

    header = -1; y_class = ""
    for row in content:

        if ">" in row:
            header = row
            header = header.replace("\n", "").replace(">", "_").replace('?', "")
            y_class = ((row.split(' '))[0]).replace('>', "").replace('?',"")


            if y_class not in ClassesN_dict.keys(): y_class = ""; continue
            directory= 'CODE AND EXPERIMENTS/DATASET/' + str(N) + 'classes_influenzaA_reduced dataset/' + y_class
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Create ", y_class,  " directory")
       
            limit = ClassesN_dict[y_class]

            if ClassesN_dict[y_class] > avg:
                limit = avg + ClassesN_dict[y_class]*(10/100)
                
            if y_class in count_dict.keys(): 
                if count_dict[y_class] < limit: count_dict[y_class] += 1
                else: y_class = ""; continue
            else: count_dict.update({y_class : 1})    

            new_fasta_class = open(directory + "/" + header + ".fasta", "w")
            new_fasta_class.write(row)
            new_fasta_class.close()
            
            

        else:
            if y_class not in ClassesN_dict.keys(): continue
            new_fasta_class = open(directory + "/" + header + ".fasta", "a")
            new_fasta_class.write(row)
            new_fasta_class.close()


def prerocessingDataset_HIV(masterpath, directory, ClassesN_dict, N):

    lines_hiv = open(masterpath, "r")
    content = lines_hiv.readlines()
    lines_hiv.close()

    header = -1; y_class = ""
    for row in content:

        if ">" in row:
            header = row
            header = header.replace("\n", "").replace(">", "_")
            y_class = ((row.split('.'))[0]).split('>')[1]

            if y_class not in ClassesN_dict.keys(): y_class = ""; continue
            
            directory= 'CODE AND EXPERIMENTS/DATASET/' + str(N) + 'classes_hiv dataset/' + y_class
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Create ", y_class,  " directory")

            new_fasta_class = open(directory + "/" + header + ".fasta", "w")
            new_fasta_class.write(row)
            new_fasta_class.close()

        else:
            if y_class not in ClassesN_dict.keys(): continue
            new_fasta_class = open(directory + "/" + header + ".fasta", "a")
            new_fasta_class.write(row)
            new_fasta_class.close()


def prerocessingDataset_Dengue(masterpath, directory, Classes4_dict):
    lines_dengue = open(masterpath, "r")
    content = lines_dengue.readlines()
    lines_dengue.close()

    header = -1; y_class = ""
    for row in content:
        if "|" in row:
            header = row
            header = header.replace("\n", "").split(">")[1].replace("|", ".").replace("/", "-")
            y_class = row.split('|')[3].split(" ")[1].replace("\n", "")

            if y_class not in Classes4_dict.keys(): y_class = ""; continue
            
            directory = 'CODE AND EXPERIMENTS/DATASET/4classes_dengue dataset/' + y_class
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Create ", y_class,  " directory")

            new_fasta_class = open(directory + "/" + header + ".fasta", "w")
            new_fasta_class.write(row)
            new_fasta_class.close()

        else:
            if y_class not in Classes4_dict.keys(): continue
            new_fasta_class = open(directory + "/" + header + ".fasta", "a")
            new_fasta_class.write(row)
            new_fasta_class.close()

def prerocessingDataset_HepatitisC(masterpath, directory, ClassesN_dict): 

    lines_hcv = open(masterpath, "r")
    content = lines_hcv.readlines()
    lines_hcv.close()
    header = -1; y_class = ""

    for row in content:
        

        if ">" in row:
            header = row
            header = header.replace("\n", "").split(">")[1].replace("|", ".").replace("/", "-")
            y_class = ((row.split('.'))[0]).split('>')[1]
            bool, y_class = classes_HepatitisC(y_class)
        
            if y_class not in ClassesN_dict.keys(): y_class = ""; continue
            
            directory= 'CODE AND EXPERIMENTS/DATASET/' +'6classes_hepatitisC dataset/' + y_class
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Create ", y_class,  " directory")

            new_fasta_class = open(directory + "/" + header + ".fasta", "w")
            new_fasta_class.write(row)
            new_fasta_class.close()

        else:
            if y_class not in ClassesN_dict.keys(): continue
            new_fasta_class = open(directory + "/" + header + ".fasta", "a")
            new_fasta_class.write(row)
            new_fasta_class.close()
    
def prerocessingDataset_HepatitisB(masterpath, directory, ClassesN_dict, N):

    lines_hepatitisB = open(masterpath, "r")
    content = lines_hepatitisB.readlines()
    lines_hepatitisB.close()

    header = -1; y_class = ""
    for row in content:

        if ">" in row:
            header = row
            header = header.replace(":","").replace(".", "").replace("\n", "").split(">")[1].replace("|", ".").replace("/", "-")
            if N == 13: # HepatitisB_2
                if 'recombinant' in row: y_class = (((row.split('|'))[2]).split('Hepatitis B Virus ')[1].split('.')[0]).split(' (')[0]
                else: y_class = (((row.split('|'))[2]).split('Hepatitis B Virus ')[1].split('.')[0]).split('genotype')[1].split(" ")[1]
            elif N == 8: # HepatitisB_1
                if 'recombinant' not in row: y_class = (((row.split('|'))[2]).split('Hepatitis B Virus ')[1].split('.')[0]).split('genotype')[1].split(" ")[1]
                else: y_class = ""
                
            if y_class not in ClassesN_dict.keys(): y_class = ""; continue
            
            directory= 'CODE AND EXPERIMENTS/DATASET/' + str(N) + 'classes_hepatitisB dataset/' + y_class
            if not os.path.exists(directory):
                os.makedirs(directory)
                print("Create ", y_class,  " directory")

            new_fasta_class = open(directory + "/" + header + ".fasta", "w")
            new_fasta_class.write(row)
            new_fasta_class.close()

        else:
            if y_class not in ClassesN_dict.keys(): continue
            new_fasta_class = open(directory + "/" + header + ".fasta", "a")
            new_fasta_class.write(row)
            new_fasta_class.close()



if __name__ == '__main__':
        # N classes that consider at least of 10% of samples in the subtype class
        flag_virus = "HIV_1"
        #flag_virus = "HIV_2"
        #flag_virus = "Dengue"
        #flag_virus = "HepatitisB_1"
        #flag_virus = "HepatitisB_2"
        #flag_virus = "HepatitisC"
        #flag_virus = "InfluenzaA"
        

        if flag_virus == 'HIV_1':
            directory='CODE AND EXPERIMENTS/DATASET/12classes_hiv dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/hiv-db.fasta'
            N = 12 # 6540
        elif flag_virus == 'HIV_2':
            directory='CODE AND EXPERIMENTS/DATASET/37classes_hiv dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/hiv-db.fasta'
            N = 37 # 7194
        elif flag_virus == 'Dengue':
            directory='CODE AND EXPERIMENTS/DATASET/4classes_dengue dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/dengue.fasta'
            N = 4 # 5079 samples
        elif flag_virus == 'HepatitisC':
            directory='CODE AND EXPERIMENTS/DATASET/9classes_hepatitisC dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/HCV.fasta'
            N = 6 # 2066 samples
        elif flag_virus == 'HepatitisB_1':
            directory='CODE AND EXPERIMENTS/DATASET/8classes_hepatitisB dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/hepatitisB_1.fasta'
            N = 8 # 6138 samples
        elif flag_virus == 'HepatitisB_2':
            directory='CODE AND EXPERIMENTS/DATASET/13classes_hepatitisB dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/hepatitisB_2.fasta'
            N = 13 # 6824 samples
        elif flag_virus == 'InfluenzaA':
            directory='CODE AND EXPERIMENTS/DATASET/56classes_influenzaA_reduced dataset'
            masterpath='CODE AND EXPERIMENTS/DATASET/influenzaA.fasta'
            N = 56 # 36549 samples
            
        if not os.path.exists(directory):
            os.makedirs(directory)

        lines_dengue = open(masterpath, "r")
        content = lines_dengue.readlines()
        lines_dengue.close()

        y_data = []
        for row in content:

            if ">" in row:
                if flag_virus == 'HIV_1' or flag_virus == 'HIV_2':
                    y_class = ((row.split('.'))[0]).split('>')[1]
                    y_data.append(y_class)
                elif flag_virus == 'Dengue':
                    y_class = row.split('|')[3].split(" ")[1].replace('\n', "")
                    y_data.append(y_class)
                elif flag_virus == 'HepatitisC':
                    y_class = ((row.split('.'))[0]).split('>')[1]
                    bool, y_class = classes_HepatitisC(y_class)
                    if y_class == '_': continue
                    if bool == True: y_data.append(y_class)
                elif flag_virus == 'HepatitisB_1':
                    if 'recombinant' in row: continue
                    else:
                        y_class = ((((row.split('|'))[2]).split('Hepatitis B Virus ')[1].split('.')[0]).split('genotype')[1].split(" ")[1])
                        y_data.append(y_class)
                elif flag_virus == 'HepatitisB_2':
                    if 'recombinant' in row: y_class = (((row.split('|'))[2]).split('Hepatitis B Virus ')[1].split('.')[0]).split(' (')[0]
                    else:
                        y_class = (((row.split('|'))[2]).split('Hepatitis B Virus')[1].split('.')[0]).split('genotype')[1].split(" ")[1]
                    y_data.append(y_class)
                elif flag_virus == 'InfluenzaA':
                    y_class = ((row.split(' '))[0]).replace('>', "").replace('?',"")
                    no_classes = ['Mixed', 'mixed']; mixed = False
                    for mix in no_classes:
                        if mix in y_class: mixed = True

                    if mixed == False: y_data.append(y_class)
                    


        cnt = Counter(y_data)
        print(cnt) # count each genomic sequence that belongs of the classes

        # sorted the name of classes for greater numbers of the genomic sequences
        Classes_dict = dict(sorted(cnt.items(), key=operator.itemgetter(1), reverse = True))
        print(Classes_dict)

        ClassesN_dict = dict(); n_classes = 0
        for cl, occ in Classes_dict.items():
            if n_classes >= N: break
            ClassesN_dict.update({cl : occ}); n_classes += 1

        print(ClassesN_dict, len(ClassesN_dict))
        
        n_samples = 0
        for num in ClassesN_dict.values():
            n_samples += num
        
        avg = round(n_samples/N)
        print('n_samples', n_samples)
        print('avg', avg)

  
        if flag_virus == 'HIV_1' or flag_virus == 'HIV_2': prerocessingDataset_HIV(masterpath, directory, ClassesN_dict, N)
        elif flag_virus == 'Dengue': prerocessingDataset_Dengue(masterpath, directory, ClassesN_dict)
        elif flag_virus == 'HepatitisC': prerocessingDataset_HepatitisC(masterpath, directory, ClassesN_dict)
        elif flag_virus == 'HepatitisB_1' or flag_virus == 'HepatitisB_2': prerocessingDataset_HepatitisB(masterpath, directory, ClassesN_dict, N)
        elif flag_virus == 'InfluenzaA': prerocessingDataset_InfluenzaA(avg, masterpath, directory, ClassesN_dict, N)
    
