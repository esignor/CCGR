import os
from collections import Counter
import operator
import argparse

# The FASTA files downloaded from the specific databases, as described in the README, 
# must already be placed in the dev/DATASET path: dengue, HCV, hepatitisB_1, hepatitisB_2, hiv-db, influenzaA

# Terminal command to activate preorocessing:
# in CCGR: python dev/VirusPreprocessingDatasets.py --virus InfluenzaA 
#                                                  (available dataset keyword options: 'HIV1', 'HIV2', 'Dengue', 'HepatitisB1', 'HepatitisB2', 'HepatitisC', 'InfluenzaA')

# VirtualPreprocessingDatasets is responsible for preprocessing the dataset's FASTA file 
# into the appropriate format required by the CCGR software


BASE_DATASET_DIR = 'dev/DATASET/'

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

def prerocessingDataset_InfluenzaA(avg, masterpath, directory, ClassesN_dict, N, formatted_flag):
    with open(masterpath, "r") as f:
        content = f.readlines()
    
    count_dict = {}
    header = -1
    y_class = ""

    for row in content:
        if ">" in row:
            header = row.strip().replace(">", "_").replace('?', "")
            y_class = row.split(' ')[0].replace('>', "").replace('?',"")

            if y_class not in ClassesN_dict.keys():
                y_class = ""
                continue

            directory = os.path.join(BASE_DATASET_DIR, f'{N}classes_{formatted_flag}_dataset', y_class)
            os.makedirs(directory, exist_ok=True)

            limit = min(ClassesN_dict[y_class], avg + round(ClassesN_dict[y_class] * 0.1))

            if y_class in count_dict:
                if count_dict[y_class] < limit:
                    count_dict[y_class] += 1
                else:
                    y_class = ""
                    continue
            else:
                count_dict[y_class] = 1

            with open(os.path.join(directory, header + ".fasta"), "w") as f:
                f.write(row)
        else:
            if y_class not in ClassesN_dict.keys():
                continue
            with open(os.path.join(directory, header + ".fasta"), "a") as f:
                f.write(row)

def prerocessingDataset_HIV(masterpath, directory, ClassesN_dict, N, formatted_flag):
    with open(masterpath, "r") as f:
        content = f.readlines()

    header = -1
    y_class = ""
    for row in content:
        if ">" in row:
            header = row.strip().replace(">", "_")
            y_class = row.split('.')[0].split('>')[1]

            if y_class not in ClassesN_dict:
                y_class = ""
                continue

            directory = os.path.join(BASE_DATASET_DIR, f'{N}classes_{formatted_flag}_dataset', y_class)
            os.makedirs(directory, exist_ok=True)

            with open(os.path.join(directory, header + ".fasta"), "w") as f:
                f.write(row)
        else:
            if y_class not in ClassesN_dict:
                continue
            with open(os.path.join(directory, header + ".fasta"), "a") as f:
                f.write(row)

def prerocessingDataset_Dengue(masterpath, directory, Classes4_dict, N, formatted_flag):
    with open(masterpath, "r") as f:
        content = f.readlines()

    header = -1
    y_class = ""
    for row in content:
        if "|" in row:
            header = row.strip().split(">")[1].replace("|", ".").replace("/", "-")
            y_class = row.split('|')[3].split(" ")[1].strip()

            if y_class not in Classes4_dict:
                y_class = ""
                continue

            directory = os.path.join(BASE_DATASET_DIR, f'{N}classes_{formatted_flag}_dataset', y_class)
            os.makedirs(directory, exist_ok=True)

            with open(os.path.join(directory, header + ".fasta"), "w") as f:
                f.write(row)
        else:
            if y_class not in Classes4_dict:
                continue
            with open(os.path.join(directory, header + ".fasta"), "a") as f:
                f.write(row)

def prerocessingDataset_HepatitisC(masterpath, directory, ClassesN_dict, N, formatted_flag): 
    with open(masterpath, "r") as f:
        content = f.readlines()

    header = -1
    y_class = ""
    for row in content:
        if ">" in row:
            header = row.strip().split(">")[1].replace("|", ".").replace("/", "-")
            y_class = row.split('.')[0].split('>')[1]
            _, y_class = classes_HepatitisC(y_class)

            if y_class not in ClassesN_dict:
                y_class = ""
                continue

            directory = os.path.join(BASE_DATASET_DIR, f'{N}classes_{formatted_flag}_dataset', y_class)
            os.makedirs(directory, exist_ok=True)

            with open(os.path.join(directory, header + ".fasta"), "w") as f:
                f.write(row)
        else:
            if y_class not in ClassesN_dict:
                continue
            with open(os.path.join(directory, header + ".fasta"), "a") as f:
                f.write(row)

def prerocessingDataset_HepatitisB(masterpath, directory, ClassesN_dict, N, formatted_flag):
    with open(masterpath, "r") as f:
        content = f.readlines()

    header = -1
    y_class = ""
    for row in content:
        if ">" in row:
            header = row.replace(":", "").replace(".", "").strip().split(">")[1].replace("|", ".").replace("/", "-")
            if N == 13:
                if 'recombinant' in row:
                    y_class = row.split('|')[2].split('Hepatitis B Virus ')[1].split('.')[0].split(' (')[0]
                else:
                    y_class = row.split('|')[2].split('Hepatitis B Virus ')[1].split('.')[0].split('genotype')[1].split(" ")[1]
            elif N == 8:
                if 'recombinant' not in row:
                    y_class = row.split('|')[2].split('Hepatitis B Virus ')[1].split('.')[0].split('genotype')[1].split(" ")[1]
                else:
                    y_class = ""

            if y_class not in ClassesN_dict:
                y_class = ""
                continue

            directory = os.path.join(BASE_DATASET_DIR, f'{N}classes_{formatted_flag}_dataset', y_class)
            os.makedirs(directory, exist_ok=True)

            with open(os.path.join(directory, header + ".fasta"), "w") as f:
                f.write(row)
        else:
            if y_class not in ClassesN_dict:
                continue
            with open(os.path.join(directory, header + ".fasta"), "a") as f:
                f.write(row)

def normalize_flag_virus(s):
    result = []
    for i, c in enumerate(s):
        if c.isupper():
            result.append(c.lower())
        else:
            result.extend(s[i:]) 
            break
    return ''.join(result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocessing virus dataset')
    parser.add_argument('--virus', type=str, required=True,
                        choices=[
                            'HIV1', 'HIV2', 'Dengue', 'HepatitisB1',
                            'HepatitisB2', 'HepatitisC', 'InfluenzaA'
                        ],
                        help='Virus dataset flag to preprocess')
    args = parser.parse_args()

    flag_virus = args.virus

    master_files = {
        'HIV1': ('hiv-db.fasta', 12),
        'HIV2': ('hiv-db.fasta', 37),
        'Dengue': ('dengue.fasta', 4),
        'HepatitisC': ('HCV.fasta', 6),
        'HepatitisB1': ('hepatitisB.fasta', 8),
        'HepatitisB2': ('hepatitisB.fasta', 13),
        'InfluenzaA': ('influenzaA.fasta', 56)
    }

    masterfile_name, N = master_files[flag_virus]

    formatted_flag = normalize_flag_virus(flag_virus.replace('_', ''))
    directory = os.path.join(BASE_DATASET_DIR, f"{N}classes_{formatted_flag}_dataset")
    masterpath = os.path.join(BASE_DATASET_DIR, masterfile_name)

    os.makedirs(directory, exist_ok=True)

    with open(masterpath, "r") as f:
        content = f.readlines()
    y_data = []
    for row in content:
        if ">" in row:
            if flag_virus in ['HIV1', 'HIV2']:
                y_class = row.split('.')[0].split('>')[1]
            elif flag_virus == 'Dengue':
                y_class = row.split('|')[3].split(" ")[1].strip()
            elif flag_virus == 'HepatitisC':
                y_class = row.split('.')[0].split('>')[1]
                valid, y_class = classes_HepatitisC(y_class)
                if not valid or y_class == '_':
                    continue
            elif flag_virus == 'HepatitisB1':
                if 'recombinant' in row:
                    continue
                y_class = row.split('|')[2].split('Hepatitis B Virus ')[1].split('.')[0].split('genotype')[1].split(" ")[1]
            elif flag_virus == 'HepatitisB2':
                if 'recombinant' in row:
                    y_class = row.split('|')[2].split('Hepatitis B Virus ')[1].split('.')[0].split(' (')[0]
                else:
                    y_class = row.split('|')[2].split('Hepatitis B Virus')[1].split('.')[0].split('genotype')[1].split(" ")[1]
            elif flag_virus == 'InfluenzaA':
                y_class = row.split(' ')[0].replace('>', "").replace('?',"")
                if y_class.lower() == 'mixed':
                    continue
            y_data.append(y_class)

    cnt = Counter(y_data)
    Classes_dict = dict(sorted(cnt.items(), key=operator.itemgetter(1), reverse=True))

    ClassesN_dict = dict(list(Classes_dict.items())[:N])
    n_samples = sum(ClassesN_dict.values())
    avg = round(n_samples / N)

    if flag_virus in ['HIV1', 'HIV2']:
        prerocessingDataset_HIV(masterpath, directory, ClassesN_dict, N, formatted_flag)
    elif flag_virus == 'Dengue':
        prerocessingDataset_Dengue(masterpath, directory, ClassesN_dict, N, formatted_flag)
    elif flag_virus == 'HepatitisC':
        prerocessingDataset_HepatitisC(masterpath, directory, ClassesN_dict, N, formatted_flag)
    elif flag_virus in ['HepatitisB1', 'HepatitisB2']:
        prerocessingDataset_HepatitisB(masterpath, directory, ClassesN_dict, N, formatted_flag)
    elif flag_virus == 'InfluenzaA':
        prerocessingDataset_InfluenzaA(avg, masterpath, directory, ClassesN_dict, N, formatted_flag)
