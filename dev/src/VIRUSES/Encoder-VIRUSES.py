import sys
sys.path.insert(1, 'dev/src/')

from VIRUSES.CCGRlib.module import * 
from VIRUSES.CCGRlib.cgr import CGR
from VIRUSES.CCGRlib.fcgr import FCGR, FCGR_RGB
from VIRUSES.CCGRlib.fcgr_pcmer_rgb import FCGR_PCMER_RGB
from VIRUSES.CCGRlib.functions_CGR_PCMER import parse_sequence

# python dev/src/VIRUSES/Encoder-VIRUSES.py --virus Coronaviruses --kmer 6 --encoding pcCCGR

# --- Mapping virus name to dataset directory ---
VIRUS_DATASET_PATHS = {
    'Coronaviruses': 'dev/DATASET/7classes_coronaviruses_dataset',
    'HIV1': 'dev/DATASET/12classes_hiv1_dataset',
    'HIV2': 'dev/DATASET/37classes_hiv2_dataset',
    'Dengue': 'dev/DATASET/4classes_dengue_dataset',
    'HepatitisC': 'dev/DATASET/6classes_hepatitisC_dataset',
    'HepatitisB1': 'dev/DATASET/8classes_hepatitisB1_dataset',
    'HepatitisB2': 'dev/DATASET/13classes_hepatitisB2_dataset',
    'InfluenzaA': 'dev/DATASET/56classes_influenzaA_dataset' 
}

# --- Argparse CLI interface ---
parser = argparse.ArgumentParser(
    description='Encode virus genome sequences into RGB images using CCGR coloring schemas'
)

parser.add_argument('--virus', type=str, required=True,
    choices=VIRUS_DATASET_PATHS.keys(),
    help='Name of the virus (i.e., Coronaviruses, HIV1, HIV2, Dengue, HepatitisC, HepatitisB1, HepatitisB2, InfluenzaA)')

parser.add_argument('--kmer', type=int,
    default=6,
    help='K-mer size for FCGR encoding')

parser.add_argument('--encoding', type=str,
    default='pcCCGR',
    choices=['kCCGR', 'pcCCGR'],
    help='Encoding color type to use.')

parser.add_argument('--threshold', type=float, nargs='+',
    default=[0, 0.5, 1],
    help='List of float thresholds T in [0, 1], used for frequency and pc-mer filtering')

parser.add_argument('--jellyfish', action='store_true',
    help='Enable jellyfish k-mer counting (if available)')

args = parser.parse_args()

# --- Start process ---
if __name__ == '__main__':
    start_time = time.time()

    jellyfish = args.jellyfish
    if not jellyfish:
        print('jellyfish multi-threader k-mers counter active')
    else:
        print('jellyfish multi-threader k-mers counter inactive')

     # Get dataset path from mapping
    virus_name = args.virus
    folder_Dataset = VIRUS_DATASET_PATHS[virus_name]

    # Output dir: e.g. CCGR/dev/src/VIRUSES/CCGR_ENCODER/Coronaviruses
    out_directory = os.path.join('dev', 'src', 'VIRUSES', 'CCGR_ENCODER')
    os.makedirs(out_directory, exist_ok=True)

    print(f'Using dataset: {folder_Dataset}')
    print(f'Output will be saved in: {out_directory}')

    my_list = os.listdir(folder_Dataset)

    w, h = 1024, 1024

    for virus_class in my_list:
        Subpath = os.path.join(folder_Dataset, virus_class)
        onlyfiles = [f for f in listdir(Subpath) if isfile(join(Subpath, f))]

        for filename in onlyfiles:
            fileFASTA = os.path.join(Subpath, filename)
            genome_viruses = parse_sequence(fileFASTA)
            title = os.path.splitext(filename)[0].replace('.', '-')

            directory_png = os.path.join(out_directory, virus_name)
            os.makedirs(directory_png, exist_ok=True)
            print('Saving to:', directory_png)

            FCGR_PCMER_RGB(2 ** (2 * args.kmer), "", 0, [], {}).build_fcgr(
                genome_viruses,
                args.kmer,
                args.encoding,
                args.threshold,
                jellyfish,
                fileFASTA,
                title,
                directory_png,
                virus_class
            )

    end_time = time.time()
    print('TOTAL TIME:', end_time - start_time)
