import sys
sys.path.insert(1, 'dev/src/')

from VIRUSES.CCGRlib.module import * 
from VIRUSES.CCGRlib.cgr import CGR
from VIRUSES.CCGRlib.fcgr import FCGR, FCGR_RGB
from VIRUSES.CCGRlib.fcgr_pcmer_rgb import FCGR_PCMER_RGB
from VIRUSES.CCGRlib.functions_CGR_PCMER import parse_sequence

# --- Argparse CLI interface ---
parser = argparse.ArgumentParser(
    description='Encode virus genome sequences into RGB images using CCGR coloring schemas.'
)

parser.add_argument('--dataset', type=str,
    default='dev/DATASET/7classes_coronaviruses dataset',
    help='Path to the dataset folder.')

parser.add_argument('--out', type=str,
    default='dev/src/VIRUSES/CCGR_ENCODER',
    help='Output directory for encoded images.')

parser.add_argument('--kmer', type=int,
    default=6,
    help='K-mer size for FCGR encoding.')

parser.add_argument('--encoding', type=str,
    default='pcCCGR',
    choices=['kCCGR', 'pcCCGR'],
    help='Encoding color type to use.')

parser.add_argument('--threshold', type=float, nargs='+',
    default=[0, 0.5, 1],
    help='List of float thresholds T in [0, 1], used for frequency and pc-mer filtering')

parser.add_argument('--jellyfish', action='store_true',
    help='Enable jellyfish k-mer counting (if available).')

args = parser.parse_args()

# --- Start process ---
if __name__ == '__main__':
    start_time = time.time()

    jellyfish = args.jellyfish
    if not jellyfish:
        print('jellyfish multi-threader k-mers counter active')
    else:
        print('jellyfish multi-threader k-mers counter inactive')

    MAIN_FOLDER = 'dev'
    out_directory = args.out
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)

    folder_Dataset = args.dataset
    list_folder_VIRUSES = [folder_Dataset]

    for folder_viruses in list_folder_VIRUSES:
        mypath = folder_viruses
        print('Dataset path:', mypath)

        # Detect virus class based on path name
        if 'hiv' in mypath and '12' in mypath:
            viruses = 'HIV1'
        elif 'hiv' in mypath and '37' in mypath:
            viruses = 'HIV2'
        elif 'dengue' in mypath:
            viruses = 'Dengue'
        elif 'coronaviruses' in mypath:
            viruses = 'Coronaviruses'
        elif 'hepatitisC' in mypath:
            viruses = 'HepatitisC'
        elif 'hepatitisB' in mypath and '8' in mypath:
            viruses = 'HepatitisB1'
        elif 'hepatitisB' in mypath and '13' in mypath:
            viruses = 'HepatitisB2'
        elif 'influenzaA' in mypath:
            viruses = 'InfluenzaA'
            

        w, h = 1024, 1024
        my_list = os.listdir(mypath)

        for i in range(len(my_list)):
            Subpath = os.path.join(mypath, my_list[i])
            onlyfiles = [f for f in listdir(Subpath) if isfile(join(Subpath, f))]

            for s in range(len(onlyfiles)):
                data = np.zeros((h, w), dtype=np.uint8)

                fileFASTA = os.path.join(Subpath, onlyfiles[s])
                genome_viruses = parse_sequence(fileFASTA)
                title = os.path.splitext(onlyfiles[s])[0].replace('.', '-')

                directory_png = os.path.join(out_directory, viruses)
                if not os.path.exists(directory_png):
                    os.makedirs(directory_png)
                print('Saving to:', directory_png)

                FCGR_VIRUSES = FCGR_PCMER_RGB(2 ** (2 * args.kmer),"",0, [],{}).build_fcgr(
                    genome_viruses,
                    args.kmer,
                    args.encoding,
                    args.threshold,
                    jellyfish,
                    fileFASTA,
                    title,
                    directory_png,
                    my_list[i]
                )

                end_time = time.time()
                print('TIME:', end_time - start_time)
