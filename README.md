# CCGR

The CCGR directory contains the implementation of the **Color Chaos Game Representation (CCGR)** approach. It provides an Econder unit to transform virological sequences into CCGR images and a Network unit consisting of the AlexNet and ResNet50 networks. 

This README describes installation, dataset preparation, software execution, and output interpretation.


## Third-party software and libraries

Color Chaos Game Representation used **Jellyfish** (https://github.com/gmarcais/Jellyfish) multi-threader k-mers counter.
In the experiments, we compared CCGR with **WalkIm** (https://github.com/SAkbari93/WalkIm).

## Software configuration

The code in **CCGR** is developed in Python **3.9**.

To run the software, the following Python libraries are required: `pandas`, `matplotlib`, `scikit-learn`, `keras`, `tensorflow`

Additional packages: `numpy`, `collections`, `concurrent.futures` (e.g. `ProcessPoolExecutor`), and `Pillow` for image processing.

We recommend using a binary package manager such as **Conda**, Anaconda, or Miniconda to create an isolated environment, independent of your operating system


## Installation

**Clone the Repository**

Download or clone this repository:

To view or download the code anonymously, please use this link:
https://anonymous.4open.science/r/CCGR-22F2

(If needed, click the button "Download Repository")

**Important:** In the commands and discussions related to this repository, `CCGR-22F2` will be abbreviated as `CCGR`.


```
cd CCGR
```

You can install CCGR with either CPU or GPU support, depending on your system and available hardware.
From TensorFlow 2.0 onwards, the tensorflow package works on both CPU and GPU.

**System Compatibility and Hardware Requirements**

This software supports both GPU-accelerated and CPU-only execution using TensorFlow.

*GPU Support:* For GPU acceleration, the following hardware and software requirements must be met:

- NVIDIA GPU (compatible with CUDA)
- CUDA Toolkit version 11.8 or higher
- cuDNN version 8.6 or higher
- NVIDIA driver version 450 or higher

These requirements align with the TensorFlow version specified in the *environment.yaml* and *requirements.txt* files (e.g., TensorFlow 2.11). Meeting these specifications enables GPU acceleration for optimal performance.

*CPU Support:* The software can run on CPU-only systems without a GPU. Recommended hardware specifications for CPU execution include:

- Modern x86_64 compatible CPU
- At least 8 GB of RAM (more recommended for large datasets)
- Python 3.9

While CPU execution is supported and suitable for testing or smaller datasets, GPU acceleration is strongly recommended for faster processing and handling of large genomic datasets.

*Note:* 
For our tests on the **CCGR** software, we used a server cluster equipped with **Nvidia A40 GPUs, 1.5 TB of RAM** (to Network units), and **Intel Xeon Platinum 8260 CPUs (2.40/3.90 GHz) with a total of 6 TB of RAM available** (to Encoder unit) for running the experiments.


🟩 **Option 1 – Install via Conda (recommended for full reproducibility)**


1. Create the environment:

```
conda env create -f environment.yaml
```

2. Activate the environment:

```
conda activate ccgr_env
```

✅ You’re now ready to run the CCGR scripts.

*Note:* 
On systems with limited memory, Conda may crash during installation. In that case, we recommend using the **pip-based installation** below.

If using a Conda environment, always make sure that you are actually working within the intended environment (e.g., *ccgr_env*). Some useful commands to verify this include checking that the pip and python executables are pointing to the Conda environment and not the system installation:

```
which python    # alternative: which pip
```
Expected (or similar) correct terminal output:
```
/home/users/anaconda3/envs/ccgr_env/bin/python
```


🟦 **Option 2 – Install via pip (lightweight alternative)**

If you prefer using `pip` instead of Conda, you can install the required dependencies using the `requirements.txt` file.

This project requires Python 3.9 to run correctly. On most Linux distributions, you can install Python 3.9 using your package manager or by compiling from source. For detailed instructions, please refer to the official Python 3.9 installation guide for *Linux*:
https://docs.python.org/3.9/using/unix.html. 
For *Windows* users, please follow the installation guide here:
https://www.python.org/downloads/release/python-390/. 

1. Create a virtual environment (recommended):

```
python3.9 -m venv ccgr_env
source ccgr_env/bin/activate  # On Windows: ccgr_env\Scripts\activate
```
In case of issues or incompatibilities to (1), it is also possible to create a Conda environment with Python 3.9 installed and proceed with step (2).

2. Install the required packages:


```
pip install -r requirements.txt
```

✅ You’re now ready to run the CCGR scripts.



## File Hierarchy

```
CCGR/
├── environment.yaml              # Conda environment definition
├── README.md                     # Project documentation
├── requirements.txt              # pip requirements (alternative to conda)
├── dev/
|   ├── _DATASET/                       # Contains CCGR-compatible datasets
│   │   └── DatasetsGuidelines.md       # Datasets documentation
|   ├── VirusPreprocessingDatasets.py   # Script to convert virological .fasta files into CCGR-compatible datasets
│   ├── _lib/
│   │   └── jellyfish/            # Contains the jellyfish binary
│   └── _src/
│       └── VIRUSES/
│           ├── CCGR_ENCODER/     # Virological images decoded in CCGR format
│           ├── CCGRlib/          # Implementation of the CCGR library
│           ├── Encoder-VIRUSES.py  # CCGR Encoder unit implementation
│           ├── NetworkUnit_AlexNet.py         # CCGR Network unit – AlexNet architecture
│           └── Networkunit_ResNet50.py        # CCGR Network unit – ResNet50 architecture
```

- `CCGR/`: root directory of the Color Chaos Game Representation (CCGR) project.

- `README`: main documentation file for the project (current file).

- `environment.yaml`: conda environment definition file for dependency setup.

- `requirements.txt`: (alternative to conda) used to install dependencies via pip.

- `dev/`: contains the implementation of the CCGR software, including preprocessing scripts, libraries, and model architecture definitions.

- `dev/DATASET`: directory containing the datasets compatible with CCGR, along with the necessary documentation (*DatasetsGuidelines.md*) to obtain them.

- `dev/VirusPreprocessingDatasets.py`: script that converts virological genome sequences in .fasta format into datasets compatible with the CCGR pipeline.

- `dev/lib/`: directory containing the Jellyfish binary used for efficient k-mer counting.

- `dev/src/`: contains the CCGR approach implementation, which includes an Encoder unit and a Network unit.

- `src/VIRUSES/CCGRlib/`: core library for the CCGR image encoding method.

- `src/VIRUSES/CCGR_ENCODER/`: stores output images generated by applying CCGR encoding to input virological sequences.

- `Encoder-VIRUSES.py`: implements the CCGR Encoder unit.

- `NetworkUnit_AlexNet.py` and `NetworkUnit_ResNet50.py`: implement the CCGR Network unit, with different architecture variants.

## Datasets

The data processed by the CCGR software belong to the virus category. It was natively developed to work on 8 datasets (belonging to 6 different virus families), namely: Coronaviruses, HIV-1, Dengue, Hepatitis C, Hepatitis B, and Influenza A.


The table below lists: **Virus Name**, which is the name of the virological dataset recognized by the CCGR software; **FASTA Name**, the name of the Fasta file containing the genomic sequences; and **Dataset Directory Name**, the name of the dataset directory after preprocessing, made compatible with CCGR.
Details regarding the sources of the Fasta files, the data extraction queries, and the preprocessing scripts are provided immediately following the table.


| Virus Name     | FASTA Name     | Dataset Directory Name             |
|----------------|----------------|------------------------------------|
| Coronaviruses  | None           | 7classes_coronaviruses_dataset     |
| HIV1           | hiv-db         | 12classes_hiv1_dataset             |
| HIV2           | hiv-db         | 37classes_hiv2_dataset             |
| Dengue         | dengue         | 4classes_dengue_dataset            |
| HepatitisC     | HCV            | 6classes_hepatitisC_dataset        |
| HepatitisB1    | hepatitisB     | 8classes_hepatitisB1_dataset       |
| HepatitisB2    | hepatitisB     | 13classes_hepatitisB2_dataset      |
| InfluenzaA     | influenzaA     | 56classes_influenzaA_dataset       | 


* **Coronaviruses:**

1. The dataset can be downloaded from the following address: [https://github.com/SAkbari93/WalkIm/tree/main/Data](https://github.com/SAkbari93/WalkIm/tree/main/Data) (authors of WalkIm)

2. Extract the `.raw` file

3. Rename the directory to `7classes_coronaviruses_dataset`

4. Place `7classes_coronaviruses_dataset` inside `CCGR/dev/DATASET`

The name assigned by the CCGR software to this dataset is *Coronaviruses*

* **HIV1:**

1. Go to the LANL Database webpage ([https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html](https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html))

2. Select the following query parameters: *virus:* HIV-1, *genomic region:* complete genome, *subtype:* any subtype, excluding problematic, *other options:* default. Press the `Search` button on the interface

3. Download the *Fasta* file and name it `hiv-db.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HIV1`

The name assigned by the CCGR software to this dataset is *HIV1*

* **HIV2:**

1. Go to the LANL Database webpage ([https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html](https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html))

2. Select the following query parameters: *virus:* HIV-1, *genomic region:* complete genome, *subtype:* any subtype, excluding problematic, *other options:* default. Press the `Search` button on the interface

3. Download the *Fasta* file and name it `hiv-db.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HIV2`

The name assigned by the CCGR software to this dataset is *HIV2*

* **Dengue:**

1. Go to the NCBI Database webpage ([https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi](https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi))

2. Select the following query parameters: *sequence type:* nucleotide, full-length sequences only, collapse identical sequences, *other options:* default. Press the `Show Results` button on the interface

3. Download the *Fasta* file and name it `dengue.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus Dengue`

The name assigned by the CCGR software to this dataset is *Dengue*

* **HepatitisC:**

1. Go to the LANL Database webpage ([https://hcv.lanl.gov/components/sequence/HCV/search/searchi.html](https://hcv.lanl.gov/components/sequence/HCV/search/searchi.html))

2. Select the following query parameters: *sequence information:* any genotype, exclude recombinants, *genomic region:* complete region, *exclude:* problematic sequence, *other options:* default. Press the `Search` button on the interface

3. Download the *Fasta* file and name it `HCV.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisC`

The name assigned by the CCGR software to this dataset is *HepatitisC*

* **HepatitisB1:**

1. Go to the HBVdb Database webpage ([https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0))

2. Download the *Fasta* file (select the options *genotype:* All, *sequence type:* Genomes) and name it `hepatitisB.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

3. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisB1`

The name assigned by the CCGR software to this dataset is *HepatitisB1*

* **HepatitisB2:**

1. Go to the HBVdb Database webpage ([https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0))

2. Download the *Fasta* file (select the options *genotype:* All, *sequence type:* Genomes) and name it `hepatitisB.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

3. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisB2`

The name assigned by the CCGR software to this dataset is *HepatitisB2*

* **InfluenzaA:**

1. Go to the NCBI Database webpage ([https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform))

2. Select the following query parameters: *sequence type:* nucleotide, *type:* A, *sequence length:* full-length only, *collection date:* from 2013-01-01 to 2023-12-31, *additionals filters:* collapse identical sequences, *other options:* default. Press the `Show Results` button on the interface

3. Download the *Fasta* file and name it `influenzaA.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus InfluenzaA`

The name assigned by the CCGR software to this dataset is *InfluenzaA*





## Software Description
The developed software includes two main components:

- **Encoder unit:** responsible for converting virological sequences into images using the Color Chaos Game Representation (CCGR) algorithm.

- **Network unit:** a deep learning module based on Convolutional Neural Networks (e.g., AlexNet and ResNet50) that classifies the generated CCGR images into their corresponding virological classes.

Before proceeding with software compilation (running the Encoder and Network unit scripts), it is **ESSENTIAL** to download the FASTA files from the appropriate databases (see the *Datasets section*) and perform data preprocessing; use the command
```python dev/VirusPreprocessingDatasets.py --help```
for more details.



## Software compilation

### Running the Encoder unit Script
To generate CCGR images, you need to run the `Encoder-VIRUSES` script located at `CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`. This script allows you to select the **dataset** (from the options described in the *Datasets section*), specify the ***k*-mer size**, choose the CCGR image **coloring scheme** (i.e., kCCGR or pcCCGR), and a **threshold float value**. The script will produce CCGR images of the selected dataset using float thresholds between 0 and 1.

To execute the Encoder unit script, navigate to the CCGR directory and run the following command in your terminal to see the available parameters and options:

```
 python dev/src/VIRUSES/Encoder-VIRUSES.py --help
```

#### Optional Input Arguments
- *--virus*: name of the virus to which the CCGR Encoder should be applied (i.e., Coronaviruses, HIV1, HIV2, Dengue, HepatitisC, HepatitisB1, HepatitisB2, InfluenzaA). See *Datasets section*.
- *--kmer*: the size of the *k*-mers to use for FCGR encoding.
- *--encoding*: the coloring scheme to apply in the image encoding (kCCGR or pcCCGR).
- *--threshold*: threshold parameter T (float between 0 and 1) that defines color assignment based on frequency and/or structural components in the CCGR image.
- *--jellyfish*: enable Jellyfish as the k-mer counting tool; if this flag is not set, the internal k-mers counter implemented within the CCGR software will be used.

To use Jellyfish, make sure to make the binary file executable with the command: `chmod +x dev/lib/jellyfish/jellyfish-binary`



#### Example Usage

The following command runs the CCGR Encoder with **coronaviruses dataset, a *k*-mer size of 6, the pcCCGR coloring scheme, a threshold of 1 (using only the structural component), and Jellyfish enabled**:

```
python dev/src/VIRUSES/Encoder-VIRUSES.py --virus Coronaviruses --kmer 6 --encoding pcCCGR --threshold 1 --jellyfish
```
The expected output is a full decoding of the coronaviruses dataset into CCGR images using the pcCCGR color scheme with a threshold of 1.


### Running the Network unit Scripts
 
To test the Network unit, it is necessary to run the `AlexNet` and `ResNet50` models located at `CCGR/dev/src/VIRUSES/NetworkUnit_AlexNet.py` and `CCGR/dev/src/VIRUSES/NetworkUnit_ResNet50.py`. The CCGR directory provides a model for each network architecture implemented and tested in CCGR. For each model, it is possible to configure a variety of parameters, such as the **coloring scheme**, the **CCGR image dataset**, and the **network parameters**.

To execute the script, navigate to the CCGR directory and run the following command in your terminal to see the available parameters and options:

- AlexNet Network unit:
```
 python dev/src/VIRUSES/NetworkUnit_AlexNet.py --help
```

- ResNet50 Network unit:
```
 python dev/src/VIRUSES/NetworkUnit_ResNet50.py --help
```

**Before proceeding with the training of an image set, the images must first be generated by the Encoder unit**. Once generated and ready for training, CCGR image sets are stored in the `CCGR_ENCODER` directory (path: `CCGR/dev/src/VIRUSES/CCGR_ENCODER`).

#### Optional Input Arguments
- *--dataset*: name of the virus dataset to which the CCGR Network unit will be applied (i.e., Coronaviruses, HIV1, HIV2, Dengue, HepatitisC, HepatitisB1, HepatitisB2, InfluenzaA).
- *--type_encoder*: image format for CCGR input (either Grayscale or RGB).
- *--kmer*: the size of *k*-mers used for FCGR in the CCGR image.
- *--threshold*: threshold parameter T (a float between 0 and 1) applied to the CCGR image set.
- *--type_encodingColour*: the coloring scheme applied to the CCGR image set (kCCGR or pcCCGR).
- *--batch_size*: number of training samples used in each iteration of the Network unit model.
- *--epochs*: number of training epochs for the Network unit model.
- *--n_task*: number of tasks to be executed in parallel during model training.


#### Training Technical Note

CCGR images generated using the pcCCGR coloring scheme utilize all three color channels and are therefore in RGB format. As a result, the ***type_encoder* value Grayscale is not supported for pcCCGR images**.

#### Example Usage

The following command starts the training of the a Network unit on pre-generated **CCGR images of the coronaviruses**, using the following image parameters: ***k*-mer size set to 6, *threshold* set to 1, and color encoding scheme** (*type_encodingColour*) set to **pcCCGR**. For training, the chosen ***batch_size* is 30, *epochs* is set to 30, and 2 tasks are executed in parallel**.


- AlexNet Network unit:
```
python dev/src/VIRUSES/NetworkUnit_AlexNet.py --dataset Coronaviruses --type_encoder RGB --kmer 6 --threshold 1 --type_encodingColour pcCCGR --batch_size 30 --epochs 30 --n_task 2
```

- ResNet50 Network unit:
```
python dev/src/VIRUSES/NetworkUnit_ResNet50.py --dataset Coronaviruses --type_encoder RGB --kmer 6 --threshold 1 --type_encodingColour pcCCGR --batch_size 30 --epochs 30 --n_task 2
```

At the end of the training, **the trained model (in .keras format)** is expected to be obtained, along with **plots and summary metrics for the training, validation, and test sets** (see the *Results Training section*).



## Results Training

At the end of a model run, two main directories are generated:

1. **Trained Models**

   **Path:** `CCGR/dev/src/VIRUSES/CCGR [Virus Name] Models`

    Each trained model from the Network unit for the specified *Virus Name* is saved in *.keras* format.

2. **Training Results**

   **Path:** `CCGR/dev/src/VIRUSES/CCGR [Virus Name] Results`

    For each training session, the following outputs are saved:

    a) **Performance summary file**

    *Contains:*

    - Confusion matrix

    - Classification report

    - Training time

    - Test and validation accuracy

    *Filename format:*
    `[type_encoder]results [Network Unit] CCGR([kmer threshold type_encodingColour]) batchsize [batch_size] epoch [epochs].txt`


    b) **Training and validation loss plot**

    *Filename format:*
    `Training-Validation Loss [Network Unit] CCGR([kmer threshold type_encodingColour]) batchsize [batch_size] epoch [epochs].png`


    **Training and validation accuracy plot**

    c) *Filename format:*
    `Training-Validation Accuracy [Network Unit] CCGR([kmer threshold type_encodingColour]) batchsize [batch_size] epoch [epochs].png`

