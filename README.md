# CCGR

The CCGR directory contains the implementation of the Color Chaos Game Representation (CCGR) approach. It provides an Econder unit to transform virological sequences into CCGR images and a Network unit consisting of the AlexNet and ResNet50 networks. 


## Third-party software and libraries:

Color Chaos Game Representation used Jellyfish (https://github.com/gmarcais/Jellyfish) multi-threader k-mers counter.
In the experiments, we compared CCGR with WalkIm (https://github.com/SAkbari93/WalkIm).

## Software configuration:

The code in **CCGR** is developed in Python **3.9**.

To run the software, the following Python libraries are required:

`pandas`, `matplotlib`, `scikit-learn`, `keras`, `tensorflow`

Additional packages: `numpy`, `collections`, `concurrent.futures` (e.g. `ProcessPoolExecutor`), and `Pillow` for image processing.

We recommend using a binary package manager such as **Conda**, Anaconda, or Miniconda to create an isolated environment, independent of your operating system


## Installation:

**Clone the Repository**

Download or clone this repository:
```
git clone https://github.com/esignor/CCGR.git
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

These requirements align with the TensorFlow version specified in the environment and requirements files (e.g., TensorFlow 2.11). Meeting these specifications enables GPU acceleration for optimal performance.

*CPU Support:* The software can run on CPU-only systems without a GPU. Recommended hardware specifications for CPU execution include:

- Modern x86_64 compatible CPU
- At least 8 GB of RAM (more recommended for large datasets)
- Python 3.9

While CPU execution is supported and suitable for testing or smaller datasets, GPU acceleration is strongly recommended for faster processing and handling of large genomic datasets.

For our tests on the CCGR software, we used a server cluster equipped with Nvidia A40 GPUs, 1.5 TB of RAM, and Intel Xeon Platinum 8260 CPUs (2.40/3.90 GHz), with a total of 6 TB of RAM available for running the experiments.


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
On systems with limited memory, Conda may crash during installation. In that case, we recommend using the pip-based installation below.

If using a Conda environment, always make sure that you are actually working within the intended environment (e.g., ccgr_env). Some useful commands to verify this include checking that the pip and python executables are pointing to the Conda environment and not the system installation:

```
which python
```
Expected correct terminal output:
```
/home/users/anaconda3/envs/ccgr_env/bin/python
```


🟦 **Option 2 – Install via pip (lightweight alternative)**

If you prefer using `pip` instead of Conda, you can install the required dependencies using the `requirements.txt` file.

This project requires Python 3.9 to run correctly. On most Linux distributions, you can install Python 3.9 using your package manager or by compiling from source. For detailed instructions, please refer to the official Python 3.9 installation guide for Linux:
https://docs.python.org/3.9/using/unix.html. For Windows users, please follow the installation guide here:
https://www.python.org/downloads/release/python-390/. 

1. Create a virtual environment (recommended):

```
python3.9 -m venv ccgr_env
source ccgr_env/bin/activate  # On Windows: ccgr_env\Scripts\activate
```
In casi di problemi o di incompatibilita' e possibile anche creare un ambiente Conda con python3.9 installato e procedere con il punto (2)
2. Install the required packages:


```
pip install -r requirements.txt
```

✅ You’re now ready to run the CCGR scripts.


## Datasets:


## File Hierarchy:

```
CCGR/
├── environment.yaml              # Conda environment definition
├── README.md                     # Project documentation
├── requirements.txt              # pip requirements (alternative to conda)
├── dev/
│   ├── _lib/
│   │   └── jellyfish/            # Contains the jellyfish binary
│   └── _src/
│       └── VIRUSES/
│           ├── CCGR_ENCODER/     # Virological images decoded in CCGR format
│           ├── CCGRlib/          # Implementation of the CCGR library
│           ├── Encoder-VIRUSES.py  # CCGR Encoder unit implementation
│           ├── AlexNet*.py         # CCGR Network unit – AlexNet variants
│           └── ResNet*.py          # CCGR Network unit – ResNet variants
```

`lib/` contains the jellyfish binary file.

`src/` contains the CCGR approach implementation, which includes an encoder unit and a network unit

`src/VIRUSES/CCGRlib/` is the CCGR library implementation.

`src/VIRUSES/CCGR_ENCODER/` stores virological images decoded using the CCGR method.

`Encoder-VIRUSES.py` implements the CCGR Encoder unit.

`AlexNet*.py` and `ResNet*.py` implement the CCGR Network unit, with different architecture variants.



## Software Description:
The developed software includes two main components:

- **Encoder Unit:** responsible for converting virological sequences into images using the Color Chaos Game Representation (CCGR) algorithm.

- **Network Unit:** a deep learning module based on Convolutional Neural Networks (e.g., AlexNet and ResNet50) that classifies the generated CCGR images into their corresponding virological classes.



## Software compilation:
To generate CCGR images, you need to run the `Encoder-VIRUSES` script located at `CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`. This script allows you to select the dataset (from the options described in the *Datasets section*), specify the *k*-mer size, and choose the CCGR image coloring scheme (i.e., kCCGR or pcCCGR). The script will produce CCGR images of the selected dataset using float thresholds between 0 and 1.

**Running the Encoder Script**

To execute the script, navigate to the CCGR directory and run the following command in your terminal to see the available parameters and options:

```
 python dev/src/VIRUSES/Encoder-VIRUSES.py --help
```

**Optional Input Arguments** 
- *--dataset*: path to the dataset directory you want to encode using the CCGR approach.
- *--out*: path to the output directory where the encoded CCGR images will be saved.
- *--kmer*: the size of the *k*-mers to use for FCGR encoding.
- *--encoding*: the coloring scheme to apply in the image encoding (kCCGR or pcCCGR).
- *--threshold*: threshold parameter T (float between 0 and 1) that defines color assignment based on frequency and/or structural components in the CCGR image.
- *--jellyfish*: enable Jellyfish as the k-mer counting tool; if this flag is not set, the internal k-mer counter implemented within the CCGR software will be used.



**Example Usage**
The following command runs the CCGR Encoder with the default coronavirus dataset, a k-mer size of 5, the pcCCGR coloring scheme, a threshold of 1 (using only the structural component), and Jellyfish enabled:

```
python dev/src/VIRUSES/Encoder-VIRUSES.py --kmer 5 --encoding pcCCGR --threshold 1 --jellyfish
```
The expected output is a full decoding of the coronavirus dataset into CCGR images using the pcCCGR color scheme with a threshold of 1.



To test the classification networks, it is necessary to run the AlexNet and ResNet50 models. CCGR directory provides a model for each possible approach implemented and tested in CCGR (kCCGR and pcCCGR with T=[0, 0.5, 1]). For each model, it is possible to set not only the coloring schema approach but also the dataset (from one of the options defined within the model) and the size of *k*-mers. 

We tested *k*-mers for sizes 4 to 10.



In order to run Jellyfish, it is necessary to make the binary file executable. In `CCGR` directory:

```chmod +x dev/lib/jellyfish/jellyfish-binary```

Set in the `Encoder-VIRUSES` file (path `CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`) `jellyfish=True`

At the end of a model run, the trained model is saved in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Models` path and the results of the training in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Results` path.