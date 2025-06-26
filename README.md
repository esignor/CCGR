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

**System Compatibility and Hardware Requirements:**

This software supports both GPU-accelerated and CPU-only execution using TensorFlow.

*GPU Support:* For GPU acceleration, the following hardware and software requirements must be met:

- NVIDIA GPU (compatible with CUDA)
- CUDA Toolkit version 11.8 or higher
- cuDNN version 8.6 or higher
- NVIDIA driver version 450 or higher

These requirements align with the TensorFlow version specified in the environment and requirements files (e.g., TensorFlow 2.13). Meeting these specifications enables GPU acceleration for optimal performance.

*CPU Support:* The software can run on CPU-only systems without a GPU. Recommended hardware specifications for CPU execution include:

- Modern x86_64 compatible CPU
- At least 8 GB of RAM (more recommended for large datasets)
- Python 3.9 or higher installed

While CPU execution is supported and suitable for testing or smaller datasets, GPU acceleration is strongly recommended for faster processing and handling of large genomic datasets.



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

*Note: On systems with limited memory, Conda may crash during installation. In that case, we recommend using the pip-based installation below.*


🟦 **Option 2 – Install via pip (lightweight alternative)**
If you prefer using `pip` instead of Conda, you can install the required dependencies using the `requirements.txt` file.

1. Create a virtual environment (recommended):

```
python -m venv ccgr_env
source ccgr_env/bin/activate  # On Windows: ccgr_env\Scripts\activate
```

2. Install the required packages:


```
pip install -r requirements.txt
```

✅ You’re now ready to run the CCGR scripts.


## File hierarchy:

```
CCGR/
|_ environment.yaml
|_ README.md
|_ requirements.txt
|_ dev/
        |_ _ lib/
                |_ jellyfish
        |_ _ src/VIRUSES/
                |_ CCGR_ENCODER/
                |_ CCGRlib/
                |_ Encoder-VIRUSES.py
                |_ AlexNet*.py
                |_ ResNet*.py
```

In `lib` is contained the jellyfish binary file

In `src` is the implementation of the CCGR approach, consisting of an encoder unit and a network unit

In `src/VIRUSES/CCGRlib`  is the implementation of the CCGR library

In `src/VIRUSES/CCGR_ENCODER` are stored the virological images decoded in CCGR

`Encoder-VIRUSES.py` is the implementation of the CCGR Encoder unit

`AlexNet*.py` and `ResNet*.py` are the implementation of the CCGR Network unit





## The developed software contains:

The Encoder unit, which decodes a virolgical sequences in image through Color Chaos Game Representation (CCGR) algorithm.
A Network unit, consisting of Deep Convolutional Neural Networks (i.e. AlexNet and ResNet50), which predicts the class to which a virological CCGR image belongs. 

## Software compilation:

To generate CCGR images, you need to run the `Encoder-VIRUSES` file (path `/CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`). You can choose the dataset (from one of the options defined within the encoder), the size of *k*-mers and the coloring schema approach of the CCGR image (i.e., kCCGR or pcCCGR). CCGR images of the selected dataset with thresholds T=0, T=0.5, and T=1 will be produced.


 To execute the `Encoder-VIRUSES.py` script, navigate to the `CCGR` directory and run the following command from the terminal:

```
python dev/src/VIRUSES/Encoder-VIRUSES.py
```

To test the classification networks, it is necessary to run the AlexNet and ResNet50 models. CCGR directory provides a model for each possible approach implemented and tested in CCGR (kCCGR and pcCCGR with T=[0, 0.5, 1]). For each model, it is possible to set not only the coloring schema approach but also the dataset (from one of the options defined within the model) and the size of *k*-mers. 

We tested *k*-mers for sizes 4 to 10.



In order to run Jellyfish, it is necessary to make the binary file executable. In `CCGR` directory:

```chmod +x dev/lib/jellyfish/jellyfish-binary```

Set in the `Encoder-VIRUSES` file (path `CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`) `jellyfish=True`

At the end of a model run, the trained model is saved in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Models` path and the results of the training in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Results` path.