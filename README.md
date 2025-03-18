# coloredFCGR

ColoredFCGR repository presents the implementation of Colour Chaos Game Representation (CCGR) and the WalkIm extension for Deep Convolutional Neural Networks (DCNNs). WalkIm: Compact image-based encoding for high-performance classification of biological sequences using simple tuning-free CNNs
Akbari Rokn Abadi S, Mohammadi A, Koohi S (2022) WalkIm: Compact image-based encoding for high-performance classification of biological sequences using simple tuning-free CNNs. PLOS ONE 17(4): e0267106. htgtps://doi.org/10.1371/journal.pone.0267106

**Software configuration:**

The code developed in “coloredFCGR” is in Python programming language version 3.9.

To compile the software correctly, installation of the following python libraries is required: pandas, matplotlib, scikit-learn, keras, and tensorflow.

Also required are the main python packages for scientific computing (numpy), modules for implementing container data types (collections), asynchronous function execution (ProcessPoolExecutor) and image module (Pillow).

It is recommended to use in binary package manager such as conda, anconda or miniconda to create a system-level environment independent of the machine's operating system.

**File Hierarchy**

coloredFCGR
|_ CODE AND EXPERIMENTS
|_ _ WalkIm
|_ _ CCGR

The WalkIm directory contains the encoder unit (in grayscale and RGB), and the Simple CNN and Complex CNN classification models as implemented in https://github.com/SAkbari93/WalkIm. coloredFCGR adds classification with AlexNet and ResNet50 networks to the WalkIm implementation.
The CCGR directory contains the encoding unit with the CCGR decoding strategy and the classification models with AlexNet and ResNet50.

**The developed software contains:**

The Encoder Unit, which decodes a virolgical sequences in image through Colour Chaos Game Representation (CCGR) algorithm.
A classification unit, consisting of deep convolutional neural networks (i.e. AlexNet and ResNet50), which predicts the class to which a virological CCGR image belongs. 

**Software compilation**

Colour Chaos Game Representation used jellyfish (resource: https://github.com/gmarcais/Jellyfish) multi-threader k-mers counter.
In order to run jellyfish, it is necessary to make the binary file executable: 

```chmod +x CODE\ AND\ EXPERIMENTS/CGR-pcmer/jellyfish/jellyfish-binar```

Set in the `Encoder-VIRUES` file (path `/coloredFCGR/CODE AND EXPERIMENTS/CCGR/Encoder-VIRUSES`) jellyfish=True