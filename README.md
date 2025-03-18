# coloredFCGR

This repository presents the implementation of Colour Chaos Game Representation (CCGR) and the WalkIm extension for Deep Convolutional Neural Networks (DCNNs). WalkIm: Compact image-based encoding for high-performance classification of biological sequences using simple tuning-free CNNs
Akbari Rokn Abadi S, Mohammadi A, Koohi S (2022) WalkIm: Compact image-based encoding for high-performance classification of biological sequences using simple tuning-free CNNs. PLOS ONE 17(4): e0267106. htgtps://doi.org/10.1371/journal.pone.0267106, https://github.com/SAkbari93/WalkIm

**Software configuration:**

The code developed in “coloredFCGR” is in Python programming language version 3.9.

To compile the software correctly, installation of the following python libraries is required: pandas, matplotlib, scikit-learn, keras, and tensorflow.

Also required are the main python packages for scientific computing (numpy), modules for implementing container data types (collections), asynchronous function execution (ProcessPoolExecutor) and image module (Pillow).

It is recommended to use in binary package manager such as conda, anconda or miniconda to create a system-level environment independent of the machine's operating system.


**The developed software contains:**

The Encoder Unit, which decodes a virolgical sequences in image through Colour Chaos Game Representation (CCGR) algorithm.
A classification unit, consisting of deep convolutional neural networks (i.e. AlexNet and ResNet50), which predicts the class to which a virological CCGR image belongs. 

**Software compilation**

In order to run jellyfish, multi-threader k-mers counter, it is necessary to make the binary file executable: `chmod +x CODE\ AND\ EXPERIMENTS/CGR-pcmer/jellyfish/jellyfish-binar`. Set in the `Encoder-VIRUES` file (path `/coloredFCGR/CODE AND EXPERIMENTS/CGR-pcmer/Encoder-VIRUSES`) settare jellyfish=True