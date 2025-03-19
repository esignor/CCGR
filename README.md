# CCGR

The CCGR directory contains the implementation of the Color Chaos Game Representation (CCGR) approach. It provides an econder unit to transform virological sequences into CCGR images and a classifier unit consisting of the AlexNet and ResNet50 networks. 

**Software configuration:**

The code developed in “coloredFCGR” is in Python programming language version 3.9.

To compile the software correctly, installation of the following python libraries is required: pandas, matplotlib, scikit-learn, keras, and tensorflow.

Also required are the main python packages for scientific computing (numpy), modules for implementing container data types (collections), asynchronous function execution (ProcessPoolExecutor) and image module (Pillow).

It is recommended to use in binary package manager such as conda, anconda or miniconda to create a system-level environment independent of the machine's operating system.

**File Hierarchy**

CCGR

|_ dev

|_ _ lib

        |_ jellyfish


|_ _ src/VIRUSES

        |_ CCGRlib

        |_ CCGR_ENCODER.py

        |_ Encoder.py

        |_ AlexNet.py

        |_ ResNet.py

In `lib` is contained the jellyfish binary file

In `src` is the implementation of the CCGR approach, consisting of encoder units and classifier units

In `src/VIRUSES/module` is the implementation of the CCGR library


In `src/VIRUSES/CCGR_ENCODER` are stored the virological images decoded in CCGR

**The developed software contains:**

The Encoder Unit, which decodes a virolgical sequences in image through Colour Chaos Game Representation (CCGR) algorithm.
A classification unit, consisting of deep convolutional neural networks (i.e. AlexNet and ResNet50), which predicts the class to which a virological CCGR image belongs. 

**Software compilation**

Colour Chaos Game Representation used jellyfish (resource: https://github.com/gmarcais/Jellyfish) multi-threader k-mers counter.
In order to run jellyfish, it is necessary to make the binary file executable: 

```chmod +x dev/lib/jellyfish/jellyfish-binar```

Set in the `Encoder-VIRUSES` file (path `/CCGR/src/VIRUSES/Encoder-VIRUSES.py`) `jellyfish=True`