# CCGR

The CCGR directory contains the implementation of the Color Chaos Game Representation (CCGR) approach. It provides an econder unit to transform virological sequences into CCGR images and a classifier unit consisting of the AlexNet and ResNet50 networks. 


**Third-party software and libraries:**

Color Chaos Game Representation used jellyfish (https://github.com/gmarcais/Jellyfish) multi-threader k-mers counter.
In the experiments, we compared CCGR with WalkIm (https://github.com/SAkbari93/WalkIm).


**Software configuration:**

The code developed in “CCGR” is in Python programming language version 3.9.

To compile the software correctly, installation of the following python libraries is required: pandas, matplotlib, scikit-learn, keras, and tensorflow.

Also required are the main python packages for scientific computing (numpy), modules for implementing container data types (collections), asynchronous function execution (ProcessPoolExecutor) and image module (Pillow).

It is recommended to use in binary package manager such as conda, anconda or miniconda to create a system-level environment independent of the machine's operating system.

**File hierarchy:**

```
CCGR

|_ dev

        |_ _ lib

                |_ jellyfish


        |_ _ src/VIRUSES

                |_ CCGRlib

                |_ CCGR_ENCODER

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





**The developed software contains:**

The Encoder unit, which decodes a virolgical sequences in image through Color Chaos Game Representation (CCGR) algorithm.
A Network unit, consisting of Deep Convolutional Neural Networks (i.e. AlexNet and ResNet50), which predicts the class to which a virological CCGR image belongs. 

**Software compilation:**

To generate CCGR images, you need to run the `Encoder-VIRUSES` file (path `/CCGR/src/VIRUSES/Encoder-VIRUSES.py`). You can choose the dataset (from one of the options defined within the encoder), the number of *k*-mers and the decoding approach of the CCGR image (i.e., kCCGR or pcCCGR). CCGR images of the selected dataset with thresholds T=0, T=0.5, and T=1 will be produced.

To test the classification networks, it is necessary to run the AlexNet and ResNet50 models. CCGR directory provides a model for each possible approach implemented and tested in CCGR (kCCGR and pcCCGR with T=[0, 0.5, 1]). For each model, it is possible to set not only the coloring schema approach but also the dataset (from one of the options defined within the model) and the size of *k*-mers. 

We tested k-mers for sizes 4 to 10.



In order to run jellyfish, it is necessary to make the binary file executable. In `CCGR` directory:

```chmod +x dev/lib/jellyfish/jellyfish-binary```

Set in the `Encoder-VIRUSES` file (path `CCGR/dev/src/VIRUSES/Encoder-VIRUSES.py`) `jellyfish=True`

At the end of a model run, the trained model is saved in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Models` path and the results of the training in `CCGR/dev/src/VIRUSES/CCGR [NAME DATASET] Results` path.