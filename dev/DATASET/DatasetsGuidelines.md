# Datasets Guidelines

The data processed by the CCGR software belong to the virus category. It was natively developed to work on 8 datasets (belonging to 6 different virus families), namely: Coronaviruses, HIV-1, Dengue, Hepatitis C, Hepatitis B, and Influenza A.


The table below lists: *Virus Name*, which is the name of the virological dataset recognized by the CCGR software; *FASTA Name*, the name of the Fasta file containing the genomic sequences; and *Dataset Directory Name*, the name of the dataset directory after preprocessing, made compatible with CCGR.
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


* **Coronavirus:**

1. The dataset can be downloaded from the following address: [https://github.com/SAkbari93/WalkIm/tree/main/Data](https://github.com/SAkbari93/WalkIm/tree/main/Data) (authors of WalkIm)

2. Extract the `.raw` file

3. Rename the directory to `7classes_coronaviruses_dataset`

4. Place `7classes_coronaviruses_dataset` inside `CCGR/dev/DATASET`

The name assigned by the CCGR software to this dataset is *Coronavirus*

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

2. Select the following query parameters: *sequence information:* any genotype, exclude recombinants ,*genomic region:* complete region, *exclude:* problematic sequence, *other options:* default. Press the `Show Results` button on the interface

3. Download the *Fasta* file and name it `HCV.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

4. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisC`

The name assigned by the CCGR software to this dataset is *HepatitisC*

* **HepatitisB1:**

1. Go to the HBVdb Database webpage ([http://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset](http://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset))

2. Download the *Fasta* file and name it `hepatitisB.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

3. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisB1`

The name assigned by the CCGR software to this dataset is *HepatitisB1*

* **HepatitisB2:**

1. Go to the HBVdb Database webpage ([http://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset](http://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset))

2. Download the *Fasta* file and name it `hepatitisB.fasta`. Place the *.fasta* file in `CCGR/dev/DATASET`

3. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus HepatitisB2`

The name assigned by the CCGR software to this dataset is *HepatitisB2*

* **InfluenzaA:**

1. Go to the NCBI Database webpage ([https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform](https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi#mainform))

2. Select the following query parameters: *sequence type:* nucleotide, *type:* A, *sequence length:* full-length only, *collection date:* from 2013-01-01 to 2023-12-31, *additionals filters:* collapse identical sequences *other options:* default. Press the `Show Results` button on the interface

3. From within CCGR, run the command `python dev/VirusPreprocessingDatasets.py --virus InfluenzaA`

The name assigned by the CCGR software to this dataset is *InfluenzaA*