# Platform for compound prioritization based on transcription factor activity
A platform for drug repurposing for MYC by prioritizing  compounds based on transcription factor activity.
The data and scripts are available in CMAP folder. 
![sbhd_figurev2](https://user-images.githubusercontent.com/48244638/57519061-272fbb80-7323-11e9-9760-8e6cc0d4c82b.png)

# DeepSIBA: Chemical Structure-based Inference of Biological Alterations
### Christos Fotis<sup>1(+)</sup>, Nikolaos Meimetis<sup>1+</sup>, Antonios Sardis<sup>1</sup>, Leonidas G.Alexopoulos<sup>1,2(*)</sup>
 #### 1. BioSys Lab, National Technical University of Athens, Athens, Greece.
#### 2. ProtATonce Ltd, Athens, Greece.

(+)Equal contributions

(*)Correspondence to: leo@mail.ntua.gr

Implementation of the paper:
> DeepSIBA: Chemical Structure-based Inference of Biological Alterations <br>
> Nikos and Chris et. al, 2020

## Abstract
Predicting whether a chemical structure shares a desired biological effect can have a significant impact for in-silico compound screening in early drug discovery.  In this study, we developed a deep learning model where compound structures are represented as graphs and then linked to their biological footprint. To make this complex problem computationally tractable, compound differences were mapped to biological effect alterations using Siamese Graph Convolutional Neural Networks. The proposed model was able to learn new representations from chemical structures and identify structurally dissimilar compounds that affect similar biological processes with high precision. Additionally, by utilizing deep ensembles to estimate uncertainty, we were able to provide more reliable and accurate predictions for chemical structures that are very different from the ones used to the train the models. Finally, we present a novel inference approach, where the trained models are used to provide an estimate of a compound’s effect on signaling pathways, using only its chemical structure.

![sbhd_figurev2](https://user-images.githubusercontent.com/48244638/57519061-272fbb80-7323-11e9-9760-8e6cc0d4c82b.png)

## Requirements

## Clone and initialize
```bash
# clone the source code on your directory
$ git clone https://github.com/jhyuklee/ReSimNet.git
$ cd ReSimNet

# make folder to save and load your data
$ cd tasks
$ mkdir -p data

# make folder to save and load your model
cd ../../..
$ mkdir -p results
```

## Data
Here you can fing all the data produced and used in this study.
The data and models folders must replace the corresponding empty folders existing on Github, in order to recreate the study. 
If there are already data and files in these folders for another study, just the contents of the downloaded folder must be copied in the data folder.
 
1.[data](http://google.com)
2.[models](http://google.com)



