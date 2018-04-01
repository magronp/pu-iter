# PU-Iter: model-based phase recovery algorithm for audio source separation

Here, you will find the code related to PU-Iter algorithm for source separation.

If you use any of the things existing in this repository, please cite the [corresponding paper](https://arxiv.org/abs/1608.01953). 

You can also find an online demo with sound examples related to this work on the [companion website](http://www.cs.tut.fi/~magron/demos/demo_PUITER.html).


## How to use

### Dataset set-up

To reproduce the experiments conducted in our paper, you will need to download the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de) and to place its content in the `dataset/DSD100` folder.

If you use this dataset, you will end up with the proper directory structure and file names, as used in the `functions/get_data_DSD` function.

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match the one from DSD100;
- or you modify the file reading function `get_data_DSD` to suit your needs.

We also provide piano notes from the [MAPS database](http://www.tsi.telecom-paristech.fr/aao/en/2010/07/08/maps-database-a-piano-database-for-multipitch-estimation-and-automatic-transcription-of-music/), as well as the corresponding loading function `get_data_MAPS_notes`. This is intended to reproduce the results of the article on piano sounds, which corresponds to the script `plot_res_piano`.


### Phase recovery

The experiments conducted in the paper rely on two phase recovery algorithms. The corresponding functions can be found in the `functions` folder, and are named `pu_iter.m` and `wiener_filters.m`. You can use those functions on any song you'd like, provided the STFT of the mixture and magnitude/variances estimates.

The script to reproduce the experiments are placed in the `scripts` folder. They will notably record audio files in the `audio_files` folder, and some metrics (SDR, SIR and SAR) in the `metrics` folder.
