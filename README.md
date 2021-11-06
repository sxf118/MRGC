# MRGC

Multi-view Robust Graph-based clustering for Cancer Subtype Identification

### Method Description

Our method first learns robust latent representations from the raw omics data to alleviate the influences of the experimental and biological noise, where a set of similarity matrices are then adaptively learned based on these new representations. Finally, a global similarity graph is obtained by exploiting the consensus structure from the graphs of each view. As a result, the three parts in our method can reinforce each other in a mutual iterative manner.

### Environment

MRGC was developed in MATLAB 2019b

### Usage

We provided a demo for users. To run this demo, please load the script 'demo.m' into your MATLAB programming environment and click 'run'.

### Dataset

All the cancer datasets used can be downloaded at http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.

### Parameters

There are three parameters in our method, i.e., 'alpha', 'beta', and the dictionary size 'base'. The default value is 0.01, 0.001 and 10, repectively. Users can change their value in 'demo.m' .

### Input and output

Users can change the input file directory and output file directory by changing the 'dataDir' variable and the 'outDir' variable in 'demo.m', respectively.
