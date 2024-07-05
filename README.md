This Interactive t-SNE analysis is a useful visualisation tool to understand the clustering nature of various sequences.

**Note:**

The analysis works best, if there are no ambiguities. If nucleotide or the input amino acid sequence has ambiguities, it gets integrated to calculations leading to a non reliable plot. Hence if there are sequences that contain multiple ambiguities (For example, X or N's) it would be useful to remove such sequences from the input multifasta file. Else, interpret the analysis with caution and subjectively. Just like any other Clustering algorithm more the number of sequences, better the results.

Dependencies:

It is best to have Anaconda installed in the system (https://www.anaconda.com/download/), which will fullfill all dependencies

Alternatively, following packages should be available in the system
1. Numpy https://numpy.org/
2. Pandas https://pandas.pydata.org/
3. Biopython https://biopython.org/
4. Scikit Learn https://scikit-learn.org/
5. Plotly https://plotly.com/

The Method uses an inbuilt **Elbow Plot** to predict the number of clusters. Once Elbow plot is generated, the script will request user for number of clusters based on idividual judgement (n <=10). Alternatively, you can ignore elbow plot values and independely run a **nbclust program** (https://www.rdocumentation.org/packages/NbClust/versions/3.0.1/topics/NbClust) and input values accordingly.

**Command**

python3 /path/to/itsne.py.

The script will ask the user for step by step input (Provide complete paths), compile and run accordingly

The following input files will be requested for:

1. aligned.fasta (Requires a multifasta file containing all sequences to be analysed that has undergone multiple sequence alignment. It is recommended to use mafft; https://github.com/GSLBiotech/mafft. However, other programs should work equally well !!!)
2. metadata.csv (Requires a .csv file that contain following headers: SequenceID, Country, Year with corresponding information. If details are unknown, enter unknown)

Output files:

An interactive .html file that can be opened with any html browser and a cluster_information.csv that will contain details on SequenceID binning to cluster, which can be used for further analysis.
