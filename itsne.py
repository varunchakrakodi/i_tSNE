print("This program is used for generating a interactive tSNE analysis of sequences")
print("Code Compiled by: Varun CN")
print("Require a multifasta file that contains all the sequences. Requires a .csv file that contains information on SequenceID, Country, Year.")

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import plotly.graph_objects as go

# Function to calculate sequence similarities
def calculate_sequence_similarities(seqs, max_seq_length):
    num_seqs = len(seqs)
    padded_seqs = []

    # Pad sequences with N (ambiguous nucleotide) to have the same length
    for seq_record in seqs:
        sequence = str(seq_record.seq)
        padding_length = max_seq_length - len(sequence)
        padded_sequence = sequence + "N" * padding_length
        padded_seqs.append(padded_sequence)

    similarities = np.zeros((num_seqs, num_seqs))

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq_i = padded_seqs[i]
            seq_j = padded_seqs[j]
            similarity = sum(1 for a, b in zip(seq_i, seq_j) if a == b) / len(seq_i)
            similarities[i, j] = similarities[j, i] = similarity

    return similarities

# Run the Main function
def main(fasta_file_path, metadata_file_path):
    # Read the FASTA file and store sequences
    sequences = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Find the maximum sequence length for padding
    max_seq_length = max(len(seq.seq) for seq in sequences)

    # Calculate sequence similarities with padding
    similarities = calculate_sequence_similarities(sequences, max_seq_length)

    # Perform t-SNE analysis
    tsne = TSNE(n_components=2, random_state=42)
    tsne_result = tsne.fit_transform(similarities)

    # Determine the optimal number of clusters using the Elbow method
    distortions = []
    max_clusters = 10  # Set a reasonable upper limit for the number of clusters
    for i in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=i, random_state=42)
        kmeans.fit(tsne_result)
        distortions.append(kmeans.inertia_)

    # Plot the Elbow curve
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(1, max_clusters + 1)), y=distortions, mode='lines+markers'))
    fig.update_layout(title='Elbow Method for Optimal Number of Clusters',
                      xaxis_title='Number of Clusters',
                      yaxis_title='Distortion',
                      showlegend=False)
    fig.show()

    # Based on the Elbow method, enter desired number of clusters
    optimal_clusters = int(input("Optimal number of clusters as per the graph: "))

    # Perform KMeans clustering with the optimal number of clusters
    kmeans = KMeans(n_clusters=optimal_clusters, random_state=42)
    clusters = kmeans.fit_predict(tsne_result)

    # Read metadata from CSV file
    metadata_df = pd.read_csv(metadata_file_path)

    # Merge metadata with cluster information
    cluster_data = {'Sequence ID': [seq.id for seq in sequences], 'Cluster': clusters}
    cluster_df = pd.DataFrame(cluster_data)
    merged_df = pd.merge(cluster_df, metadata_df, on='Sequence ID')

    # Save the merged information to a CSV file
    merged_df.to_csv('cluster_information.csv', index=False)

    # Prompt user to input title for the plot
    plot_title = input("Enter the title for the plot: ")

    # Plot the t-SNE results with different colors for each cluster
    fig = go.Figure()

    for i in range(optimal_clusters):
        indices = np.where(clusters == i)[0]  # Get the indices where the cluster is i
        x_values = tsne_result[indices, 0]    # Extract x values for the cluster
        y_values = tsne_result[indices, 1]    # Extract y values for the cluster
        cluster_metadata = merged_df.loc[merged_df['Cluster'] == i]  # Get metadata for the cluster
        text_values = [f'Sequence ID: {row["Sequence ID"]}<br>Country: {row["Country"]}<br>Year: {row["Year"]}' for idx, row in cluster_metadata.iterrows()]
        fig.add_trace(go.Scatter(x=x_values, y=y_values, mode='markers',
                                 marker=dict(color=i, size=10),
                                 text=text_values,
                                 name=f'Cluster {i}'))

    fig.update_layout(title=plot_title,
                      xaxis_title="Component 1",
                      yaxis_title="Component 2",
                      showlegend=True)

    fig.show()

fasta_file_path = input("Provide path to multifasta file: ")
metadata_file_path = input("Provide path to metadata CSV file: ")

main(fasta_file_path, metadata_file_path)
