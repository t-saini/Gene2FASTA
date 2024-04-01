# Gene2FASTA
**Author:** Tanvir Saini

Gene2FASTA is a Python script designed to retrieve gene information, DNA sequences, and homologs from the MyGene and Ensembl databases via REST API. It allows users to input a gene name and generates a FASTA file containing the gene's DNA sequence along with the longest Open Reading Frame (ORF) translated into an amino acid sequence. Additionally, it extracts homologous genes from various species and creates a list of unique species with their homologs.

## Functionality Overview

1. **Fetching Gene Information:** Retrieves gene information from the MyGene database using the provided gene name.
2. **Converting to Ensembl ID:** Converts the retrieved gene information to the corresponding Ensembl ID.
3. **Retrieving DNA Sequence:** Fetches the DNA sequence associated with the Ensembl ID from the Ensembl database.
4. **Finding Longest ORF:** Identifies the longest Open Reading Frame (ORF) within the DNA sequence.
5. **Translating DNA to Amino Acids:** Translates the longest ORF into an amino acid sequence.
6. **Finding Homologs:** Searches for homologous genes using the Ensembl ID.
7. **Sorting Homology Results:** Sorts and extracts unique species from the homology results.
8. **Writing to Files:** Outputs a FASTA file containing both DNA and amino acid sequences, along with a text file listing unique homologous species.

## Workflow Overview

1. **Input Gene Name:** User provides the name of the gene of interest as a command-line argument when executing the script.
2. **Fetch Gene Information (MyGene):** Utilizes the provided gene name to query the MyGene database and retrieve relevant gene information, such as the Entrez ID.
3. **Convert to Ensembl ID:** Converts the retrieved Entrez ID to the corresponding Ensembl ID using the MyGene database.
4. **Retrieve DNA Sequence (Ensembl):** Fetches the DNA sequence associated with the obtained Ensembl ID from the Ensembl database.
5. **Find Longest ORF (Open Reading Frame):** Identifies the longest Open Reading Frame (ORF) within the DNA sequence obtained from Ensembl.
6. **Translate DNA to Amino Acids:** Translates the longest ORF into an amino acid sequence using Biopython's translation functionality.
7. **Find Homologs (Ensembl):** Searches for homologous genes using the Ensembl ID obtained earlier.
8. **Sort Homology Results:** Sorts and extracts unique species from the homology results, filtering out the Homo sapiens species.
9. **Output Results:** Generates a FASTA file containing both DNA and amino acid sequences of the gene. Creates a text file listing unique homologous species found.

## Prerequisites
Before using Gene2FASTA, ensure you have the following dependencies installed:

- Python 3.x
- Biopython library

You can install Biopython using pip:
```
pip install biopython
```


## How to Run
```
To use Gene2FASTA, run the script with the desired gene name as the positional argument:
python3 main.py "GENE_NAME"
```


**Optional Arguments:**

## Logging (Optional):
Optionally enables verbose logging to provide detailed information about the execution process.

To use logging use the following:

```
python3 main.py "GENE_NAME" --logging
```
