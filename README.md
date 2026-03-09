# M23 Metalloendopeptidase Domain Analysis Pipeline

A Nextflow (DSL2) computational pipeline designed to identify and analyze protein sequences containing the M23 metalloendopeptidase domain (PFAM: PF01551) within bacterial genome assemblies. The pipeline searches for the M23 domain, identifies other co-occurring domains, annotates their genomic locations, and generates visual charts of domain architectures.

## Pipeline Steps

The workflow executes the following processes:

* **extractorfs**: Extracts Open Reading Frames (ORFs, $\ge$ 900 bp) from uncompressed genomic sequences and translates them into protein sequences using a custom C++ application.
* **uniquetrans**: Creates a non-redundant (unique) set of protein sequences while preserving cluster abundance information in separate TSV files (`unique_clusts.tsv` and `unique_lengths.tsv`).
* **splitfasta**: Splits the unique sequences file into a specified number of chunks for parallel processing.
* **hmmfetch**: Extracts the Peptidase M23 domain model (PF01551) from the provided PFAM database.
* **hmmselsearch**: Searches the non-redundant protein sequence chunks for the M23 domain.
* **extracttrans**: Extracts full protein sequences that contain the M23 domain from the initial sequence sets.
* **hmmallsearch**: Performs a secondary search on the M23-positive protein sequences to find all other domains available in the PFAM database.
* **integrate**: Merges the identified M23-containing protein sequences and their corresponding domain search results into comprehensive summary files using the `collectFile()` operator.
* **annotdom**: Generates a GFF3 annotation file based on the search results, indicating the exact locations of all found domains on the protein sequences.
* **archchart**: Renders an HTML schematic representing the domain architectures of all analyzed M23-containing proteins, incorporating cluster frequency data.

## Requirements

* **Nextflow** (DSL2 enabled).
* **Conda**: All Python scripts and tools are executed within a unified Conda environment defined in `envs/powb-pyhmmer.yml`.

## Usage

To run the pipeline, use the following Nextflow command with the `-with-conda` argument:

```bash
nextflow run main.nf -with-conda
