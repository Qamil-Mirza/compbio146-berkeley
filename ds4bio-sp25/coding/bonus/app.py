import streamlit as st
import requests
from Bio.Data import CodonTable
from Bio import Entrez
from Bio import SeqIO
from io import StringIO
from urllib.parse import quote  # URL encoding for species search

# Set the email for api access to NCBI Database
Entrez.email = "qamilmirza@berkeley.edu"


# Function to fetch gene names for a species from NCBI
def get_gene_names(species_name):
    """
    Input: species_name (str) - the name of the species to search for
    Output: gene_names (list) - a list of gene names found for the species
    """
    try:
        search_query = f"{species_name}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_query, retmax=20)  # Fetch max 20 genes
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            return []
        
        gene_ids = record["IdList"]
        gene_names = []
        
        for gene_id in gene_ids:
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            
            if "DocumentSummarySet" in summary_record:
                doc_summary = summary_record["DocumentSummarySet"]["DocumentSummary"]
                if doc_summary:
                    gene_names.append(doc_summary[0]["Name"])
        
        return sorted(set(gene_names))  # Remove duplicates and sort
    except Exception as e:
        st.error(f"Error fetching gene names: {e}")
        return []

# Function to get the genetic code table ID for a species
def get_species_codon_table(species_name):
    """
    Input: species_name (str) - the name of the species to search for
    Output: table_id (int) - the genetic code table ID for the species
    """
    encoded_species = quote(species_name)
    search_url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{encoded_species}"
    
    try:
        response = requests.get(search_url)
        response.raise_for_status()
        data = response.json()
        
        if isinstance(data, list) and len(data) > 0:
            taxon_info = data[0]
            return int(taxon_info.get("mitochondrialGeneticCode", 1))  # Ensure integer
        else:
            st.error("No taxon information found for this species.")
            return None
    except requests.RequestException as e:
        st.error(f"Error fetching codon table information: {e}")
        return None

# Function to fetch a DNA sequence from NCBI given a species and gene name
def fetch_gene_sequence(species_name, gene_name):
    """
    Input: 
    species_name (str) - the name of the species to search for
    gene_name (str) - the name of the gene to search for
            
    Output:
    seq_record (str) - the DNA sequence of the gene
    """
    try:
        search_query = f"{species_name}[Organism] AND {gene_name}[Gene]"
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            st.error("No gene IDs found for the provided query.")
            return None
        
        gene_id = record["IdList"][0]
        fetch_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()
        
        fasta_io = StringIO(fasta_data)
        seq_record = SeqIO.read(fasta_io, "fasta")
        return str(seq_record.seq)
    
    except Exception as e:
        st.error(f"Error fetching gene sequence: {e}")
        return None

# Function to translate a DNA sequence using a given codon table
def translate_sequence(seq, table_id=1):
    """
    Input:
    seq (str) - the DNA sequence to translate
    table_id (int) - the genetic code table ID to use for translation

    Output:
    protein (str) - the translated protein sequence

    Note: If the table_id is not found, the Universal Code (1) is used by default.
    """
    if table_id not in CodonTable.unambiguous_dna_by_id:
        st.warning(f"Genetic Code ID {table_id} not found, defaulting to Universal Code (1).")
        table_id = 1

    codon_table = CodonTable.unambiguous_dna_by_id[table_id]
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) != 3:
            continue
        if codon in codon_table.stop_codons:
            
            protein += "*"
        else:
            protein += codon_table.forward_table.get(codon, "X")
    return protein


# UI RELATED CODE BELOW #
st.title("Species-Specific Genetic Code Translator 🧬")
st.write("Retrieve a DNA sequence for a species and translate it using the correct codon table.")

# Input fields for species and gene names
species_name = st.text_input("Enter a species name:", "Saccharomyces cerevisiae")

# persist gene_names in the sessions cacheeeee 0_0
if "gene_names" not in st.session_state:
    st.session_state.gene_names = []

# Button to fetch gene names, saving them to session state
if st.button("Fetch Available Genes"):
    st.session_state.gene_names = get_gene_names(species_name)
    if not st.session_state.gene_names:
        st.error("No genes found for this species. Try another one.")

# If gene names have been fetched, show a selectbox; otherwise allow manual entry
if st.session_state.gene_names:
    gene_name = st.selectbox("Select a Gene:", st.session_state.gene_names)
else:
    gene_name = st.text_input("Enter a gene name manually:")

# Button to fetch DNA sequence
if st.button("Fetch DNA Sequence"):
    if not species_name or not gene_name:
        st.error("Please enter/select both a species name and a gene name.")
    else:
        dna_sequence = fetch_gene_sequence(species_name, gene_name)
        if dna_sequence:
            st.success(f"Successfully retrieved {gene_name} gene for {species_name}!")
            st.text_area("Fetched DNA Sequence:", dna_sequence, height=150)
            
            # Fetch the correct codon table from the ENA API
            table_id = get_species_codon_table(species_name)
            if table_id is None:
                st.error("Could not find the genetic code for this species.")
            else:
                st.write(f"Genetic Code Table ID: {table_id}")
                protein = translate_sequence(dna_sequence.upper(), table_id)
                st.subheader("Translated Protein Sequence:")
                st.code(protein, language="text")
        else:
            st.error("Could not retrieve DNA sequence. Try another gene name.")
