import streamlit as st
import requests
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from io import StringIO
from urllib.parse import quote  # URL encoding for species search

# Set the email for API access to the NCBI Database
Entrez.email = "qamilmirza@berkeley.edu"

# ----------------------------- #
#        Helper Functions       #
# ----------------------------- #

def get_gene_names(species_name):
    """
    Fetch gene names for a species from NCBI.
    
    Args:
        species_name (str): e.g., "Homo sapiens"
        
    Returns:
        list: Sorted list of gene names.
    """
    try:
        search_query = f"{species_name}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_query, retmax=30)
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
        return sorted(set(gene_names))
    except Exception as e:
        st.error(f"Error fetching gene names: {e}")
        return []


def get_species_codon_table(species_name):
    """
    Retrieve the species-specific mitochondrial codon table from the ENA API.
    
    Args:
        species_name (str): e.g., "Homo sapiens"
        
    Returns:
        int: Codon table ID (default is 1 if not found)
    """
    encoded_species = quote(species_name)
    search_url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{encoded_species}"
    try:
        response = requests.get(search_url)
        response.raise_for_status()
        data = response.json()
        if isinstance(data, list) and len(data) > 0:
            taxon_info = data[0]
            return int(taxon_info.get("mitochondrialGeneticCode", 1))
        else:
            st.error("No taxon information found for this species.")
            return None
    except requests.RequestException as e:
        st.error(f"Error fetching codon table information: {e}")
        return None


def get_mrna_variants(species_name, gene_name):
    """
    Fetch all available mRNA variants from NCBI for a given species and gene.
    
    Args:
        species_name (str): e.g., "Homo sapiens"
        gene_name (str): e.g., "INS"
        
    Returns:
        dict: Dictionary mapping FASTA header (variant description) to processed mRNA sequence.
    """
    try:
        # Use a generic mRNA filter
        search_query = f"{species_name}[Organism] AND {gene_name}[Gene] AND mRNA[Filter]"
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=20)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            return {}
        
        gene_ids = record["IdList"]
        variants = {}  # key: header (variant description), value: sequence (with T's replaced by U's)
        for gene_id in gene_ids:
            fetch_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()
            
            fasta_lines = fasta_data.strip().split('\n')
            # Extract the header (without the '>' character)
            header = fasta_lines[0].strip()[1:]
            # Concatenate sequence lines and convert T's to U's
            sequence_only = ''.join(fasta_lines[1:]).replace(" ", "").replace("T", "U").replace("t", "u")
            variants[header] = sequence_only
        
        return variants
    except Exception as e:
        st.error(f"Error fetching mRNA variants: {e}")
        return {}


def translate_sequence(mrna_sequence, codon_table_id=1):
    """
    Translate an mRNA sequence into a protein sequence using the specified codon table.
    
    Args:
        mrna_sequence (str): The mRNA sequence.
        codon_table_id (int): Genetic code table ID.
        
    Returns:
        str: The translated protein sequence.
    """
    try:
        mrna_seq = Seq(mrna_sequence)
        protein_seq = mrna_seq.translate(table=codon_table_id)
        return str(protein_seq)
    except Exception as e:
        st.error(f"Error during translation: {e}")
        return ""


# ----------------------------- #
#         Streamlit UI          #
# ----------------------------- #

st.title("Species-Specific Genetic Code Translator 🧬")
st.write("Fetch mRNA variants from NCBI and choose a variant to translate using the species-specific codon table.")

# Let the user choose the input method (only NCBI fetching in this branch)
input_method = st.radio("Select input method:", ("Fetch from NCBI", "Upload FASTA file"))

if input_method == "Fetch from NCBI":
    st.subheader("Fetch mRNA from NCBI")
    species_name = st.text_input("Enter a species name:", "Homo sapiens")
    
    # Persist gene names in session state
    if "gene_names" not in st.session_state:
        st.session_state.gene_names = []
    
    if st.button("Fetch Available Genes"):
        st.session_state.gene_names = get_gene_names(species_name)
        if not st.session_state.gene_names:
            st.error("No genes found for this species. Try another one.")
    
    # If gene names are available, show a selectbox; otherwise allow manual entry
    if st.session_state.gene_names:
        gene_name = st.selectbox("Select a Gene:", st.session_state.gene_names)
    else:
        gene_name = st.text_input("Enter a gene name manually:")
    
    if st.button("Fetch mRNA Variants"):
        if not species_name or not gene_name:
            st.error("Please enter/select both a species name and a gene name.")
        else:
            variants = get_mrna_variants(species_name, gene_name)
            if variants:
                st.session_state.variants = variants
                st.success(f"Retrieved {len(variants)} mRNA variant(s) for {gene_name} in {species_name}.")
            else:
                st.error("No mRNA variants retrieved. Try another gene name.")
    
    # If variants are fetched, let the user choose one for translation
    if "variants" in st.session_state and st.session_state.variants:
        variant_keys = list(st.session_state.variants.keys())
        selected_variant = st.selectbox("Select an mRNA variant to translate:", variant_keys)
        selected_sequence = st.session_state.variants[selected_variant]
        
        st.text_area("Selected mRNA Sequence:", selected_sequence, height=150)
        st.write(f"Sequence Length: {len(selected_sequence)} nucleotides")
        
        # Retrieve the species-specific codon table (if available)
        table_id = get_species_codon_table(species_name)
        if table_id is None:
            st.error("Could not determine the genetic code for this species. Defaulting to table ID 1.")
            table_id = 1
        else:
            st.write(f"Genetic Code Table ID (from ENA): {table_id}")
        # Allow user to override the codon table if needed
        # table_id = st.number_input("Set Codon Table ID (if needed):", min_value=1, max_value=30, value=table_id)
        
        if st.button("Translate Selected mRNA"):
            protein = translate_sequence(selected_sequence.upper(), table_id)
            st.subheader("Translated Protein Sequence:")
            st.code(protein, language="text")
            
elif input_method == "Upload FASTA file":
    st.subheader("Upload a FASTA File")
    uploaded_file = st.file_uploader("Choose a FASTA file", type=["fasta", "fa", "txt"])
    
    if uploaded_file is not None:
        file_content = uploaded_file.read().decode("utf-8")
        st.text_area("Uploaded FASTA File Content:", file_content, height=150)
        # Parse the FASTA content
        fasta_io = StringIO(file_content)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        
        if not records:
            st.error("No FASTA records found in the file.")
        else:
            # If multiple records exist, let the user select one
            if len(records) > 1:
                record_ids = [record.id for record in records]
                selected_id = st.selectbox("Select a FASTA record to translate", record_ids)
                record = next((rec for rec in records if rec.id == selected_id), records[0])
            else:
                record = records[0]
            
            st.write(f"Selected Record: **{record.id}**")
            st.write(f"Sequence Length: {len(record.seq)} nucleotides")
            
            # Process sequence: replace T's with U's if needed (assume DNA input)
            processed_sequence = str(record.seq).replace("T", "U").replace("t", "u")
            st.text_area("Processed mRNA Sequence:", processed_sequence, height=150)
            
            # Let user set a codon table ID manually (or use default of 1)
            table_id = st.number_input("Enter Codon Table ID:", min_value=1, max_value=30, value=1)
            if st.button("Translate Uploaded mRNA"):
                protein = translate_sequence(processed_sequence.upper(), table_id)
                st.subheader("Translated Protein Sequence:")
                st.code(protein, language="text")
