import streamlit as st
import requests
import pandas as pd
from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from io import StringIO
from urllib.parse import quote  # URL encoding for species search
import random

# Set the email for API access to the NCBI Database
Entrez.email = "qamilmirza@berkeley.edu"

# Global variables
genetic_code_names = {
    1: "Standard Genetic Code",
    2: "Vertebrate Mitochondrial Code",
    3: "Yeast Mitochondrial Code",
    4: "Mold, Protozoan, and Coelenterate Mitochondrial Code and Mycoplasma/Spiroplasma Code",
    5: "Invertebrate Mitochondrial Code",
    6: "Ciliate, Dasycladacean, and Hexamita Nuclear Code",
    9: "Echinoderm and Flatworm Mitochondrial Code",
    10: "Euplotid Nuclear Code",
    11: "Bacterial, Archaeal, and Plant Plastid Code",
    12: "Alternative Yeast Nuclear Code",
    13: "Ascidian Mitochondrial Code",
    14: "Alternative Flatworm Mitochondrial Code",
    15: "Blepharisma Nuclear Code",
    16: "Chlorophycean Mitochondrial Code",
    21: "Trematode Mitochondrial Code",
    22: "Scenedesmus obliquus Mitochondrial Code",
    23: "Thraustochytrium Mitochondrial Code",
    24: "Rhabdopleuridae Mitochondrial Code",
    25: "Candidatus Division SR1 and Gracilibacteria Code",
    26: "Pachysolen tannophilus Nuclear Code",
    27: "Karyorelict Nuclear Code",
    28: "Condylostoma Nuclear Code",
    29: "Mesodinium Nuclear Code",
    30: "Peritrich Nuclear Code",
    31: "Blastocrithidia Nuclear Code",
    32: "Balanophoraceae Plastid Code",
    33: "Cephalodiscidae Mitochondrial UAA-Tyr Code",
}

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
        # Search the ncbi gene database for gene names pertaining to the species but only take the first 30 for brevity
        search_query = f"{species_name}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_query, retmax=30)
        record = Entrez.read(handle)
        handle.close()

        # If we cannot find anything, then just return empty list
        if not record["IdList"]:
            return []

        # otgerwise get the gene ids
        gene_ids = record["IdList"]

        # Then collect all the gene names using the gene ids to query the database
        gene_names = []
        for gene_id in gene_ids:
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            if "DocumentSummarySet" in summary_record:
                doc_summary = summary_record["DocumentSummarySet"]["DocumentSummary"]
                if doc_summary:
                    gene_names.append(doc_summary[0]["Name"])

        # I used a set to prevent duplicates and then sorted the list
        return sorted(set(gene_names))
    except Exception as e:
        st.error(f"Error fetching gene names: {e}")
        return []


def get_species_codon_table_id(species_name):
    """
    Retrieve the species-specific mitochondrial codon table ID from the ENA API.

    Args:
        species_name (str): e.g., "Homo sapiens"

    Returns:
        int: Codon table ID (default is 1 if not found)
    """
    encoded_species = quote(species_name)
    search_url = (
        f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{encoded_species}"
    )
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
        search_query = (
            f"{species_name}[Organism] AND {gene_name}[Gene] AND mRNA[Filter]"
        )
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=20)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return {}

        gene_ids = record["IdList"]
        variants = (
            {}
        )  # key: header (variant description), value: sequence (with T's replaced by U's)
        for gene_id in gene_ids:
            fetch_handle = Entrez.efetch(
                db="nucleotide", id=gene_id, rettype="fasta", retmode="text"
            )
            fasta_data = fetch_handle.read()
            fetch_handle.close()

            # Here I'm just formatting the Fasta data to
            fasta_lines = fasta_data.strip().split("\n")
            header = fasta_lines[0].strip()[1:]
            sequence_only = (
                "".join(fasta_lines[1:])
                .replace(" ", "")
                .replace("T", "U")
                .replace("t", "u")
            )
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
        # Here i convert to a sequence object so I can use the translate method
        mrna_seq = Seq(mrna_sequence)
        protein_seq = mrna_seq.translate(table=codon_table_id)
        return str(protein_seq)
    except Exception as e:
        st.error(f"Error during translation: {e}")
        return ""


def reverse_translation_count(protein, table_id):
    """
    Given a protein sequence and a species-specific codon table ID,
    return the total number of distinct DNA sequences that could encode the protein.

    This function uses the species-specific codon dictionary (including stop codons)
    to compute the product of the number of codons for each amino acid.

    Args:
        protein (str): Translated protein sequence.
        table_id (int): The species-specific codon table ID.

    Returns:
        int: Total number of possible DNA sequences.
    """
    # Get the species-specific codon table
    species_table = CodonTable.unambiguous_dna_by_id[table_id]
    codon_dict = species_table.forward_table

    # Build a mapping from amino acid -> count of codons that encode it
    aa_to_count = {}
    for codon, aa in codon_dict.items():
        aa_to_count[aa] = aa_to_count.get(aa, 0) + 1

    # For stop codons, use the stop codons list
    stop_count = len(species_table.stop_codons)

    # Solve the combinatorial problem
    total_options = 1
    for aa in protein:
        if aa == "*":
            total_options *= stop_count
        else:
            total_options *= aa_to_count.get(aa, 0)
    return total_options


def reverse_translate_random(protein, target_table_id):
    """
    Given a protein sequence and a target species codon table ID,
    return an mRNA sequence where each amino acid is encoded by
    a randomly chosen codon from the target species codon table.

    Args:
        protein (str): Translated protein sequence.
        target_table_id (int): Codon table ID for the target species.

    Returns:
        str: mRNA sequence (with U's) reverse-translated using target species codons.
    """
    target_table = CodonTable.unambiguous_dna_by_id[target_table_id]
    aa_to_codons = {}
    for codon, aa in target_table.forward_table.items():
        aa_to_codons.setdefault(aa, []).append(codon)
    aa_to_codons["*"] = target_table.stop_codons

    mrna_seq = ""
    for aa in protein:
        if aa not in aa_to_codons:
            raise ValueError(f"Invalid amino acid in protein sequence: {aa}")
        mrna_seq += random.choice(aa_to_codons[aa])
    # Convert DNA codons (with T's) to mRNA codons (with U's)
    return mrna_seq.replace("T", "U")


def display_codon_table(table_id):
    """
    Display the species-specific codon table as a scrollable DataFrame.

    Args:
        table_id (int): The species-specific codon table ID.
    """
    species_table = CodonTable.unambiguous_dna_by_id[table_id]
    codon_mapping = species_table.forward_table
    stop_codons = species_table.stop_codons

    table_data = []
    for codon, aa in codon_mapping.items():
        table_data.append({"Codon": codon, "Amino Acid": aa})
    for codon in stop_codons:
        table_data.append({"Codon": codon, "Amino Acid": "*"})

    df = pd.DataFrame(table_data)
    df = df.sort_values(by="Codon")

    st.subheader(f"{genetic_code_names[table_id]} Codon Table")
    st.dataframe(df, height=300)


# ----------------------------- #
#         Streamlit UI          #
# ----------------------------- #

st.title("Species-Specific Genetic Code Translator 🧬")
st.caption("Developed by: Janis Prak, Alyssa Chew, Skye Leng, Abby Huang, Qamil Mirza")
st.write(
    "The species-specific genetic code translator allows you to fetch mRNA sequences from NCBI to translate to protein sequences. The species-specific codon table is retrieved from NCBI based on the mitochondrial genetic code ID for the species as defined on the ENA API."
)

st.divider()

st.subheader("Fetch mRNA from NCBI")
st.markdown("Enter a species name and a gene name to fetch mRNA variants. Or, fetch available gene names from NCBI but note homo sapiens gene name fetching is not yet supported in this version!")
st.markdown("""
**Some examples** in the form `(species, gene name)`:

1. (Homo sapiens, INS)  
2. (Mus musculus, Ins2)
3. (Saccharomyces cerevisiae, CYC1)
""")
species_name = st.text_input("Enter a species name:", "Saccharomyces cerevisiae")

# Let the user choose how to provide the gene name.
gene_input_mode = st.radio(
    "Choose Gene Input Method", ("Enter Manually", "Fetch Available Genes")
)

if gene_input_mode == "Fetch Available Genes":
    # If the user chooses to fetch genes, let them click a button to fetch them.
    if st.button("Fetch Available Genes"):
        st.session_state.gene_names = get_gene_names(species_name)
        if not st.session_state.gene_names:
            st.error("No genes found for this species. Try another one.")
    # If gene names are available, display a selectbox.
    if "gene_names" in st.session_state and st.session_state.gene_names:
        gene_name = st.selectbox("Select a Gene:", st.session_state.gene_names)
    else:
        gene_name = ""  # No gene selected yet.
elif gene_input_mode == "Enter Manually":
    # When manual entry is selected, clear any previously fetched gene names.
    st.session_state.gene_names = []
    gene_name = st.text_input("Enter a gene name manually:")

if st.button("Fetch mRNA Variants"):
    if not species_name or not gene_name:
        st.error("Please enter/select both a species name and a gene name.")
    else:
        variants = get_mrna_variants(species_name, gene_name)
        if variants:
            st.session_state.variants = variants
            st.success(
                f"Retrieved {len(variants)} mRNA variant(s) for {gene_name} in {species_name}."
            )
        else:
            st.error("No mRNA variants retrieved. Try another gene name.")


# If variants are fetched, let the user choose one for translation
if "variants" in st.session_state and st.session_state.variants:
    variant_keys = list(st.session_state.variants.keys())
    selected_variant = st.selectbox(
        "Select an mRNA variant to translate:", variant_keys
    )
    selected_sequence = st.session_state.variants[selected_variant]

    st.text_area("Selected mRNA Sequence:", selected_sequence, height=150)
    st.write(f"Sequence Length: {len(selected_sequence)} nucleotides")

    # Retrieve the species-specific codon table ID
    table_id = get_species_codon_table_id(species_name)
    if table_id is None:
        st.error(
            "Could not determine the genetic code for this species. Defaulting to table ID 1."
        )
        table_id = 1

    code_name = genetic_code_names.get(table_id, "Unknown")
    st.write(f"Genetic Code Table ID: {table_id} - {code_name}")

    # Display the codon table so the user can see the mapping.
    display_codon_table(table_id)

    if st.button("Translate Selected mRNA"):
        protein = translate_sequence(selected_sequence.upper(), table_id)
        st.session_state["translated_protein"] = protein


if "translated_protein" in st.session_state:
    st.subheader("Translated Protein Sequence:")
    st.code(st.session_state["translated_protein"], language="text")
    st.markdown("🎉 Voila! You have the translated protein sequence from the selected mRNA sequence.")
    st.divider()
    st.subheader("Reverse Translation Analysis (Bonus Feature)")
    st.write("While our main focus was only on translating mRNA to protein for any species, we decided to include reverse translation for fun since we had the necessary information to do so! Consider this a bonus feature and have fun with it!")
    rev_count = reverse_translation_count(
        st.session_state["translated_protein"], table_id
    )

    st.subheader("Number Of Distinct DNA Sequences")
    st.markdown("Here, we calculate the total number of distinct DNA sequences that could encode the protein sequence. This calculation is done by multiplying the number of codons that encode each amino acid.")
    st.text_area(
        "Total Number of Distinct DNA Sequences:",
        str(rev_count),
        height=150,
    )
    st.subheader("Reverse Translation With Another Species Codon Table")
    st.write("In this section, you can reverse translate the protein sequence to an mRNA sequence using the codon table of another species.")
    target_species = st.text_input("Enter a target species name:", "Mus musculus")
    if st.button("Reverse Translate Protein", key="reverse_translate_fetch"):
        target_species_table_id = get_species_codon_table_id(target_species)
        if target_species_table_id is None:
            st.error(
                "Could not determine the genetic code for the target species. Defaulting to table ID 1."
            )
            target_species_table_id = 1

        target_code_name = genetic_code_names.get(target_species_table_id, "Unknown")
        st.write(
            f"Genetic Code Table ID for {target_species}: {target_species_table_id} - {target_code_name}"
        )
        display_codon_table(target_species_table_id)
        reverse_mrna = reverse_translate_random(
            st.session_state["translated_protein"], target_species_table_id
        )
        st.text_area(
            f"Reverse-Translated mRNA Sequence to {target_species}:",
            reverse_mrna,
            height=150,
        )

# ----------------------------- #
#         Acknowledgements      #
# ----------------------------- #
st.divider()
st.subheader("Data Sources & Acknowledgements")
st.markdown("""  
- **NCBI Gene and Nucleotide Databases**:  
  This app uses data retrieved via the NCBI Entrez API. See [NCBI's website](https://www.ncbi.nlm.nih.gov/) 
  and [Entrez Programming Utilities (E-utilities)](https://www.ncbi.nlm.nih.gov/books/NBK25501/) 
  for more information and usage policies.

- **EBI ENA API**:  
  Mitochondrial genetic code information is fetched from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) 
  via their Taxonomy REST API. Please refer to the [ENA documentation](https://www.ebi.ac.uk/ena/portal/api/) 
  for further details and terms of use.
""")
