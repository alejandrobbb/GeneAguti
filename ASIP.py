from Bio import Entrez, SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import os

#Su correo tiene que estar asociado a NCBI para poder usar la API
Entrez.email = "correo@gmai.com"  

def fetch_asip_sequence(species_name):
    # La query que uso es especifica para ese gen
    # dependiendo el gen con el vayan a trabajar pueden buscar informacion en el mismo NCBI
    # hay muchas formas de armar estas consultas y se pueden obtener datos muy robustos
    query = f"ASIP[Gene] AND {species_name}[Organism] AND mRNA"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        return None
    
    seq_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()
    
    filename = f"{species_name.replace(' ', '_')}_ASIP.fasta"
    with open(filename, "w") as output_handle:
        SeqIO.write(seq_record, output_handle, "fasta")
    
    return seq_record

species_list = [
    "Canis lupus familiaris",  # Perro
    "Mus musculus",            # Raton
    "Equus caballus"           # Caballo
]

sequences = {}
for species in species_list:
    record = fetch_asip_sequence(species)
    if record:
        sequences[species] = record


def compare_sequences(seq1, seq2, name1, name2, max_length=3000):
    # Recortamos las secuencias largas
    if len(seq1) > max_length:
        seq1 = seq1[:max_length]
    if len(seq2) > max_length:
        seq2 = seq2[:max_length]

    # yo uso alineamiento local para limitar el uso de memoria  
    alignments = pairwise2.align.localxx(seq1.seq, seq2.seq, one_alignment_only=True)

    if not alignments:
        print(f"No se pudo alinear las secuencias de {name1} y {name2}.")
        return

    alignment = alignments[0]
    identity = alignment[2] / alignment[4] * 100
    print(f"\nComparación: {name1} vs. {name2}")
    print(f"Identidad: {identity:.2f}%")
    print(format_alignment(*alignment))

# Output
species_keys = list(sequences.keys())
for i in range(len(species_keys)):
    for j in range(i + 1, len(species_keys)):
        compare_sequences(
            sequences[species_keys[i]],
            sequences[species_keys[j]],
            species_keys[i],
            species_keys[j]
        )

### los espacios y guiones son principalmente delaciones, secuencias que no estan o son diferentes entre las especies





