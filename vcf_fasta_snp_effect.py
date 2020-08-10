# Create summary tables from VCF file and corresponding FASTA file

import sys
import csv
import pandas as pd

vcf_file = open(sys.argv[1], 'r')
fasta_file = open(sys.argv[2], 'r')
output_prefix = sys.argv[3]


def create_list_of_dicts(file):
    clean_subset = []

    for line in file:  # Remove rows starting with '#'
        if not line.startswith('#'):
            clean_subset.append(line)

    colnames = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "HEADERS", "FORMAT"
        ]
    clean_subset_with_header = []

    for line in csv.DictReader(clean_subset, fieldnames=colnames, delimiter='\t'):
        clean_subset_with_header.append(dict(line))

    return clean_subset_with_header


list_of_dicts = create_list_of_dicts(vcf_file)


def create_seq_dict(file):
    # Create dictionary of header:sequence pairs
    fasta = {}
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            sequence_header = line[1:]
            if sequence_header not in fasta:
                fasta[sequence_header] = []
            continue
        sequence = line
        fasta[sequence_header].append(sequence)

    headers = []
    sequences = []

    for key, value in fasta.items():  # Handle multiline fasta
        conc_seq = "".join(value)
        headers.append(str(key))
        sequences.append(str(conc_seq))

    dictionary = dict(zip(headers, sequences))
    return dictionary


dict_of_seqs = create_seq_dict(fasta_file)


def start_nuc(position):
    # Determine first nucleotide of codon
    start_position = []
    if position % 3 == 0:
        start_position = position - 2
    elif position % 3 == 2:
        start_position = position - 1
    elif position % 3 == 1:
        start_position = position

    return start_position


def nuc_pos_in_codon(raw_position):
    # Determine nucleotide position in codon
    pos_in_codon = []
    if raw_position % 3 == 0:
        pos_in_codon = 2
    elif raw_position % 3 == 2:
        pos_in_codon = 1
    elif raw_position % 3 == 1:
        pos_in_codon = 0

    return pos_in_codon


def translation(nucleotide_seq, start):
    # Translate nucleotide codon to amino acid
    threes = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

    codon = nucleotide_seq[start - 1:start + 2]
    amino_acid = threes.get(codon)
    return amino_acid


def get_alt_codon(codon, position, alt_position):
    # Get alternative codon
    alt_codon = codon[:position] + alt_position + codon[position + 1:]
    return alt_codon


def get_snp_descriptor(list_of_dicts_obj):
    # Create table describing each SNP
    snp_table = []
    snp_fasta = ''

    for dictionary in list_of_dicts_obj:
        if str(dictionary.get('INFO').split(";")[-1]) == 'TYPE=snp' and (
                len(dictionary.get('ALT')) == 1):
            chrom = dictionary.get("CHROM")
            seq = str(dict_of_seqs.get(chrom))
            pos1 = int(dictionary.get("POS"))
            pos2 = start_nuc(pos1)
            triple = seq[pos2 - 1:pos2 + 2]
            new_codon = get_alt_codon(triple, nuc_pos_in_codon(pos1), dictionary.get('ALT'))
            ref_aa = str(translation(triple, 1))
            alt_aa = str(translation(new_codon, 1))
            if ref_aa == alt_aa:
                snp_type = 'S'
            else:
                snp_type = 'N'

            snp_table.append([
                chrom,
                len(seq),
                dictionary.get('REF'),
                dictionary.get("FORMAT").split(":")[-5],
                dictionary.get('ALT'),
                dictionary.get("FORMAT").split(":")[-3],
                triple,
                new_codon,
                str(dictionary.get('POS')),
                str(nuc_pos_in_codon(pos1) + 1),
                ref_aa,
                alt_aa,
                snp_type
                ])

            # Create fasta file of each sequence containing SNP
            snp_fasta += (">" + str(chrom) + "\n" + seq + "\n")

    f = open(str(output_prefix) + ".fa", 'w')
    f.write(str(snp_fasta))
    f.close()

    col_names = [
        'SEQ_NAME', 'LEN', 'REF', 'REF_COUNT', 'ALT', 'ALT_COUNT',
        'REF_CODON', 'ALT_CODON', 'POS', 'POS_IN_CODON', 'REF_AA', 'ALT_AA',
        'SNP_TYPE'
        ]

    snp_table_df = pd.DataFrame(snp_table, columns=col_names)
    snp_table_df.to_csv(str(output_prefix) + ".csv", sep='\t', index=False)


get_snp_descriptor(list_of_dicts)
