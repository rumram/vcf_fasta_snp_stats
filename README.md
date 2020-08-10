# vcf_fasta_snp_stats
Get synonymous, non-synonymous SNPs and track nucleotide and amino acid sequences modifications given VCF and corresponding FASTA files.
The Python script provides description of each SNP detected within VCF file and predicts how the SNP affects the nucleotide and amino acid sequence.

## Prerequisities
- pandas python library

## Usage
In order to run script provide path to VCF and FASTA file. Also provide output files prefix as a third argument.
```
vcf_fasta_snp_effect.py path/to/vcf/file.vcf path/to/fasta/file.fa output_prefix
```

### Example:
```
vcf_fasta_snp_effect.py /home/rumram/file.vcf /home/rumram/file.fa test
```

The above command will result with two files:
- FASTA file containing only sequences were at least one SNP was detected (test.fa),
- summary table in csv file format (test.csv).

### Summary table header:
| SEQ_NAME | LEN |	REF |	REF_COUNT |	ALT |	ALT_COUNT |	REF_CODON |	ALT_CODON |	POS |	POS_IN_CODON |	REF_AA |	ALT_AA |	SNP_TYPE |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Sequence1 |	957 |	A |	13 |	T |	10 |	CCA |	CCT |	276 |	3 |	P |	P |	S |
| Sequence1	| 957 |	T |	6 |	C |	6 |	ATA |	ACA |	311 |	2 |	I |	T |	N |
| Sequence2 |	402 |	C |	89 |	G |	75 |	CAA |	GAA |	331 |	1 |	Q |	E |	N |
| Sequence3 |	771 |	T |	28 |	C |	22 |	GAT |	GAC |	747 |	3 |	D |	D |	S |
| Sequence4 |	588 |	A |	10 |	G |	6 |	ATT |	GTT |	262 |	1 |	I |	V |	N |
| Sequence4 |	1806 |	C |	227 |	T |	60 |	CGC |	TGC |	529 |	1 |	R |	C |	N |
| Sequence5 |	321 |	T |	108 |	C |	87 |	TCT |	TCC |	75 |	3 |	S |	S |	S |


### Headers description:
- 'SEQ_NAME': name of FASTA sequence,
- 'LEN': length of nucleotide sequence,
- 'REF': reference nucleotide,
- 'REF_COUNT': reference nucleotide counts,
- 'ALT': alternative nucleotide,
- 'ALT_COUNT': alternative nucleotide counts,
- 'REF_CODON': reference codon,
- 'ALT_CODON': alternative codon,
- 'POS': location of SNP within sequence,
- 'POS_IN_CODON': position of SNP within codon (1-3),
- 'REF_AA': reference amino acid,
- 'ALT_AA': amino acid change as a result of alternative nucleotide occurrence,
- 'SNP_TYPE': indicates if SNP is synonymous (S) or non-synonymous (N).
