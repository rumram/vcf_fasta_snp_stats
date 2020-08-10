# vcf_fasta_snp_stats
Get synonymous, non-synonymous SNPs and track nucleotide and amino acid sequences modifications given VCF and corresponding FASTA files.

## Usage
In order to run script provide path to VCF and FASTA file. Also provide output files prefix as a third argument.
```
vcf_fasta_snp_effect.py path/to/vcf/file.vcf path/to/fasta/file.fa output_prefix
```

The above command will result with two files:
- FASTA file containig only sequences were at least one SNP were detected
- summary table in csv file format


```
SequenceName	LEN	REF	REF_COUNT	ALT	ALT_COUNT	REF_CODON	ALT_CODON	POS	POS_IN_CODON	REF_AA	ALT_AA	SNP_TYPE
Sequence1	957	A	13	T	10	CCA	CCT	276	3	P	P	S
Sequence1	957	T	6	C	6	ATA	ACA	311	2	I	T	N
Sequence2	402	C	89	G	75	CAA	GAA	331	1	Q	E	N
Sequence3	771	T	28	C	22	GAT	GAC	747	3	D	D	S
Sequence4	588	A	10	G	6	ATT	GTT	262	1	I	V	N
Sequence5	1806	C	227	T	60	CGC	TGC	529	1	R	C	N
Sequence6	321	T	108	C	87	TCT	TCC	75	3	S	S	S
```