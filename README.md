# prepare_gg_parts.py

Standalone Python pipeline that fetches gene sequences from NCBI, codon-optimizes them for *S. cerevisiae*, removes restriction sites, designs Golden Gate assembly primers, and writes all outputs to disk.

## Dependencies

- [Biopython](https://biopython.org/) (`Bio.Entrez`, `Bio.SeqIO`, `Bio.Seq`, `Bio.SeqRecord`)
- Python standard library: `argparse`, `csv`, `datetime`, `os`, `re`, `sys`, `time`

## Usage

```bash
python prepare_gg_parts.py genes.tsv --email you@lab.edu --output ./output_sequences
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `input` | Yes | — | Input TSV file (see format below) |
| `--email` | Recommended | `user@example.com` | NCBI Entrez email address |
| `--output` / `-o` | No | `./output_sequences` | Output directory for FASTA files and reports |
| `--tsv` | No | `<output>/gg_parts_summary.tsv` | Custom path for summary TSV |
| `--quiet` / `-q` | No | False | Suppress per-step logging |

## Input TSV Format

Tab-separated, header row required.

| Column | Required | Description |
|--------|----------|-------------|
| `gene_name` | Yes | Label used for output files (e.g. `XYL1_K270R`) |
| `accession` | Yes | NCBI nucleotide or protein accession; semicolons separate fallback accessions tried in order |
| `mutation` | No | Point mutations e.g. `K270R` or `K270R;A123V` |
| `notes` | No | Free text, copied into FASTA headers |

Example:

```
gene_name	accession	mutation	notes
XYL1_K270R	X59465	K270R	S. stipitis xylose reductase
XYL2	XM_001385144		S. stipitis xylitol dehydrogenase
```

## Pipeline Steps

For each gene in the input TSV, the pipeline runs the following steps:

1. **Fetch & validate** — fetches the CDS from NCBI (nucleotide or protein accession); validates that the record matches the gene name via eSearch and record qualifiers; back-translates protein accessions to nucleotide CDS when needed
2. **Apply mutations** — applies point mutations (e.g. `K270R`) to the CDS at the amino acid level
3. **Codon optimize** — replaces each codon with the most frequent *S. cerevisiae* synonym using the Kazusa codon usage table (taxid 4932); codons below 10% relative frequency are excluded
4. **Remove restriction sites** — scans for BsaI and BsmBI recognition sites and removes any found via synonymous codon substitution
5. **Verify translation** — confirms the optimized nucleotide sequence still encodes the correct protein (exact amino acid match)
6. **Design Golden Gate primers** — adds BsaI tails to generate specific 4-nt overhangs after digestion:
   - 5' overhang: `TATG` (inclusive of ATG start codon)
   - 3' overhang: `ATCC`
   - Binding region: 18–30 nt, extended until Tm ≥ 60°C (Wallace rule)
7. **Write FASTA files** — two FASTA files per gene: original CDS and codon-optimized CDS

### Primer structure

```
5'─CGGTCTCA·T·ATG──[GOI]──ATCC·TGAGACCG─5'
        BsaI ┘ └ spacer
```

## Outputs

All files are written to the `--output` directory.

| File | Description |
|------|-------------|
| `<gene_name>_original.fasta` | Original CDS fetched from NCBI |
| `<gene_name>_optimized.fasta` | Codon-optimized, restriction-site-free CDS |
| `gg_parts_summary.tsv` | One row per gene: `gene_name`, `species`, `accession`, `fwd_primer`, `rev_primer`, `codon_optimized_sequence`, `pcr_product` |
| `pipeline_report.txt` | Human-readable run report with per-gene details: codons changed, sites removed, primer Tms, validation status, mutations applied |
