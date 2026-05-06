#!/usr/bin/env python3
"""
prepare_gg_parts.py — Fetch, codon-optimize, and design Golden Gate primers

PIPELINE
--------
1. Read a TSV gene list (gene_name, accession, [mutation], [notes]).
2. Fetch each CDS from NCBI; validate the gene name matches the record.
3. Apply optional point mutations (e.g. K270R).
4. Codon-optimize for S. cerevisiae (Kazusa table, taxid 4932).
5. Remove all BsaI / BsmBI recognition sites via synonymous substitution.
6. Design PCR primers that introduce BsaI sites:

     5'─CGGTCTCA·T·ATG──[GOI]──ATCC·TGAGACCG─5'
           BsaI ┘ └ spacer
     4-nt overhangs after BsaI digestion: TATG (5') and ATCC (3')

7. Write output TSV and per-gene FASTA files.

INPUT TSV (tab-separated, header required)
------------------------------------------
  gene_name   accession       mutation   notes
  XYL1_K270R  X59465          K270R      S. stipitis xylose reductase
  XYL2        XM_001385144               S. stipitis xylitol dehydrogenase

  gene_name   Required. Label for output files.
  accession   Required. NCBI nucleotide or protein accession; semicolon-
              separated list tried in order (first success wins).
  mutation    Optional. Point mutations e.g. K270R or K270R;A123V.
  notes       Optional. Copied to FASTA headers.

OUTPUT TSV COLUMNS
------------------
  gene_name | species | accession | fwd_primer | rev_primer |
  codon_optimized_sequence | pcr_product

USAGE
-----
  python prepare_gg_parts.py genes.tsv
  python prepare_gg_parts.py genes.tsv --output ./parts --email you@lab.edu
"""

import argparse
import csv
import datetime
import os
import re
import sys
import time
from io import StringIO

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_EMAIL = "user@example.com"
ENTREZ_DELAY  = 0.4    # seconds between NCBI requests
OUTPUT_DIR    = "./output_sequences"

# BsaI Golden Gate primer tails
#   Fwd: BsaI (CGGTCTC) + spacer A + extra T → after cut: TATG 5' overhang
#   Rev: BsaI (CGGTCTC) + spacer A + GGAT (RC of ATCC) → after cut: ATCC 5' overhang
GG_FWD_TAIL   = "CGGTCTCAT"     # 9 nt  (appended before binding region)
GG_REV_TAIL   = "CGGTCTCAGGAT"  # 12 nt (appended before RC binding region)

# RC of GG_REV_TAIL — appended to 3' end of PCR-product top strand
_GG_PCR_SUFFIX = "ATCCTGAGACCG"

# Primer binding parameters
BIND_MIN   = 18    # minimum nt to try for binding region
BIND_MAX   = 30    # maximum nt (hard cap)
TM_TARGET  = 60.0  # °C, minimum Tm for binding portion

# ─────────────────────────────────────────────────────────────────────────────
# S. cerevisiae codon usage — Kazusa, taxid 4932 (occurrences per 1 000 codons)
# ─────────────────────────────────────────────────────────────────────────────

SC_CODON_FREQ = {
    "TTT": 13.0, "TTC": 20.8,
    "TTA": 14.1, "TTG": 27.2, "CTT": 12.3, "CTC":  6.4,
    "CTA": 14.1, "CTG": 10.5,
    "ATT": 18.2, "ATC": 17.2, "ATA": 17.7,
    "ATG": 20.9,
    "GTT": 18.8, "GTC": 12.0, "GTA": 11.8, "GTG": 11.2,
    "TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG":  8.6,
    "AGT": 14.2, "AGC":  9.8,
    "CCT": 13.5, "CCC":  6.8, "CCA": 18.3, "CCG":  5.3,
    "ACT": 20.3, "ACC": 12.7, "ACA": 17.9, "ACG":  6.7,
    "GCT": 21.2, "GCC": 12.6, "GCA": 16.3, "GCG":  6.2,
    "TAT": 12.2, "TAC": 15.3,
    "TAA":  1.0, "TAG":  0.3, "TGA":  0.5,
    "CAT": 10.9, "CAC":  7.8,
    "CAA": 27.3, "CAG": 12.1,
    "AAT": 17.6, "AAC": 25.1,
    "AAA": 41.9, "AAG": 30.8,
    "GAT": 37.6, "GAC": 20.2,
    "GAA": 45.0, "GAG": 19.2,
    "TGT":  8.1, "TGC":  5.9,
    "TGG": 10.4,
    "CGT":  6.4, "CGC":  2.6, "CGA":  3.4, "CGG":  2.0,
    "AGA": 21.3, "AGG":  9.2,
    "GGT": 23.9, "GGC":  9.8, "GGA": 10.8, "GGG":  6.0,
}

AA_TO_CODONS = {
    'F': ['TTT','TTC'], 'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
    'I': ['ATT','ATC','ATA'], 'M': ['ATG'],
    'V': ['GTT','GTC','GTA','GTG'],
    'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
    'P': ['CCT','CCC','CCA','CCG'], 'T': ['ACT','ACC','ACA','ACG'],
    'A': ['GCT','GCC','GCA','GCG'], 'Y': ['TAT','TAC'],
    '*': ['TAA','TAG','TGA'],
    'H': ['CAT','CAC'], 'Q': ['CAA','CAG'], 'N': ['AAT','AAC'],
    'K': ['AAA','AAG'], 'D': ['GAT','GAC'], 'E': ['GAA','GAG'],
    'C': ['TGT','TGC'], 'W': ['TGG'],
    'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
    'G': ['GGT','GGC','GGA','GGG'],
}


def _build_best_codon_map(min_rel_freq=0.10):
    best = {}
    for aa, codons in AA_TO_CODONS.items():
        freqs = {c: SC_CODON_FREQ.get(c, 0) for c in codons}
        total = sum(freqs.values())
        eligible = {c: f for c, f in freqs.items()
                    if total == 0 or (f / total) >= min_rel_freq}
        pool = eligible if eligible else freqs
        best[aa] = max(pool, key=lambda c: pool[c])
    best['*'] = 'TAA'
    return best


BEST_CODON  = _build_best_codon_map()
CODON_TO_AA = {c: aa for aa, codons in AA_TO_CODONS.items() for c in codons}

# ─────────────────────────────────────────────────────────────────────────────
# Restriction site patterns — BsaI and BsmBI, both strands
# ─────────────────────────────────────────────────────────────────────────────

RE_SITES = {
    'BsaI_fwd':  re.compile(r'GGTCTC'),
    'BsaI_rev':  re.compile(r'GAGACC'),
    'BsmBI_fwd': re.compile(r'CGTCTC'),
    'BsmBI_rev': re.compile(r'GAGACG'),
}

# ─────────────────────────────────────────────────────────────────────────────
# Sequence utilities
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def translate(nt_seq: str) -> str:
    """Translate CDS to protein, stopping at first stop codon."""
    return str(Seq(nt_seq).translate(to_stop=True))


def back_translate(protein: str) -> str:
    """Back-translate protein using best S. cerevisiae codon per AA; add TAA stop."""
    return ''.join(BEST_CODON.get(aa, 'NNN') for aa in protein) + 'TAA'


def find_re_sites(seq: str) -> list:
    hits = []
    for name, pat in RE_SITES.items():
        for m in pat.finditer(seq):
            ctx_s = max(0, m.start() - 10)
            ctx_e = min(len(seq), m.end() + 10)
            hits.append({
                'enzyme':  name,
                'start':   m.start(),
                'end':     m.end(),
                'context': seq[ctx_s:ctx_e],
            })
    return hits


def tm_wallace(seq: str) -> float:
    """Wallace-rule Tm for the binding portion of a primer (°C)."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    n  = len(seq)
    if n == 0:
        return 0.0
    if n < 14:
        return float(2 * at + 4 * gc)
    return 64.9 + 41.0 * (gc - 16.4) / n

# ─────────────────────────────────────────────────────────────────────────────
# Point mutations
# ─────────────────────────────────────────────────────────────────────────────

def apply_mutation(cds: str, mutation: str) -> tuple:
    """
    Apply a point mutation formatted as <OrigAA><1-based position><NewAA> (e.g. K270R).
    Returns (mutated_cds, annotation_string).
    """
    m = re.fullmatch(r'([A-Z])(\d+)([A-Z])', mutation.strip())
    if not m:
        raise ValueError(f"Unrecognised mutation '{mutation}' — expected format e.g. K270R")
    orig_aa, pos_str, new_aa = m.group(1), m.group(2), m.group(3)
    pos = int(pos_str)       # 1-based protein position
    nt_s = (pos - 1) * 3
    nt_e = nt_s + 3

    if nt_e > len(cds):
        raise ValueError(
            f"Mutation {mutation}: position {pos} exceeds CDS length ({len(cds)//3} codons)")

    current_codon = cds[nt_s:nt_e].upper()
    current_aa    = CODON_TO_AA.get(current_codon, '?')
    if current_aa != orig_aa:
        raise ValueError(
            f"Mutation {mutation}: codon {pos} is {current_codon} "
            f"({current_aa}), expected {orig_aa}")

    new_codon = BEST_CODON[new_aa]
    mutated   = cds[:nt_s] + new_codon + cds[nt_e:]
    annotation = f"{mutation} ({current_codon}→{new_codon} at nt {nt_s+1}-{nt_e})"
    return mutated, annotation

# ─────────────────────────────────────────────────────────────────────────────
# Codon optimization
# ─────────────────────────────────────────────────────────────────────────────

def codon_optimize(cds: str) -> tuple:
    """
    Replace each codon with the highest-frequency S. cerevisiae synonym.
    Preserves ATG start; replaces stop codon with TAA.
    Returns (optimized_nt, n_codons_changed, opt_score_pct).
    """
    if len(cds) % 3 != 0:
        raise ValueError(f"CDS length {len(cds)} is not divisible by 3")

    aa_seq   = str(Seq(cds).translate())   # includes stop codon character
    codons   = [cds[i:i+3].upper() for i in range(0, len(cds), 3)]
    out      = []
    changed  = 0
    best_cnt = 0

    for i, (codon, aa) in enumerate(zip(codons, aa_seq)):
        best = 'ATG' if i == 0 else BEST_CODON.get(aa, codon)
        if codon == best:
            best_cnt += 1
        else:
            changed += 1
        out.append(best)

    opt_nt = ''.join(out)
    score  = round(100 * best_cnt / len(codons), 1)
    return opt_nt, changed, score

# ─────────────────────────────────────────────────────────────────────────────
# Restriction site removal
# ─────────────────────────────────────────────────────────────────────────────

def remove_re_sites(cds: str, gene_name: str) -> tuple:
    """
    Iteratively remove all BsaI/BsmBI sites via synonymous single-codon substitutions.
    Returns (edited_cds, list_of_edits).
    """
    seq   = cds.upper()
    edits = []

    for _ in range(100):
        hits = find_re_sites(seq)
        if not hits:
            break

        hit   = hits[0]
        fixed = False

        for offset in range(hit['end'] - hit['start']):
            nt_pos    = hit['start'] + offset
            codon_idx = nt_pos // 3
            cs, ce    = codon_idx * 3, codon_idx * 3 + 3

            if ce > len(seq):
                continue

            current = seq[cs:ce]
            aa      = CODON_TO_AA.get(current)
            if aa is None or aa == '*':
                continue
            if codon_idx == 0 and current == 'ATG':
                continue    # never touch the start codon

            # Rank alternatives by frequency (descending); skip current
            alts = sorted(
                [(c, SC_CODON_FREQ.get(c, 0))
                 for c in AA_TO_CODONS.get(aa, []) if c != current],
                key=lambda x: x[1], reverse=True,
            )

            win_s = max(0, hit['start'] - 20)
            win_e = min(len(seq), hit['end'] + 20)

            for alt_codon, freq in alts:
                candidate = seq[:cs] + alt_codon + seq[ce:]

                # Site must be destroyed and no new sites in ±20 nt window
                site_gone = not any(
                    p.search(candidate[max(0, hit['start']-2):hit['end']+2])
                    for p in RE_SITES.values()
                )
                if (site_gone and
                        len(find_re_sites(candidate[win_s:win_e])) <=
                        len(find_re_sites(seq[win_s:win_e]))):
                    seq = candidate
                    edits.append({
                        'enzyme':   hit['enzyme'],
                        'site_pos': f"{hit['start']}–{hit['end']}",
                        'codon':    codon_idx + 1,
                        'change':   f"{current}→{alt_codon} ({aa}, {freq:.1f}/1000)",
                    })
                    fixed = True
                    break

            if fixed:
                break

        if not fixed:
            raise RuntimeError(
                f"{gene_name}: cannot remove {hit['enzyme']} at "
                f"{hit['start']}–{hit['end']}. Context: {hit['context']}")

    remaining = find_re_sites(seq)
    if remaining:
        raise RuntimeError(
            f"{gene_name}: {len(remaining)} site(s) remain after editing: {remaining}")
    return seq, edits

# ─────────────────────────────────────────────────────────────────────────────
# Golden Gate primer design
# ─────────────────────────────────────────────────────────────────────────────

def design_primers(cds_nt: str) -> tuple:
    """
    Design BsaI Golden Gate primers for a codon-optimized CDS.

    Forward primer: 5'-CGGTCTCAT-[binding]-3'
      Binding starts at ATG (nt 0 of CDS).
      After BsaI digestion → TATG 5' overhang (T from tail + ATG start codon).

    Reverse primer: 5'-CGGTCTCAGGAT-[binding_RC]-3'
      Binding is RC of the 3' end of CDS (includes stop codon).
      After BsaI digestion → ATCC 5' overhang (RC of GGAT).

    PCR product (top strand): CGGTCTCAT·[full CDS]·ATCCTGAGACCG

    Returns (fwd_primer, rev_primer, pcr_product).
    """
    if not cds_nt.upper().startswith('ATG'):
        raise ValueError("CDS does not start with ATG — cannot design forward primer")

    def _pick_binding(seq_slice, from_5prime=True):
        """Extend from BIND_MIN to BIND_MAX until Tm ≥ TM_TARGET."""
        for n in range(BIND_MIN, BIND_MAX + 1):
            bind = seq_slice[:n] if from_5prime else seq_slice[-n:]
            if tm_wallace(bind) >= TM_TARGET:
                return bind
        return seq_slice[:BIND_MAX] if from_5prime else seq_slice[-BIND_MAX:]

    fwd_bind = _pick_binding(cds_nt, from_5prime=True)
    rev_bind = reverse_complement(_pick_binding(cds_nt, from_5prime=False))

    fwd_primer  = GG_FWD_TAIL + fwd_bind
    rev_primer  = GG_REV_TAIL + rev_bind
    pcr_product = GG_FWD_TAIL + cds_nt + _GG_PCR_SUFFIX

    return fwd_primer, rev_primer, pcr_product

# ─────────────────────────────────────────────────────────────────────────────
# NCBI fetch and validation
# ─────────────────────────────────────────────────────────────────────────────

def _ncbi_fetch(db: str, accession: str, rettype: str) -> str | None:
    for attempt in range(3):
        try:
            h    = Entrez.efetch(db=db, id=accession, rettype=rettype, retmode='text')
            data = h.read()
            h.close()
            time.sleep(ENTREZ_DELAY)
            return data
        except Exception as e:
            print(f"    [retry {attempt+1}] {db}/{accession}: {e}", file=sys.stderr)
            time.sleep(1.5)
    return None


def _is_protein_accession(accession: str) -> bool:
    acc = accession.split('.')[0]
    if acc.startswith(('NP_', 'XP_', 'YP_', 'WP_', 'AP_', 'SP_')):
        return True
    if re.fullmatch(r'[A-Z]{3}\d{5,}', acc):    # GenBank protein (e.g. ABX42514)
        return True
    if re.fullmatch(r'[A-NR-Z]\d[A-Z][A-Z0-9]{2}\d', acc):   # UniProt
        return True
    if re.fullmatch(r'[OPQ]\d[A-Z0-9]{3}\d', acc):            # UniProt OPQ
        return True
    return False


def _extract_from_gb(gb_text: str, gene_hint: str | None = None) -> dict | None:
    """
    Parse GenBank text. Return dict with cds_nt, protein, species, validated.
    validated = True when gene_hint was found in a CDS qualifier.
    """
    records = list(SeqIO.parse(StringIO(gb_text), 'genbank'))
    if not records:
        return None

    for rec in records:
        species   = rec.annotations.get('organism', 'unknown organism')
        cds_feats = [f for f in rec.features if f.type == 'CDS']
        if not cds_feats:
            continue

        # Prefer a CDS whose qualifiers mention gene_hint
        if gene_hint:
            for feat in cds_feats:
                q = feat.qualifiers
                names = (q.get('gene', []) + q.get('locus_tag', []) +
                         q.get('product', []))
                if any(gene_hint.lower() in n.lower() for n in names):
                    cds_nt  = str(feat.extract(rec.seq)).upper()
                    protein = q.get('translation', [''])[0]
                    if len(cds_nt) > 90:
                        return dict(cds_nt=cds_nt, protein=protein,
                                    species=species, validated=True)

        # Fallback: first CDS in the record
        feat    = cds_feats[0]
        cds_nt  = str(feat.extract(rec.seq)).upper()
        protein = feat.qualifiers.get('translation', [''])[0]
        if len(cds_nt) > 90:
            return dict(cds_nt=cds_nt, protein=protein,
                        species=species, validated=False)
    return None


def _esearch_validate(gene_hint: str, accession: str) -> tuple:
    """
    Use eSearch to check whether <accession> is associated with <gene_hint>.
    Returns (search_validated: bool, message: str).
    """
    # Nucleotide search
    term = f'"{accession}"[Accession] AND "{gene_hint}"[Gene Name]'
    try:
        h   = Entrez.esearch(db='nucleotide', term=term, retmax=1)
        res = Entrez.read(h)
        h.close()
        time.sleep(ENTREZ_DELAY)
        if int(res.get('Count', 0)) > 0:
            return True, f"eSearch confirmed: '{gene_hint}' found for accession '{accession}'"
    except Exception:
        pass

    # Protein search fallback
    try:
        h   = Entrez.esearch(db='protein', term=term, retmax=1)
        res = Entrez.read(h)
        h.close()
        time.sleep(ENTREZ_DELAY)
        if int(res.get('Count', 0)) > 0:
            return True, f"eSearch confirmed (protein db): '{gene_hint}' found for '{accession}'"
    except Exception:
        pass

    return False, (f"eSearch: no records matching gene '{gene_hint}' + "
                   f"accession '{accession}' — verify input")


def fetch_and_validate(gene_name: str, accession: str) -> dict | None:
    """
    Fetch the CDS for accession and validate that gene_name matches the record.

    Returns a dict with keys:
      cds_nt, protein, species, acc_used, back_translated,
      search_validated, record_validated, validation_msg
    or None on complete failure.
    """
    gene_hint = re.split(r'[_\-]', gene_name)[0]   # e.g. XYL1 from XYL1_K270R

    # ── Step A: eSearch validation ────────────────────────────────────────────
    search_ok, search_msg = _esearch_validate(gene_hint, accession)

    # ── Step B: Fetch record ──────────────────────────────────────────────────
    is_prot = _is_protein_accession(accession)

    if not is_prot:
        gb = _ncbi_fetch('nucleotide', accession, 'gb')
        if gb:
            result = _extract_from_gb(gb, gene_hint)
            if result:
                prot = result['protein'] or translate(result['cds_nt'])
                return dict(
                    cds_nt=result['cds_nt'],
                    protein=prot,
                    species=result['species'],
                    acc_used=accession,
                    back_translated=False,
                    search_validated=search_ok,
                    record_validated=result['validated'],
                    validation_msg=search_msg,
                )

    # ── Step C: Protein accession ─────────────────────────────────────────────
    prot_fasta = _ncbi_fetch('protein', accession, 'fasta')
    if prot_fasta:
        recs = list(SeqIO.parse(StringIO(prot_fasta), 'fasta'))
        if recs and len(recs[0].seq) > 30:
            aa   = str(recs[0].seq)
            desc = recs[0].description.lower()
            rec_ok = gene_hint.lower() in desc

            # Try to find a linked nucleotide CDS
            try:
                link_h    = Entrez.elink(dbfrom='protein', db='nucleotide',
                                         id=accession, linkname='protein_nuccore')
                link_data = Entrez.read(link_h)
                link_h.close()
                time.sleep(ENTREZ_DELAY)
                linked = [lk['Id'] for lk in link_data[0]['LinkSetDb'][0]['Link']]
                for nt_id in linked[:3]:
                    gb = _ncbi_fetch('nucleotide', nt_id, 'gb')
                    if gb:
                        res = _extract_from_gb(gb, gene_hint)
                        if res and abs(len(translate(res['cds_nt'])) - len(aa)) <= 2:
                            return dict(
                                cds_nt=res['cds_nt'].upper(),
                                protein=aa,
                                species=res['species'],
                                acc_used=f"{accession}+nt{nt_id}",
                                back_translated=False,
                                search_validated=search_ok,
                                record_validated=rec_ok or res['validated'],
                                validation_msg=search_msg,
                            )
            except Exception:
                pass

            # Back-translate as last resort
            bt      = back_translate(aa)
            sp_m    = re.search(r'\[(.+?)\]', recs[0].description)
            species = sp_m.group(1) if sp_m else 'unknown organism'
            return dict(
                cds_nt=bt,
                protein=aa,
                species=species,
                acc_used=f"{accession}_bt",
                back_translated=True,
                search_validated=search_ok,
                record_validated=rec_ok,
                validation_msg=search_msg,
            )

    return None

# ─────────────────────────────────────────────────────────────────────────────
# Input parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_tsv(path: str) -> list:
    entries = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            name  = (row.get('gene_name') or row.get('name') or '').strip()
            accs  = (row.get('accession') or row.get('accessions') or '').strip()
            mut   = (row.get('mutation') or row.get('mutations') or '').strip()
            notes = (row.get('notes') or row.get('note') or '').strip()
            if not name or not accs:
                continue
            entries.append({
                'gene_name':  name,
                'accessions': [a.strip() for a in accs.split(';') if a.strip()],
                'mutations':  [m.strip() for m in mut.split(';') if m.strip()],
                'notes':      notes,
            })
    return entries

# ─────────────────────────────────────────────────────────────────────────────
# Per-gene pipeline
# ─────────────────────────────────────────────────────────────────────────────

def process_gene(entry: dict, output_dir: str, verbose: bool = True) -> dict:
    gene_name  = entry['gene_name']
    accessions = entry['accessions']
    mutations  = entry['mutations']

    def log(msg):
        if verbose:
            print(msg)

    log(f"\n{'─'*68}")
    log(f"  {gene_name}  [{', '.join(accessions)}]")
    log(f"{'─'*68}")

    result = dict(
        gene_name='', species='', accession='',
        fwd_primer='', rev_primer='',
        codon_optimized_sequence='', pcr_product='',
        status='FAILED', error='',
    )
    result['gene_name'] = gene_name

    # ── 1. Fetch and validate ─────────────────────────────────────────────────
    log("  [1] Fetching and validating sequence...")
    fetch = None
    for acc in accessions:
        log(f"      Trying {acc}...")
        fetch = fetch_and_validate(gene_name, acc)
        if fetch:
            break

    if not fetch:
        result['error'] = f"Could not fetch sequence from: {accessions}"
        log(f"  ERROR: {result['error']}")
        return result

    cds      = fetch['cds_nt']
    species  = fetch['species']
    acc_used = fetch['acc_used']
    bt_flag  = fetch['back_translated']
    result['species']   = species
    result['accession'] = acc_used

    log(f"  Species  : {species}")
    log(f"  Accession: {acc_used}")
    log(f"  CDS      : {len(cds)} nt")
    if bt_flag:
        log("  [!] Protein accession — CDS back-translated from protein sequence")

    # Validation report
    search_ok = fetch['search_validated']
    record_ok = fetch['record_validated']
    if search_ok and record_ok:
        log(f"  Validation: PASS — eSearch and record both confirm gene/accession match ✓")
    elif search_ok or record_ok:
        confirmed_by = "eSearch" if search_ok else "record qualifiers"
        log(f"  Validation: PASS (confirmed by {confirmed_by}) ✓")
        log(f"  {fetch['validation_msg']}")
    else:
        log(f"  Validation: WARNING — gene name not confirmed by eSearch or record.")
        log(f"  {fetch['validation_msg']}")
        log(f"  Verify that accession '{acc_used}' is the correct sequence for '{gene_name}'.")

    # Verify CDS is in-frame
    if len(cds) % 3 != 0:
        result['error'] = f"CDS length {len(cds)} not divisible by 3"
        log(f"  ERROR: {result['error']}")
        return result

    # ── 2. Apply mutations ────────────────────────────────────────────────────
    working    = cds
    mut_annots = []
    if mutations:
        log(f"  [2] Applying {len(mutations)} mutation(s)...")
        for mut in mutations:
            try:
                working, annot = apply_mutation(working, mut)
                mut_annots.append(annot)
                log(f"      {annot}")
            except ValueError as e:
                result['error'] = str(e)
                log(f"  ERROR: {e}")
                return result

    ref_protein = translate(working)

    # ── 3. Codon optimize ─────────────────────────────────────────────────────
    log("  [3] Codon optimizing for S. cerevisiae...")
    try:
        opt_nt, n_changed, opt_score = codon_optimize(working)
    except ValueError as e:
        result['error'] = str(e)
        log(f"  ERROR: {e}")
        return result
    log(f"      Codons changed: {n_changed} | Opt score: {opt_score}%")

    # ── 4. Remove restriction sites ───────────────────────────────────────────
    log("  [4] Scanning for BsaI / BsmBI sites...")
    pre_hits = find_re_sites(opt_nt)
    if pre_hits:
        log(f"      Found {len(pre_hits)} site(s): "
            f"{[(h['enzyme'], h['start']) for h in pre_hits]}")
        try:
            opt_nt, site_edits = remove_re_sites(opt_nt, gene_name)
            for ed in site_edits:
                log(f"      Removed {ed['enzyme']} at {ed['site_pos']}: "
                    f"codon {ed['codon']} {ed['change']}")
        except RuntimeError as e:
            result['error'] = str(e)
            log(f"  ERROR: {e}")
            return result
    else:
        site_edits = []
        log("      No sites found — sequence is clean ✓")

    # ── 5. Verify translation ─────────────────────────────────────────────────
    log("  [5] Verifying translation identity...")
    opt_protein = translate(opt_nt)
    if opt_protein != ref_protein:
        # Locate first mismatch for diagnostic
        for i, (a, b) in enumerate(zip(ref_protein, opt_protein)):
            if a != b:
                result['error'] = (
                    f"Translation mismatch at AA {i+1}: "
                    f"expected {a}, got {b} (context: ...{ref_protein[max(0,i-3):i+4]}...)")
                log(f"  ERROR: {result['error']}")
                return result
        result['error'] = (f"Protein length mismatch: "
                           f"{len(ref_protein)} aa vs {len(opt_protein)} aa")
        log(f"  ERROR: {result['error']}")
        return result

    remaining = find_re_sites(opt_nt)
    if remaining:
        result['error'] = f"Restriction sites remain: {[r['enzyme'] for r in remaining]}"
        log(f"  ERROR: {result['error']}")
        return result

    log(f"      {len(opt_protein)} aa verified, no restriction sites ✓")

    # ── 6. Design Golden Gate primers ─────────────────────────────────────────
    log("  [6] Designing Golden Gate primers...")
    try:
        fwd, rev, pcr = design_primers(opt_nt)
    except ValueError as e:
        result['error'] = str(e)
        log(f"  ERROR: {e}")
        return result

    fwd_bind_len = len(fwd) - len(GG_FWD_TAIL)
    rev_bind_len = len(rev) - len(GG_REV_TAIL)
    log(f"      Fwd ({len(fwd)} nt, Tm bind ≈ {tm_wallace(fwd[len(GG_FWD_TAIL):]):.0f}°C): {fwd}")
    log(f"      Rev ({len(rev)} nt, Tm bind ≈ {tm_wallace(rev[len(GG_REV_TAIL):]):.0f}°C): {rev}")
    log(f"      PCR product: {len(pcr)} nt")
    log(f"      5' overhang: TATG  |  3' overhang: ATCC")

    # ── 7. Write FASTA files ──────────────────────────────────────────────────
    orig_header = (f"{species} | {acc_used} | {len(cds)} nt | original"
                   + (" | back_translated" if bt_flag else ""))
    opt_header  = (f"{species} | {acc_used} | Sc_codon_optimized | "
                   f"codons_changed:{n_changed} | sites_removed:{len(site_edits)} | "
                   f"opt_score:{opt_score}%"
                   + (f" | {'; '.join(mut_annots)}" if mut_annots else "")
                   + (" | back_translated" if bt_flag else ""))

    _write_fasta(os.path.join(output_dir, f"{gene_name}_original.fasta"),
                 f"{gene_name}_original", orig_header, cds)
    _write_fasta(os.path.join(output_dir, f"{gene_name}_optimized.fasta"),
                 f"{gene_name}_optimized", opt_header, opt_nt)

    result.update(dict(
        status='OK',
        fwd_primer=fwd,
        rev_primer=rev,
        codon_optimized_sequence=opt_nt,
        pcr_product=pcr,
        # metadata for report
        orig_nt_len=len(cds),
        opt_nt_len=len(opt_nt),
        protein_len=len(opt_protein),
        codons_changed=n_changed,
        opt_score=opt_score,
        sites_removed=len(site_edits),
        site_edits=site_edits,
        mut_annots=mut_annots,
        back_translated=bt_flag,
        fwd_tm=round(tm_wallace(fwd[len(GG_FWD_TAIL):]), 1),
        rev_tm=round(tm_wallace(rev[len(GG_REV_TAIL):]), 1),
        search_validated=fetch['search_validated'],
        record_validated=fetch['record_validated'],
        validation_msg=fetch['validation_msg'],
    ))
    return result


def _write_fasta(path: str, seq_id: str, description: str, sequence: str):
    rec = SeqRecord(Seq(sequence), id=seq_id, description=description)
    with open(path, 'w') as f:
        SeqIO.write(rec, f, 'fasta')

# ─────────────────────────────────────────────────────────────────────────────
# Output
# ─────────────────────────────────────────────────────────────────────────────

TSV_FIELDS = [
    'gene_name', 'species', 'accession',
    'fwd_primer', 'rev_primer',
    'codon_optimized_sequence', 'pcr_product',
]


def write_summary_tsv(results: list, path: str):
    ok = [r for r in results if r['status'] == 'OK']
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=TSV_FIELDS, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(ok)
    print(f"\nOutput TSV : {os.path.abspath(path)}  ({len(ok)}/{len(results)} genes OK)")

def write_report(results: list, path: str, input_file: str):
    """Write a human-readable pipeline report modelled on pipeline_report.txt."""
    date_str = datetime.date.today().isoformat()
    n_ok     = sum(1 for r in results if r['status'] == 'OK')
    sep      = '=' * 70
    dash     = '-' * 70

    with open(path, 'w') as f:
        def w(line=''):
            f.write(line + '\n')

        # ── Header ────────────────────────────────────────────────────────────
        w('GOLDEN GATE PARTS PREPARATION PIPELINE REPORT')
        w(sep)
        w(f'Generated  : {date_str}')
        w(f'Input file : {os.path.basename(input_file)}')
        w(f'Genes      : {n_ok}/{len(results)} processed successfully')
        w(f'Codon usage: Kazusa S. cerevisiae (taxid 4932)')
        w(f'RE scanned : BsaI (GGTCTC/GAGACC), BsmBI (CGTCTC/GAGACG)')
        w(f'GG format  : 5\'-CGGTCTCAT·[GOI]·ATCCTGAGACCG-3\'')
        w(f'Overhangs  : 5\' TATG (incl. ATG start)  |  3\' ATCC')
        w()

        # ── Summary table ─────────────────────────────────────────────────────
        w('SUMMARY TABLE')
        w(dash)
        w(f"{'Gene':<16} {'Accession':<24} {'AA':>4} {'NT':>6} "
          f"{'Chgd':>5} {'Sites':>5} {'Score':>7}  Status")
        w(dash)
        for r in results:
            if r['status'] == 'OK':
                w(f"{r['gene_name']:<16} {r['accession']:<24} "
                  f"{r['protein_len']:>4} {r['opt_nt_len']:>6} "
                  f"{r['codons_changed']:>5} {r['sites_removed']:>5} "
                  f"{r['opt_score']:>6.1f}%  OK")
            else:
                w(f"{r['gene_name']:<16} {'–':<24} {'–':>4} {'–':>6} "
                  f"{'–':>5} {'–':>5} {'–':>7}  FAILED")
        w()

        # ── Detailed results ──────────────────────────────────────────────────
        w('DETAILED RESULTS')
        w(sep)
        for r in results:
            w()
            w(f"Gene: {r['gene_name']}")
            w(f"  Species   : {r['species']}")
            w(f"  Accession : {r['accession']}")

            if r['status'] != 'OK':
                w(f"  Status    : FAILED — {r['error']}")
                continue

            w(f"  Status    : OK")
            w(f"  Source CDS: {r['orig_nt_len']} nt"
              + (" (back-translated from protein)" if r['back_translated'] else ""))
            w(f"  Optimized : {r['opt_nt_len']} nt, {r['protein_len']} aa")
            w(f"  Opt score : {r['opt_score']}%  |  Codons changed: {r['codons_changed']}")
            w(f"  RE sites  : {r['sites_removed']} removed (0 remain)")

            # Mutations
            if r['mut_annots']:
                w(f"  Mutations : {'; '.join(r['mut_annots'])}")

            # Site-removal detail
            if r['site_edits']:
                w(f"  Site edits:")
                for ed in r['site_edits']:
                    w(f"    {ed['enzyme']} at {ed['site_pos']} — "
                      f"codon {ed['codon']}: {ed['change']}")

            # Validation
            search_ok = r.get('search_validated', False)
            record_ok = r.get('record_validated', False)
            if search_ok and record_ok:
                val_status = 'PASS (eSearch + record)'
            elif search_ok or record_ok:
                val_status = f"PASS ({'eSearch' if search_ok else 'record qualifiers'})"
            else:
                val_status = 'WARNING — not confirmed; verify accession'
            w(f"  Validation: {val_status}")

            # Primers
            w(f"  Fwd primer: {r['fwd_primer']}")
            w(f"    ({len(r['fwd_primer'])} nt | binding Tm ≈ {r['fwd_tm']}°C | "
              f"overhang: TATG)")
            w(f"  Rev primer: {r['rev_primer']}")
            w(f"    ({len(r['rev_primer'])} nt | binding Tm ≈ {r['rev_tm']}°C | "
              f"overhang: ATCC)")
            w(f"  PCR product: {len(r['pcr_product'])} nt")
        w()

        # ── Notes ─────────────────────────────────────────────────────────────
        w('NOTES')
        w(sep)
        w()
        w('Codon optimization:')
        w('  For each amino acid the highest-frequency S. cerevisiae codon was')
        w('  selected (Kazusa taxid 4932, relative frequency threshold ≥ 10%).')
        w('  Stop codons use TAA (most common in yeast). The ATG start codon is')
        w('  always preserved. BsaI/BsmBI sites are removed by synonymous')
        w('  substitution choosing the next-highest-frequency alternative.')
        w()
        w('Primer design:')
        w('  Forward: 5\'-CGGTCTCAT-[binding]-3\'')
        w('    CGGTCTC = BsaI recognition | A = spacer | T completes TATG overhang.')
        w('    Binding starts at ATG (codon 1). Length adjusted from 18 to 30 nt')
        w('    until binding-region Tm ≥ 60°C (Wallace rule).')
        w('  Reverse: 5\'-CGGTCTCAGGAT-[binding_RC]-3\'')
        w('    GGAT = RC of ATCC overhang. Binding is RC of 3\' end of CDS')
        w('    (includes stop codon). Same Tm target.')
        w('  PCR product top strand: CGGTCTCATATG-[GOI]-ATCCTGAGACCG')
        w('  After BsaI digestion: TATG 5\' overhang (5\') and ATCC 5\' overhang (3\').')
        w()

        # Gene-specific notes for back-translated sequences
        bt_genes = [r for r in results if r['status'] == 'OK' and r['back_translated']]
        if bt_genes:
            w('Back-translated genes:')
            for r in bt_genes:
                w(f"  {r['gene_name']}: native bacterial/non-standard CDS not used.")
                w(f"    Protein fetched from {r['accession'].split('_bt')[0]}; "
                  f"CDS back-translated using best yeast codons.")
                w(f"    Original and optimized FASTA are identical (0 codon changes).")
            w()

        mut_genes = [r for r in results if r['status'] == 'OK' and r['mut_annots']]
        if mut_genes:
            w('Point mutations applied:')
            for r in mut_genes:
                for annot in r['mut_annots']:
                    w(f"  {r['gene_name']}: {annot}")
            w()

    print(f"Report     : {os.path.abspath(path)}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Fetch → codon-optimize → GG primer design pipeline.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input',
                        help='Input TSV (gene_name, accession, [mutation], [notes])')
    parser.add_argument('--output', '-o', default=OUTPUT_DIR,
                        help=f'Output directory for FASTA files (default: {OUTPUT_DIR})')
    parser.add_argument('--tsv', default=None,
                        help='Output TSV path (default: <output>/gg_parts_summary.tsv)')
    parser.add_argument('--email', default=DEFAULT_EMAIL,
                        help='NCBI Entrez email address')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress per-step output')
    args = parser.parse_args()

    Entrez.email = args.email
    os.makedirs(args.output, exist_ok=True)
    tsv_out = args.tsv or os.path.join(args.output, 'gg_parts_summary.tsv')

    genes = parse_tsv(args.input)
    if not genes:
        print("No genes found in input TSV.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(genes)} gene(s) from '{args.input}'...")

    results = []
    for entry in genes:
        r = process_gene(entry, args.output, verbose=not args.quiet)
        results.append(r)

    # ── Summary table ─────────────────────────────────────────────────────────
    print(f"\n{'═'*80}")
    print("SUMMARY")
    print(f"{'═'*80}")
    print(f"{'Gene':<20} {'Species':<32} {'AA':>5}  {'Primers':>7}  Status")
    print(f"{'─'*80}")
    for r in results:
        if r['status'] == 'OK':
            sp  = r['species'][:30] if len(r['species']) > 30 else r['species']
            aa  = len(r['codon_optimized_sequence']) // 3 - 1
            plen = f"{len(r['fwd_primer'])}/{len(r['rev_primer'])}"
            print(f"{r['gene_name']:<20} {sp:<32} {aa:>5}  {plen:>7}  OK")
        else:
            print(f"{r['gene_name']:<20} {'–':<32} {'–':>5}  {'–':>7}  "
                  f"FAILED: {r['error']}")

    report_out = os.path.join(args.output, 'pipeline_report.txt')
    write_summary_tsv(results, tsv_out)
    write_report(results, report_out, args.input)
    print(f"FASTAs     : {os.path.abspath(args.output)}/")

    n_ok = sum(1 for r in results if r['status'] == 'OK')
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == '__main__':
    main()
