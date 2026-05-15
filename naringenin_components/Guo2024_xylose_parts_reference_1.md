# Guo et al. 2024 — Xylose-responsive parts (verified from Supplementary Table S4)

## Source
Guo S, Du J, Li D, Xiong J, Chen Y. (2025). Versatile xylose and arabinose genetic
switches development for yeasts. *Metabolic Engineering* 87:21–36.
DOI: 10.1016/j.ymben.2024.11.004

Sequences extracted directly from the supplementary docx file (Table S4),
preserving the authors' color annotations.

## Author's annotation key (from Table S4 footer)
- **Blue** (color hex `0070C0`)  — TATA-box
- **Red** (color hex `FF0000`)   — binding sites of regulators (xlnO / xylO operators)
- **Green** (color hex `00B050`) — transcription start site (TSS) region
- **Bold**                       — kozak sequence (`TAAATAAAA` preceding ATG)
- **lower case**                 — Golden Gate scars (`tttt` at 5′, `aatg` at 3′ as kozak-start)
- **underline**                  — core promoter

## Important note from the inspection
The authors color-coded the **xlnO** operators (XlnR binding sites) in red consistently
in all xylose-responsive promoters. But in the dual-mode promoters (Pxln.2b-xyl,
Pxln.4-xyl), the **xylO** operators (XylR binding sites, transferred from Pxyl) appear
as black text between the TATA and TSS, NOT separately colored red.

The xylO operator consensus is `AGTTAGTTTGTTTAAACAACAAACT` (25 bp, an inverted repeat),
present as a **tandem pair** in Pxyl and inherited as a tandem pair in Pxln.2b-xyl
and Pxln.4-xyl. The colored OPERATOR tag below only marks the xlnO sites; the xylO
positions are noted with a `[xylO:...]` annotation I added by sequence matching to
the Pxyl consensus.

---

## Architecture summary

| Promoter      | xlnO sites | xylO sites | TATA       | TSS region              | Length |
|---------------|-----------|-----------|------------|-------------------------|--------|
| Pxln.2b       | 2 (red)   | 0         | TATAAATA   | AGAATATCAAGCTACAAAAA    | 382 nt |
| Pxyl          | 0         | 2         | TATAAAA    | AGAATATCAAGCTACAAAAA    | 177 nt |
| **Pxln.2b-xyl** | **2**   | **2**     | **TATAAATA** | **AGAATATCAAGCTACAAAAA** | **377 nt** |
| Pxln.4        | 4 (red)   | 0         | TATAAATA   | AGAATATCAAGCTACAAAAA    | 422 nt |
| Pxln.4-xyl    | 4         | 2         | TATAAATA   | AGAATATCAAGCTACAAAAA    | 417 nt |
| PxlnD         | (native)  | 0         | (native)   | (native)                | 588 nt |
| Pxylp.a       | (native)  | 0         | (native)   | (native)                | 729 nt |

Construction logic for the dual-mode (Pxln.X-xyl) promoters: take the upstream
UAS+xlnO region of Pxln.X, replace the downstream core promoter region (between
TATA and TSS) with the xylO-containing region from Pxyl. This places XylR
operators **between TATA and TSS**, which is exactly the position where
LacI-family repressors block transcription initiation.

---

## Sequences

All sequences are oriented 5′ → 3′ and include the YTK-compatible Golden Gate
flanking overhangs: `tttt` at the 5′ end, and `aatg` (the start of the kozak/ATG
fusion site) at the 3′ end. The annotated coding sequences for the two
transcription factors are also included, formatted the same way for direct
drop-in as YTK Type 3a parts.

### Operator consensus (for primer design / verification)

- **xlnO (XlnR binding site)**: `GAATTTAGGCTAAAGAAAGATC` (22 bp; appears multiple
  times in tandem, separated by short spacers, upstream of the TATA-box)
- **xylO (XylR binding site)**: `AGTTAGTTTGTTTAAACAACAAACT` (25 bp inverted
  repeat; appears as a tandem pair between TATA and TSS in Pxyl, Pxln.2b-xyl,
  and Pxln.4-xyl)


===========================================================================
>Pxln.2b-xyl  [PRIMARY DUAL-MODE PROMOTER — chassis default]
===========================================================================
Length: 377 nt
Architecture: 2× xlnO (XlnR-binding, activator gate) upstream of TATA
              + 2× xylO (XylR-binding, repressor gate) between TATA and TSS
Behavior: Both gates must be open for transcription.
          - On glucose: XlnR inactive (no recruitment) AND XylR bound to xylO
            (repression). Both gates closed. Promoter silent.
          - On xylose: XlnR-P recruits Pol II AND XylR released from xylO.
            Both gates open. >4000-fold induction per Guo 2024 Fig 3.

## Plain sequence (lowercase = Golden Gate scars):

ttttGGGATCCTGAGAATTTAGGCTAAAGAAAGATCGCTAATAGCCGACAAACGGACTCCAAAAAAAAAAA
GGAGGCACAGGACAAACGCAGCACCTGCGTCATTCACGACACCGGGCAGAATTTAGGCTAAAGAAAGATCA
ATAAGAGAATTTCAGATTGAGAGAATGAAAAAAAAAAAAAAAAAAAAGGCAGTCTTGGACAATAGTTATAT
CACGTTGACTCAGGAAGTGGTGGTGGAGACTGCCACATATAAATAAGTTAGTTTGTTTAAACAACAAACTT
TTTTCATTTCTTTTGTTTCCCCTTCTCTTCTTTTAGTTAGTTTGTTTAAACAACAAACTAGAATATCAAGC
TACAAAAATAAATAAAAAaatg

## Feature breakdown (5′ → 3′):

  [Golden Gate 5' scar]   tttt
  [linker]                GGGATCCTGA
  [xlnO #1]               GAATTTAGGCTAAAGAAAGATC
  [spacer]                GCTAATAGCCGACAAACGGACTCCAAAAAAAAAAAGGAGGCACA
                          GGACAAACGCAGCACCTGCGTCATTCACGACACCGGGCA
  [xlnO #2]               GAATTTAGGCTAAAGAAAGATC
  [UAS/intergenic]        AATAAGAGAATTTCAGATTGAGAGAATGAAAAAAAAAAAAAAAAAAAA
                          GGCAGTCTTGGACAATAGTTATATCACGTTGACTCAGGAAGTGGTGGTGGA
                          GACTGCCACA
  [TATA box]              TATAAATA
  [xylO #1]               AGTTAGTTTGTTTAAACAACAAACT
  [spacer]                TTTTTCATTTCTTTTGTTTCCCCTTCTCTTCTTTT
  [xylO #2]               AGTTAGTTTGTTTAAACAACAAACT
  [TSS region]            AGAATATCAAGCTACAAAAA
  [kozak (Type 3a scar)]  TAAATAAAAA
  [Golden Gate 3' scar]   aatg


===========================================================================
>Pxln.2b  [XlnR-only single-mode — for activator-only architecture]
===========================================================================
Length: 382 nt
Architecture: 2× xlnO upstream of TATA. No xylO. Higher leak on glucose but
              simpler regulatory architecture than the dual-mode.

## Plain sequence:

ttttGGGATCCTGAGAATTTAGGCTAAAGAAAGATCGCTAATAGCCGACAAACGGACTCCAAAAAAAAAAA
GGAGGCACAGGACAAACGCAGCACCTGCGTCATTCACGACACCGGGCAGAATTTAGGCTAAAGAAAGATCA
ATAAGAGAATTTCAGATTGAGAGAATGAAAAAAAAAAAAAAAAAAAAGGCAGTCTTGGACAATAGTTATAT
CACGTTGACTCAGGAAGTGGTGGTGGAGACTGCCACATATAAATAGAGTGCCAGTAGCGACTTTTTTCACA
CTCGAAATACTCTTACTACTGCTCTCTTGTTGTTTTTATCACTTCTTGTTTCTTCTTGGTAAATAGAATAT
CAAGCTACAAAAATAAATAAAAAaatg

## Feature breakdown (note this is identical to Pxln.2b-xyl up to and
   including TATA; differs in the post-TATA region which has no xylO):

  [Golden Gate 5' scar]   tttt
  [linker]                GGGATCCTGA
  [xlnO #1]               GAATTTAGGCTAAAGAAAGATC
  [spacer 1]              GCTAATAGCCGACAAACGGACTCCAAAAAAAAAAAGGAGGCACA
                          GGACAAACGCAGCACCTGCGTCATTCACGACACCGGGCA
  [xlnO #2]               GAATTTAGGCTAAAGAAAGATC
  [UAS/intergenic]        AATAAGAGAATTTCAGATTGAGAGAATGAAAAAAAAAAAAAAAAAAAA
                          GGCAGTCTTGGACAATAGTTATATCACGTTGACTCAGGAAGTGGTGGTGGA
                          GACTGCCACA
  [TATA box]              TATAAATA
  [core/spacer]           GAGTGCCAGTAGCGACTTTTTTCACACTCGAAATACTCTTACTACTGCTCT
                          CTTGTTGTTTTTATCACTTCTTGTTTCTTCTTGGTAAAT
  [TSS region]            AGAATATCAAGCTACAAAAA
  [kozak]                 TAAATAAAAA
  [Golden Gate 3' scar]   aatg


===========================================================================
>Pxyl  [XylR-only single-mode — for repressor-only architecture]
===========================================================================
Length: 177 nt
Architecture: No xlnO. 2× xylO between TATA and TSS. The shortest of the
              xylose-responsive set; weakest absolute output but cleanest OFF.

## Plain sequence:

ttttCAACGGCCTAGCATGTGATTAATTAATTATTTTGTTTTTTTTTTGCAGTATAAAAAGTTAGTTTGTT
TAAACAACAAACTTTTTTCATTTCTTTTGTTTCCCCTTCTCTTCTTTTAGTTAGTTTGTTTAAACAACAAA
CTAGAATATCAAGCTACAAAAATAAATAAAAaatg

## Feature breakdown:

  [Golden Gate 5' scar]   tttt
  [linker]                CAACGGCCTAGCATGTGATTAATTAATTATTTTGTTTTTTTTTTGCAG
  [TATA box]              TATAAAA
  [xylO #1]               AGTTAGTTTGTTTAAACAACAAACT
  [spacer]                TTTTTCATTTCTTTTGTTTCCCCTTCTCTTCTTTT
  [xylO #2]               AGTTAGTTTGTTTAAACAACAAACT
  [TSS region]            AGAATATCAAGCTACAAAAA
  [kozak]                 TAAATAAAA
  [Golden Gate 3' scar]   aatg


===========================================================================
>Pxln.4  [stronger XlnR-only variant with 4× xlnO]
===========================================================================
Length: 422 nt (extracted from supp)
Architecture: 4× xlnO upstream of TATA. Higher max output than Pxln.2b but
              more leakiness on glucose. Used in P. pastoris (per Fig S17/S18).


===========================================================================
>Pxln.4-xyl  [dual-mode 4× xlnO + 2× xylO]
===========================================================================
Length: 417 nt
Architecture: 4× xlnO upstream of TATA + 2× xylO between TATA and TSS.
              Higher max output than Pxln.2b-xyl. Roughly equivalent dynamic
              range (4000+ fold). Choice between Pxln.2b-xyl and Pxln.4-xyl
              depends on whether basal expression (lower in 2b-xyl) or max
              output (higher in 4-xyl) matters more for the target gene.


===========================================================================
>XlnR (Aspergillus nidulans) — codon-optimized CDS
===========================================================================
Length: 2633 nt (including aatg 5' scar and taaa 3' scar)
Protein: 875 aa (after removing scars)
Encodes: XlnR transcriptional activator (Zn2Cys6 zinc cluster family)
Function: Activated by xylose binding (phosphorylation in native context);
          recruits RNA Pol II at the xlnO operator.

## Plain CDS (lowercase = YTK Type 3a Golden Gate scars):

aatgTCGCAATCCCAGTCTCAGACGATCGGGCTTGACACCCTCGCCGAGGGCTCGCAATATGTGCTAGAGC
AGCTGCAGTTATCGCGAGAGGGCGGCAACTCTGAGAACAACTCTACTTTCAAGCCATCCTCCGTCCGCGAC
TCTTTAGCTGAAGCCCGTTCCATGATCCGCAAGAACTCTTCGTCAGCGCCTGTCCGCCGGAGAATCAGTCG
GGCTTGTGACCAATGTAACCAACTCCGGACAAAGTGCGACGGGCAGAATCCGTGCGCGCATTGTATAGAAT
TCGGTTTAACATGCGAATACGCTCGGGAGCGGAAGAAACGGGGCAAGGCTTCAAAAAAAGATATTGCTGCT
GCTGCTGCCGCTGCAGGACATCAGGGAGGCATGGGTAACCGATCACCTACAGATAGACGCCTGTCGCAGGA
GCCAGGCGGCCGGTACGATTCCGTTCTTGAGGCATCGCGCGTTCAATCACATCTACCTGCGAACGGTTTGT
CTAGTATTCACAACACTCAAGCGGCGCACTCGCAGCCGCCGCTGGGGTCGGCCCTTGATGCTCTACATTTG
AACCATTTTACTCAGCTGAACGAGTCGGGCCGCTCCCAGATGCCCGTGTCAGACCTTCGATCACTCCAAAT
CCTCCACAACAACCCTCGCTCTCCATCCGCTCTTCCGCACGGCCTAAACGCCTATAATGATAACACATTCT
CGTTGCTGAACTCGCAAGAGCCGAATACGACTTCACTCAATCACTTCCGACTCGGAAACTCCACGGATAAC
CCGTCGGCTCAGTTTTTAGGCCTCTCACCTCCCGCTCAATCCCCAGGATGGCTTCCGCTGCCATCACCGTC
GCCCGCAAACTTCCCTTCATTTCCCATGGCTCCGTTCTCTGGAACCAGTCTACGTTACCCCGTCCTTCAGC
CGGTTCTTCCACATATCGCCTCAATAATTCCTCAGTCTCTCGCTTGCGACCTTCTTGACCTTTATTTTACG
AGCTCCTCCTCTTCACACCTATCACCCCAGTCTCCGTACGTTGTTGGTTATATCTTTCGCAAACAGTCGTT
CCTCCACCCAACAAAGCCGCGCGTGTGCTCGCCTGGGCTATTAGCGAGTATGCTCTGGGTAGGCGCACAAA
CAAGCGATGCGCCGTTCCTGACGTCTCCGCCCTCAGCGCGCGGTCGGGTATGTCAGAAGCTACTAGAATTG
ACGATAGGGCTGCTGCGCCCGCTCATTCACGGACCGGCGCTCGGAGAAGCCTCGCCAAATTATGCCGCAAA
TATGGTAATAAACGGCGTCGCGCTCGGTGGTTTTGGTGTCTCGATGGACCAACTAGGAGCCCAAAGTACTG
CGACAGGAGCTGTGGATGATGTCGCCACGTACGTCCATCTCGCTACGGTAGTATCTGCCAGTGAATACAAA
GCCGCAAGTATGCGCTGGTGGACGGCCGCCTGCTCACTCGCCCGAGAGCTGAAACTTGGCCGCGAGCTGCC
CCCCAACGCATCGCAACCGGGTCAGGACGGCGAGCGAGAAAACGAGGGCGACAATCCATCAAAGCGAAACC
AGTCATTGCACGGCGGAAACTCTAATGTCAACGTCACAGAAGAAGAACGAGAGGAGCGCCGGCGTCTGTGG
TGGCTACTGTACGCTACCGATCGCCATCTGGCTTTATGCTATAATAGGCCGCTTACTCTGTTGGATAAAGA
GTGTTCGCAGCTACTGCAGCCCATGAATGATGATTTATGGCAGGCAGGAGATTTTCCAGCTGCGACCTACC
GCGCTGTCGGCCCGCCCATTGAATGCACAGGCCACAGCATGTTCGGCTACTTCCTGCCATTGATGACGATA
CTTGGGGGAATTATCGATCTCCAACAAGCGCGAGAACATCCACGGTACGGACTTACCTTCCGCAGCGGCCC
CGATCTAGATCAGTACATCATGGCCATAACCCAGCAACTTGACGCCTACGGGCAGAGCCTAAAAGACTTCG
AAGCACGATATATAAATAGCCTCGCCCTAGCAGAGAATGAACCGCCCGAGAACCCGCACATTGATCACCTC
AGCCCATCCGGCCGGTCCAGCAGTACGGTTGGCTCGCGCGTCAACGAGTCAATCGTCCATACTAAGATGGT
AGTCGCCTACGGCACCCACATCATGCACGTCTTGTATGTTCTCCTAGCGGGTAAATGGGACCCCATCAACC
TCCTTGAGGACCATGATATGTGGATCTCTTCCGAGTCATTCCTCGCGGCCATGAGCCACGCTGTCGGCGCC
GCAGAGGCAGCAGCCGATATTCTCGAGTACGATCCGGATTTGAGCTTCATGCCGTTCTTTTTCGGCATTTA
CCTGCTCCAGGGAAGTTTTCTACTCTTGCTTGCTGCAGATAAGCTACAGGGGGATGCGAATCCGAGCGTCG
TCCGCGCTTGTGAGACTATCGTGCGCGCACATGAGGCTTGTGTCGTTACGCTGAATACAGAGTATCAGCGA
ACATTCCGCAAAGTGATGCGCTCCGCTCTTGCCCAGGTTCGGGGACGCGTTCCAGATGATTTTGGTGAGCA
GCAGCAGCGCCGGCGGGAAGTGCTGTCTCTTTACCGCTGGACCGGTGATGGGACGGGGCTCGCGCTCTCTT
GAtaaa


===========================================================================
>XylR (Bacillus licheniformis) — codon-optimized CDS
===========================================================================
Length: 1196 nt (including aatg 5' scar and taaa 3' scar)
Protein: 397 aa (after removing scars)
Encodes: XylR transcriptional repressor (LacI family).
Function: Binds the xylO operator in absence of xylose. Xylose binding causes
          allosteric release from the operator, derepressing the promoter.
          Note this is the B. licheniformis version (used by Guo 2024),
          NOT the C. crescentus version (used by Vanee 2017 / US 9,506,097).

## Plain CDS (lowercase = YTK Type 3a Golden Gate scars):

aatgAACACCGCTGATCAGGCTCTTGTTAAAAAAATGAATAAGGCGCTAATCTTTGAGCAAATCATTGAGA
ACGGCCCCGTCTCCAGGGCAAAACTAAGCGAGATAACCGGATTAAATAAGGCTACCGTCTCAAGCCAGGTA
AGCAGTCTACTTCAGAAGGATATAATTTACGAAACCGGCCCTGGCGAAAGTTCTGGGGGAAGGAGGCCTGT
TATGCTGAAATTTAATCGTAAAGCCGGCTACGCCGTAGGTGTAGATGTGGGAACTAATTACATTATCGTCG
CTCTAACCGACTTGGAAGGCCATTTAATAGAACAATTCGAAAGGACATTAGATGAGGAGGACATCCAGGCT
ACTGAGGAAGCGTTAATTGAGCTGACCGGCCTGGCAGTAGATAAAATACCGCCCTCTCCTTTCGGTCTAAC
AGGTATCGGTGTCTGCGTGCCCGGACTTGTCGATAATGAACGTCATGTCGTGTTTACCCCTAACAAGCCCA
TTCATTTGATTCCAATAAAGGAAAAACTAGAGGAGAGGTTTGGAGTCCCGATCTTAATCGAAAACGAAGCC
AACGCCGGGGCAGTCGGCGAAAAGGAATATGGGGAAGGCGGCCAACTGGAACACGCAGTGTTTGTCTCCAT
TAATACGGGGATAGGATTGGGAATTTTGATGAATGGCAAGTTGTTTAGGGGTGTGCAGGGATTTTCAGGGG
AAGCAGGCCACATGAGTATTCATTTCGATGGCCCACTATGCCGTTGCGGGAACCGTGGTTGTTGGGAGTTG
TACGCGAGTGAGAAGGCCGTTTTTTCTCATTACGCGGCCAACAGCGGCGCGCAACTATATGAAACTGTAAA
GGAATTAGCCGACAGAGGTGACCCCGGCATGATGGAGACGTTCGAACGTTTTGGATTCCATATTGGGATCG
GATTGTTGAATATATTGAAAACGCTGAATCCGGATACAATTATACTGAGAAATACGATTGTTGAGTCATAC
CCTAGCATAGTCGACGCAATAAAGAAAACCATTGCCTCTCGTTCTGCCGCAGAGGCGTTAAGCAATTATCA
CCTTAAAATCAGCACTCTAGGCAGGACGGCATCTGCGCTAGGGATGTCCAGCTTAGTGACGGAGCGTTTCT
TAGAAAGGTTTATGAATGAGCGTTTTGGATCACCAAAAAAGAAGCGTAAGGTGTAAtaaa

Note: The TF column in Table S4 lists this as "XylR -NLS(Bacillus licheniformis)"
when paired with Pxyl. **Verified by in silico translation**: this CDS encodes
the native B. licheniformis XylR (388 aa) C-terminally fused via a GS linker to
the SV40 large T-antigen NLS (PKKKRKV), giving a 397 aa final protein. The
C-terminus reads `...MNERF` (native XylR end) + `GS` (linker) + `PKKKRKV` (SV40
NLS) + stop. The NLS is therefore already encoded in this sequence — no
additional NLS fusion is needed when using this CDS.

By contrast, the XlnR CDS above does NOT have an engineered NLS; as a fungal
Zn2Cys6 transcription factor it has native nuclear targeting and apparently
doesn't need one. Verified by translation: XlnR C-terminus is `...DGTGLALS*`,
no NLS motif.

---

## How to use these in your build

For the YTK/MYT MoClo system:
1. The 5' `tttt` and 3' `aatg` flanks are designed as Golden Gate scars
   compatible with the Lee YTK Type 2 (promoter) overhangs, with the 3' end
   landing exactly at the kozak/ATG fusion site for a downstream Type 3a CDS.
2. The TF CDSs (XlnR, XylR) have `aatg...taaa` flanks — these are Type 3a
   (full ORF) overhangs in YTK convention.
3. Synthesize as gBlocks or pre-cloned in pUC19 backbones. The promoters are
   small enough (≤420 nt) to be ordered as single gBlocks. The XlnR CDS at
   2633 nt may need to be split into two gBlocks or ordered as a clonal gene.
4. BsaI/BsmBI scanning: I have NOT verified these sequences are free of
   internal BsaI/BsmBI sites. Run the standard scan before ordering.

## Citation

When citing these parts in protocols or publications:
Guo S, Du J, Li D, Xiong J, Chen Y. Versatile xylose and arabinose genetic
switches development for yeasts. Metabolic Engineering 2025;87:21-36.
doi:10.1016/j.ymben.2024.11.004. Supplementary Table S4.
