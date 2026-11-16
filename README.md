# Feather and snakeskin cfDNA fragmentomics

Code for Sacco et al., *CpG Methylation and Gene Expression Inference from DNA
Fragmentation in Feathers and Snakeskin*. Bams: SRA BioProject PRJNA1482836.

Every script starts from duplicate-removed BAM files (reads adapter-trimmed and merged
with SeqPrep2, aligned with `bwa aln`, duplicates removed with Picard MarkDuplicates; see
the STAR Methods). Sample sheets list BAM file names only; edit them to point at your
copies.

| Folder | Makes |
|---|---|
| `01_fragment_lengths/` | Figure 1A, 1B |
| `02_cpg_termini/` | Figure 1C, Figure S2 |
| `03_expression_septiles/` | Figure 2, Table S2 |
| `04_gsea/` | Figure 3, Figure S3, Figure S4, Table S4 |
| `05_ctcf/` | Figure S1 |
| `06_embryonic_rnaseq/` | Figure S5 |

## Environment

```bash
pixi install          # Python, samtools, bedtools, R; versions pinned in pixi.lock
pixi shell
```
FIMO (MEME suite 5.5.9) is also needed for `05_ctcf/run_fimo.sh`. The R annotation
packages used by `04_gsea/export_chicken_go.R` fetch their data in a post-install script,
which pixi only runs after `pixi config set --local run-post-link-scripts insecure`; this
is needed only to regenerate `GOBP_chicken.gmt`, which is included.

## External inputs

| Input | Source |
|---|---|
| `galGal6.fa`, galGal4 `refGene.txt.gz`, galGal6 `rmsk` and `cpgIslandExt` tables, thaSir1 `cpgIslandExt` | UCSC Genome Browser |
| `GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna` (thaSir1) | NCBI |
| `thaSir.genome.fa.out` (snake CR1 annotation) | Gable et al. 2024, *Genome Biol Evol* |
| `Galgal6_refseq.csv` (galGal6 RefSeq genes, refGene layout) | UCSC Table Browser |
| `Gallus_gallus.GRCg6a.104.gtf` | Ensembl release 104 |
| `GSE146956_expression_all_genes.fpkm.txt` (feather follicle FPKM, Zheng et al. 2020) | GEO GSE146956 |
| Embryonic skin RNA-seq coverage bigWigs (Liang et al.) | GEO GSE136224 (URLs in `06_embryonic_rnaseq/samples.tsv`) |

## Running

**Fragment lengths** (`01_fragment_lengths/`)
```bash
python fragment_lengths.py samples.tsv
```

**CpG termini and CG/NCG by compartment** (`02_cpg_termini/`)
```bash
bash make_compartment_beds.sh                       # 12-bp CpG windows per compartment
bash count_termini.sh SCS_DF1.rmdup.bam galGal6.fa fragkon/SCS_DF1 \
    chicken.whole_genome.bed chicken.CR1.bed chicken.CpG_island.bed   # once per library
python cg_ncg_compartments.py samples.tsv           # Fig 1C, Fig S2
```

**Expression septiles** (`03_expression_septiles/`)
```bash
python septiles.py --gtf Gallus_gallus.GRCg6a.104.gtf --fpkm GSE146956_expression_all_genes.fpkm.txt \
    --fragkons ../02_cpg_termini/fragkon/SCS_*.whole_genome.fragkon --bams SCS_*.rmdup.bam
```

**GSEA** (`04_gsea/`)
```bash
python gsea.py --bams feather_bams.tsv --genes Galgal6_refseq.csv --gmt GOBP_chicken.gmt \
    --fpkm GSE146956_expression_all_genes.fpkm.txt     # writes gsea_out/
```

**CTCF motifs** (`05_ctcf/`)
```bash
bash run_fimo.sh galGal6.fa GCF_001077635.1_Thamnophis_sirtalis-6.0_genomic.fna
python ctcf_profiles.py --feather-fimo fimo_out/fimo.tsv --feather-bams feather_bams.tsv \
    --snake-fimo fimo_snake.tsv --snake-bams snake_bams.tsv
```

**Embryonic skin RNA-seq** (`06_embryonic_rnaseq/`; needs the GSEA rankings)
```bash
python correlogram.py --refgene refGene.txt.gz --bigwigs <bigwig_dir> \
    --down ../04_gsea/gsea_out/ranking_down.tsv --flight ../04_gsea/gsea_out/ranking_flight.tsv
```
