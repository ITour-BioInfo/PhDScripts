#!/usr/bin/env python3
"""
Unified MAPT junction extractor + E10 inclusion calculator + batch orchestrator.

Subcommands:
  extract-juncs   Extract splice junctions from a BAM in a target region into regtools-like BED12.
  compute-e10     Compute MAPT exon-10 inclusion (Ψ, 4R) and 4R/3R ratio from *.MAPT.juncs.bed files.
  batch           For many BAMs: run extract-juncs on each, then compute-e10 once across all outputs.

Dependencies:
  pysam, pandas, gtfparse  (optional: polars; if present, it's converted to pandas transparently)

Examples:
  # 1) Extract junctions
  python3 mapt_unified.py extract-juncs --bam A.bam --region 17:45890000-46050000 --out A.MAPT.juncs.bed

  # 2) Compute Ψ and 4R/3R
  python3 mapt_unified.py compute-e10 MAPT_only.gtf ./juncs --canonical --tol 2 --min-support 20

  # 3) Batch: run both steps over a directory of BAMs, then compute-e10
  python3 mapt_unified.py batch --bam-dir ./bams --gtf MAPT_only.gtf --region 17:45890000-46050000 \
                                --outdir ./out --canonical --tol 2 --min-support 20
"""
from __future__ import annotations
import sys, os, argparse, collections, subprocess
from pathlib import Path

# ---------- Utilities ----------
def eprint(*a, **kw):
    print(*a, file=sys.stderr, **kw)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

# ---------- Subcommand: extract-juncs (from pysam) ----------
def extract_junctions(
    bam: str,
    region: str = "17:45890000-46050000",
    out: str | None = None,
    min_anchor: int = 8,
    min_intron: int = 50,
    max_intron: int = 500_000,
    force_strand: str = ".",
) -> str:
    """
    Parse CIGARs, detect N (splice) ops, enforce anchor sizes and intron span, and emit BED12 with two blocks.
    """
    import pysam  # local import to avoid forcing dependency on unrelated subcommands

    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    sam = pysam.AlignmentFile(bam, "rb")

    counts = collections.Counter()
    for read in sam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.cigartuples:
            continue

        ref_pos = read.reference_start
        ops = read.cigartuples
        for i, (op, lenop) in enumerate(ops):
            if op == 3:  # N (splice)
                j_start = ref_pos
                j_end = ref_pos + lenop
                # left anchor
                left_anchor = 0
                for j in range(i - 1, -1, -1):
                    op_b, len_b = ops[j]
                    if op_b in (0, 7, 8):  # M, =, X
                        left_anchor = len_b
                        break
                # right anchor
                right_anchor = 0
                for j in range(i + 1, len(ops)):
                    op_f, len_f = ops[j]
                    if op_f in (0, 7, 8):
                        right_anchor = len_f
                        break

                intron_len = j_end - j_start
                if (
                    left_anchor >= min_anchor
                    and right_anchor >= min_anchor
                    and (min_intron <= intron_len <= max_intron)
                ):
                    counts[(chrom, j_start, j_end, left_anchor, right_anchor, force_strand)] += 1

            if op in (0, 2, 3, 7, 8):  # ref-consuming ops
                ref_pos += lenop

    sam.close()

    outpath = out or (Path(bam).stem + ".MAPT.juncs.bed")
    with open(outpath, "w") as outfh:
        for (chrom, jstart, jend, left_a, right_a, strand), cnt in sorted(counts.items()):
            chromStart = jstart - left_a
            chromEnd = jend + right_a
            blockCount = 2
            intron_len = jend - jstart
            blockSizes = f"{left_a},{right_a}"
            blockStarts = f"0,{left_a + intron_len}"
            name = f"junc_{jstart}_{jend}"
            outfh.write(
                "\t".join(
                    map(
                        str,
                        [
                            chrom,
                            chromStart,
                            chromEnd,
                            name,
                            cnt,
                            strand,
                            chromStart,
                            chromEnd,
                            0,
                            blockCount,
                            blockSizes,
                            blockStarts,
                        ],
                    )
                )
                + "\n"
            )
    eprint(f"Wrote {outpath}")
    return outpath

# ---------- Subcommand: compute-e10 (Ψ and 4R/3R) ----------
def maybe_to_pandas(df):
    try:
        import polars as pl
        if isinstance(df, pl.DataFrame):
            return df.to_pandas()
    except Exception:
        pass
    return df

def load_mapt_exons(gtf_path: str, use_canonical: bool = False):
    import pandas as pd
    from gtfparse import read_gtf

    gtf = read_gtf(gtf_path)
    gtf = maybe_to_pandas(gtf)
    gtf = gtf[gtf.get("gene_name") == "MAPT"]
    if gtf.empty:
        sys.exit("ERROR: No MAPT rows in GTF. Check build / gene_name / 'chr' prefix.")

    if use_canonical:
        canon_mask = (
            (gtf.get("transcript_name") == "MAPT-201")
            | (gtf.get("tag", "").astype(str).str.contains("Ensembl_canonical", na=False))
            | (gtf.get("tag", "").astype(str).str.contains("MANE", na=False))
        )
        canon = gtf[(gtf["feature"].isin(["gene", "transcript", "exon"])) & canon_mask]
        if not canon.empty:
            gtf = canon

    exons = gtf[gtf["feature"] == "exon"].copy()
    if exons.empty:
        sys.exit("ERROR: MAPT GTF has no exon rows (after filtering).")
    return exons

def infer_exon10_triplet(exons):
    import pandas as pd
    strand = exons["strand"].mode().iloc[0]
    rows = []
    for tx, df in exons.groupby("transcript_id"):
        df = df.sort_values("start" if strand == "+" else "end")
        es = df[["start", "end", "exon_id"]].to_numpy()
        for i in range(len(es) - 1):
            up_s, up_e, _ = es[i]
            dn_s, dn_e, _ = es[i + 1]
            rows.append((int(up_s), int(up_e), int(dn_s), int(dn_e)))
    import pandas as pd
    adj = pd.DataFrame(rows, columns=["up_s", "up_e", "dn_s", "dn_e"])
    edges = adj.value_counts(["up_s", "up_e", "dn_s", "dn_e"]).rename("n").reset_index()

    cands = []
    for _, r in edges.iterrows():
        u_s, u_e, dv_s, dv_e = int(r.up_s), int(r.up_e), int(r.dn_s), int(r.dn_e)
        vw = edges[(edges.up_s == dv_s) & (edges.up_e == dv_e)]
        for _, rr in vw.iterrows():
            w_s, w_e = int(rr.dn_s), int(rr.dn_e)
            skip = edges[
                (edges.up_s == u_s)
                & (edges.up_e == u_e)
                & (edges.dn_s == w_s)
                & (edges.dn_e == w_e)
            ]
            if len(skip):
                exon_len = dv_e - dv_s + 1
                support = int(r.n) + int(rr.n) + int(skip.n.iloc[0])
                cands.append(
                    {
                        "u_e": u_e,
                        "v_s": dv_s,
                        "v_e": dv_e,
                        "w_s": w_s,
                        "exon_len": exon_len,
                        "support": support,
                    }
                )
    c = pd.DataFrame(cands)
    if c.empty:
        sys.exit("ERROR: Could not infer exon-10 cassette from GTF (try --canonical).")
    c["len_score"] = -abs(c["exon_len"] - 93)  # exon 10 ~ 93 bp
    best = c.sort_values(["len_score", "support"], ascending=[False, False]).iloc[0]
    u_e, v_s, v_e, w_s = map(int, [best["u_e"], best["v_s"], best["v_e"], best["w_s"]])
    eprint(f"DEBUG: exon10 coords (u_e={u_e}, v_s={v_s}, v_e={v_e}, w_s={w_s})")
    return u_e, v_s, v_e, w_s

def read_regtools_like_bed(fp: Path):
    import pandas as pd
    rows = []
    with open(fp) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("track"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            chromStart, chromEnd = int(f[1]), int(f[2])
            cnt = int(f[4])
            blockSizes = [int(x) for x in f[10].strip().split(",") if x != ""]
            if len(blockSizes) < 2:
                continue
            j_start = chromStart + blockSizes[0]
            j_end = chromEnd - blockSizes[1]
            rows.append((j_start, j_end, cnt))
    import pandas as pd
    return pd.DataFrame(rows, columns=["j_start", "j_end", "count"])

def sum_near(df, js: int, je: int, tol: int) -> int:
    m = df[(df["j_start"].between(js - tol, js + tol)) & (df["j_end"].between(je - tol, je + tol))]
    return int(m["count"].sum())

def compute_e10(
    gtf: str,
    junc_dir: str,
    canonical: bool = False,
    tol: int = 2,
    min_support: int = 20,
    pseudocount: float = 0.5,
    out_csv: str = "MAPT_3R4R_per_sample.csv",
) -> Path:
    import pandas as pd

    exons = load_mapt_exons(gtf, use_canonical=canonical)
    u_e, v_s, v_e, w_s = infer_exon10_triplet(exons)

    jdir = Path(junc_dir)
    beds = sorted(jdir.glob("*.MAPT.juncs.bed"))
    if not beds:
        sys.exit(f"ERROR: No *.MAPT.juncs.bed files found in {jdir}")

    out_rows = []
    for bed in beds:
        df = read_regtools_like_bed(bed)
        c1 = sum_near(df, u_e, v_s, tol)  # 9->10
        c2 = sum_near(df, v_e, w_s, tol)  # 10->11
        cs = sum_near(df, u_e, w_s, tol)  # 9->11 (skip)
        I = (c1 + c2) / 2.0
        S = float(cs)
        total = I + S
        if total < min_support:
            psi = (I + pseudocount) / (total + 2 * pseudocount)
        else:
            psi = (I / total) if total > 0 else 0.5
        ratio = psi / (1.0 - psi) if psi < 1.0 else float("inf")

        out_rows.append(
            {
                "sample": bed.name.replace(".MAPT.juncs.bed", ""),
                "incl_left": int(c1),
                "incl_right": int(c2),
                "skip": int(cs),
                "psi_4R": round(psi, 6),
                "ratio_4R_over_3R": ("inf" if ratio == float("inf") else round(ratio, 6)),
                "total_support": int(c1 + c2 + cs),
            }
        )

    out = pd.DataFrame(out_rows).sort_values("sample")
    out.to_csv(out_csv, index=False)
    print(out.to_string(index=False))
    eprint(f"Wrote {out_csv}")
    return Path(out_csv)

# ---------- Subcommand: batch (orchestrate) ----------
def batch_pipeline(
    bam_dir: str | None,
    bam_list: str | None,
    gtf: str,
    region: str,
    outdir: str = "out",
    min_anchor: int = 8,
    min_intron: int = 50,
    max_intron: int = 500_000,
    force_strand: str = ".",
    canonical: bool = False,
    tol: int = 2,
    min_support: int = 20,
    pseudocount: float = 0.5,
):
    outdir = Path(outdir)
    junc_dir = outdir / "junctions"
    ensure_dir(junc_dir)

    # Collect BAMs
    bams = []
    if bam_list:
        with open(bam_list) as fh:
            for ln in fh:
                ln = ln.strip()
                if ln and not ln.startswith("#"):
                    bams.append(Path(ln))
    elif bam_dir:
        bams = sorted(Path(bam_dir).glob("*.bam"))
    else:
        sys.exit("Provide --bam-dir or --bam-list")

    if not bams:
        sys.exit("No BAMs found.")

    # Extract junctions for each BAM
    for bam in bams:
        out_bed = junc_dir / (bam.stem + ".MAPT.juncs.bed")
        extract_junctions(
            bam=str(bam),
            region=region,
            out=str(out_bed),
            min_anchor=min_anchor,
            min_intron=min_intron,
            max_intron=max_intron,
            force_strand=force_strand,
        )

    # Compute E10 metrics
    out_csv = outdir / "MAPT_3R4R_per_sample.csv"
    compute_e10(
        gtf=gtf,
        junc_dir=str(junc_dir),
        canonical=canonical,
        tol=tol,
        min_support=min_support,
        pseudocount=pseudocount,
        out_csv=str(out_csv),
    )

# ---------- CLI ----------
def build_parser():
    ap = argparse.ArgumentParser(
        description="Unified MAPT junction extraction and E10 (4R) quantification."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    # extract-juncs
    p1 = sub.add_parser("extract-juncs", help="Extract junctions to regtools-like BED12")
    p1.add_argument("--bam", required=True)
    p1.add_argument("--region", default="17:45890000-46050000")
    p1.add_argument("--out", default=None)
    p1.add_argument("--min-anchor", type=int, default=8)
    p1.add_argument("--min-intron", type=int, default=50)
    p1.add_argument("--max-intron", type=int, default=500000)
    p1.add_argument("--force-strand", choices=["+", "-", "."], default=".")

    # compute-e10
    p2 = sub.add_parser("compute-e10", help="Compute Ψ (4R) and 4R/3R ratio from junction BEDs")
    p2.add_argument("gtf", help="MAPT-only GTF (GRCh38; seqnames consistent with BAMs/BEDs)")
    p2.add_argument(
        "junc_dir",
        help="Directory with *.MAPT.juncs.bed files (from 'extract-juncs' or regtools)",
    )
    p2.add_argument("--canonical", action="store_true", help="Prefer canonical MAPT transcript if present")
    p2.add_argument("--tol", type=int, default=2, help="Tolerance (bp) when matching junction coords")
    p2.add_argument("--min-support", type=int, default=20, dest="min_support")
    p2.add_argument("--pseudocount", type=float, default=0.5)
    p2.add_argument("--out-csv", default="MAPT_3R4R_per_sample.csv")

    # batch
    p3 = sub.add_parser("batch", help="Run extraction for many BAMs, then compute Ψ/ratio once")
    src = p3.add_mutually_exclusive_group(required=True)
    src.add_argument("--bam-dir", help="Directory with *.bam")
    src.add_argument("--bam-list", help="Text file with BAM paths (one per line)")
    p3.add_argument("--gtf", required=True, help="MAPT-only GTF")
    p3.add_argument("--region", default="17:45890000-46050000")
    p3.add_argument("--outdir", default="out")
    p3.add_argument("--min-anchor", type=int, default=8)
    p3.add_argument("--min-intron", type=int, default=50)
    p3.add_argument("--max-intron", type=int, default=500000)
    p3.add_argument("--force-strand", choices=["+", "-", "."], default=".")
    p3.add_argument("--canonical", action="store_true")
    p3.add_argument("--tol", type=int, default=2)
    p3.add_argument("--min-support", type=int, default=20, dest="min_support")
    p3.add_argument("--pseudocount", type=float, default=0.5)

    return ap

def main():
    ap = build_parser()
    args = ap.parse_args()

    if args.cmd == "extract-juncs":
        extract_junctions(
            bam=args.bam,
            region=args.region,
            out=args.out,
            min_anchor=args.min_anchor,
            min_intron=args.min_intron,
            max_intron=args.max_intron,
            force_strand=args.force_strand,
        )

    elif args.cmd == "compute-e10":
        compute_e10(
            gtf=args.gtf,
            junc_dir=args.junc_dir,
            canonical=args.canonical,
            tol=args.tol,
            min_support=args.min_support,
            pseudocount=args.pseudocount,
            out_csv=args.out_csv,
        )

    elif args.cmd == "batch":
        batch_pipeline(
            bam_dir=getattr(args, "bam_dir", None),
            bam_list=getattr(args, "bam_list", None),
            gtf=args.gtf,
            region=args.region,
            outdir=args.outdir,
            min_anchor=args.min_anchor,
            min_intron=args.min_intron,
            max_intron=args.max_intron,
            force_strand=args.force_strand,
            canonical=args.canonical,
            tol=args.tol,
            min_support=args.min_support,
            pseudocount=args.pseudocount,
        )

if __name__ == "__main__":
    main()
