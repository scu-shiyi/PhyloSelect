
# integration.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Sequence
from Bio import Phylo
from io import StringIO
from statsmodels.stats.knockoff_regeffects import CorrelationEffects
from phyloselect.integration.absrel import ABSRELRunner, ABSRELBranchResult
import numpy as np
import statsmodels.api as sm
from ete3 import Tree
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO


PathLike = Union[str, Path]

def _read_newick_text(tree_newick: PathLike) -> str:
    """
    Accept either:
      - a path to a tree file
      - a raw newick string

    Also handles tree files that start with a PHYLIP-style header line like:
      "18  1"
    by skipping lines until a line containing '(' is found.
    """
    p = Path(str(tree_newick))
    if p.exists():
        text = p.read_text(encoding="utf-8", errors="ignore").strip()
    else:
        text = str(tree_newick).strip()

    # If it's already a single-line newick, return
    if "(" in text and ";" in text and text.lstrip().startswith("("):
        return text

    # Otherwise, try to find the first line that looks like newick
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    for ln in lines:
        if "(" in ln and ";" in ln:
            return ln
    # fallback: join and hope
    joined = "".join(lines)
    if "(" in joined and ";" in joined:
        return joined
    raise ValueError("无法从 tree_newick 解析 Newick。请提供标准 newick（含括号和分号）。")


def _build_vcv_from_tree(newick_text: str) -> pd.DataFrame:
    """
    Brownian-motion VCV:
      V[i,i] = distance(root, tip_i)
      V[i,j] = distance(root, MRCA(tip_i, tip_j))
    """
    tree = Phylo.read(StringIO(newick_text), "newick")
    tips = [t.name for t in tree.get_terminals()]
    if any(t is None or str(t).strip() == "" for t in tips):
        raise ValueError("树上存在没有名字的 tip（terminal）。请检查 newick 的 tip 名称。")

    n = len(tips)
    vcv = np.zeros((n, n), dtype=float)

    # Precompute root->tip distances
    root = tree.root
    dist_root_tip = {t: tree.distance(root, t) for t in tips}

    # Fill VCV
    for i, ti in enumerate(tips):
        vcv[i, i] = dist_root_tip[ti]
        for j in range(i + 1, n):
            tj = tips[j]
            mrca = tree.common_ancestor(ti, tj)
            cov = tree.distance(root, mrca)
            vcv[i, j] = cov
            vcv[j, i] = cov

    return pd.DataFrame(vcv, index=tips, columns=tips)

def run_absrel_and_parse(
    alignment: PathLike,
    tree: PathLike,
    *,
    out_dir: PathLike = "results/absrel",
    code: str = "Universal",
    branches: str = "All",
    multiple_hits: str = "None",
    srv: str = "No",
    timeout_sec: Optional[int] = None,
) -> Dict[str, Any]:

    alignment = Path(alignment)
    tree = Path(tree)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)


    out_name = alignment.stem

    json_path = out_dir / f"{out_name}.ABSREL.json"
    runner = ABSRELRunner()
    # 1) run
    json_path = runner.run(
        alignment=alignment,
        tree=tree,
        output_json=json_path,
        code=code,
        branches=branches,
        multiple_hits=multiple_hits,
        srv=srv,
        timeout_sec=timeout_sec,
    )

    # 2) parse
    df = runner.to_dataframe(json_path, tips_only=True)
    rename_map = {}
    if "LRT" in df.columns and "lrt" not in df.columns:
        rename_map["LRT"] = "lrt"
    if "Uncorrected P-value" in df.columns and "p_value_raw" not in df.columns:
        rename_map["Uncorrected P-value"] = "p_value_raw"
    if "Corrected P-value" in df.columns and "q_value_fdr" not in df.columns:
        rename_map["Corrected P-value"] = "q_value_fdr"
    if "p_value" in df.columns and "p_value_raw" not in df.columns:
        rename_map["p_value"] = "p_value_raw"
    if "q_value" in df.columns and "q_value_fdr" not in df.columns:
        rename_map["q_value"] = "q_value_fdr"


    if rename_map:
        df = df.rename(columns=rename_map)

    if "p_value_raw" not in df.columns:
        df["p_value_raw"] = np.nan
    if "q_value_fdr" not in df.columns:
        df["q_value_fdr"] = np.nan

    if df["q_value_fdr"].notna().any():
        df["selected_branch"] = df["q_value_fdr"] <= 0.05
    else:
        df["selected_branch"] = df["p_value_raw"].notna() & (df["p_value_raw"] <= 0.05)

    return df

def _fmt_effect(x):
    """效应值/omega/dN/dS：保留 4 位有效数字"""
    if pd.isna(x):
        return ""
    return f"{x:.4g}"

def _fmt_lrt(x):
    """LRT：保留 3 位小数"""
    if pd.isna(x):
        return ""
    return f"{x:.3f}"

def _fmt_t(x):
    """t 值：保留 3 位小数"""
    if pd.isna(x):
        return ""
    return f"{x:.3f}"

def _fmt_pq(x):
    """P/Q 值：小于 0.001 显示 <0.001，否则保留 3 位小数"""
    if pd.isna(x):
        return ""
    if x < 0.001:
        return "<0.001"
    return f"{x:.3f}"

def _make_exclusion_reason(row):
    reasons = []
    if row.get("flag_ds_small", False):
        reasons.append("small dS")
    if row.get("flag_omega_cap", False):
        reasons.append("extreme ω")
    return "; ".join(reasons)



def validate_alignment_env_taxa(alignment, env_csv, taxon_col="taxon"):
    aln_taxa = {rec.id.strip() for rec in SeqIO.parse(str(alignment), "fasta") if rec.id}
    env = pd.read_csv(env_csv)
    env_taxa = set(env[taxon_col].astype(str).str.strip())
    aln_taxa.discard("")
    env_taxa.discard("")
    if aln_taxa != env_taxa:
        raise ValueError("Mismatch between alignment taxa and environment taxa")
    
def run_pgls_reporting(
        absrel_df,
        env_csv,
        tree_newick,
        outdir,
        *,
        log_transform: bool = True,
        alpha: float = 0.05,
        min_n:int = 5,
        ):
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df = absrel_df.copy()

    rename = {}
    if "LRT" in df.columns and "lrt" not in df.columns:
        rename["LRT"] = "lrt"
    if "Uncorrected P-value" in df.columns and "p_value_raw" not in df.columns:
        rename["Uncorrected P-value"] = "p_value_raw"
    if "Corrected P-value" in df.columns and "q_value_fdr" not in df.columns:
        rename["Corrected P-value"] = "q_value_fdr"
    if "p_value" in df.columns and "p_value_raw" not in df.columns:
        rename["p_value"] = "p_value_raw"
    if "q_value" in df.columns and "q_value_fdr" not in df.columns:
        rename["q_value"] = "q_value_fdr"
    df = df.rename(columns=rename)

    taxon_col = "taxon"
    omega_col = "omega_weighted"

    for c in [taxon_col, omega_col]:
        if c not in df.columns:
            raise ValueError(f"absrel_df 缺少列: {c}")

    for c in ["lrt", "p_value_raw", "q_value_fdr"]:
        if c not in df.columns:
            df[c] = np.nan
    keep_col = 'keep'
    if keep_col not in df.columns:
        df[keep_col] = True
    if "selected_branch" not in df.columns:
        if df["q_value_fdr"].notna().any():
            df["selected_branch"] = df["q_value_fdr"] <= 0.05
        else:
            df["selected_branch"] = df["p_value_raw"].notna() & (df["p_value_raw"] <= 0.05)

    # -------- dtype cleaning --------
    df[taxon_col] = df[taxon_col].astype(str)
    df[omega_col] = pd.to_numeric(df[omega_col], errors="coerce")
    df["lrt"] = pd.to_numeric(df["lrt"], errors="coerce")
    df["p_value_raw"] = pd.to_numeric(df["p_value_raw"], errors="coerce")
    df["q_value_fdr"] = pd.to_numeric(df["q_value_fdr"], errors="coerce")
    df[keep_col] = df[keep_col].astype(bool)

    # ==================================================
    # Table S1: full aBSREL branch results
    # ==================================================
    table_s1_cols = [
        "taxon",
        "omega_weighted",
        "omega_max",
        "baseline_omega",
        "dn",
        "ds",
        "omega_dn_ds",
        "lrt",
        "p_value_raw",
        "q_value_fdr",
        "selected_branch",
        keep_col,
        "flag_ds_small",
        "flag_omega_cap",
    ]
    table_s1_cols = [c for c in table_s1_cols if c in df.columns]
    table_s1 = df[table_s1_cols].copy()

    # 新增更直观的排除原因列
    table_s1["exclusion_reason"] = table_s1.apply(_make_exclusion_reason, axis=1)

    table_s1 = table_s1.sort_values(
        by=["selected_branch", "q_value_fdr", "p_value_raw", "taxon"],
        ascending=[False, True, True, True],
        na_position="last"
    ).reset_index(drop=True)

    table_s1_export = table_s1.rename(columns={
        "taxon": "Taxon",
        "omega_weighted": "ω_weighted",
        "omega_max": "ω_max",
        "baseline_omega": "Baseline ω",
        "dn": "dN",
        "ds": "dS",
        "omega_dn_ds": "dN/dS",
        "lrt": "LRT",
        "p_value_raw": "P-value",
        "q_value_fdr": "Q-value",
        "selected_branch": "Significant",
        keep_col: "Retained for analysis",
        "exclusion_reason": "Exclusion reason",
    })

    # 不再导出原始 flag 列
    drop_cols = [c for c in ["flag_ds_small", "flag_omega_cap"] if c in table_s1_export.columns]
    if drop_cols:
        table_s1_export = table_s1_export.drop(columns=drop_cols)

    for col in ["ω_weighted", "ω_max", "Baseline ω", "dN", "dS", "dN/dS"]:
        if col in table_s1_export.columns:
            table_s1_export[col] = table_s1_export[col].map(_fmt_effect)

    if "LRT" in table_s1_export.columns:
        table_s1_export["LRT"] = table_s1_export["LRT"].map(_fmt_lrt)
    if "P-value" in table_s1_export.columns:
        table_s1_export["P-value"] = table_s1_export["P-value"].map(_fmt_pq)
    if "Q-value" in table_s1_export.columns:
        table_s1_export["Q-value"] = table_s1_export["Q-value"].map(_fmt_pq)

    f_s1 = output_dir / "TableS1_aBSREL_full_branch_results.csv"
    table_s1_export.to_csv(f_s1, index=False, encoding="utf-8")

    # =========================================================
    # Table 1: significant aBSREL branches
    # =========================================================
    table_1 = df[[taxon_col, "lrt", "p_value_raw", "q_value_fdr", "selected_branch"]].copy()

    if table_1["q_value_fdr"].notna().any():
        table_1["selected_branch"] = table_1["q_value_fdr"] <= 0.05
    else:
        table_1["selected_branch"] = table_1["p_value_raw"].notna() & (table_1["p_value_raw"] <= 0.05)

    table_1 = table_1.sort_values(
        by=["selected_branch", "q_value_fdr", "p_value_raw", taxon_col],
        ascending=[False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)

    table_1_export = table_1.rename(columns={
        "taxon": "Taxon",
        "lrt": "LRT",
        "p_value_raw": "P-value",
        "q_value_fdr": "Q-value",
        "selected_branch": "Significant",
    })

    table_1_export["LRT"] = table_1_export["LRT"].map(_fmt_lrt)
    table_1_export["P-value"] = table_1_export["P-value"].map(_fmt_pq)
    table_1_export["Q-value"] = table_1_export["Q-value"].map(_fmt_pq)

    f_t1 = output_dir / "Table1_aBSREL_significant_branches.csv"
    table_1_export.to_csv(f_t1, index=False)

    # =========================================================
    # Internal merged table for PGLS (not exported)
    # =========================================================
    env = pd.read_csv(env_csv)
    if taxon_col not in env.columns:
        raise ValueError(f"Missing column in environment table: {taxon_col}")
    env[taxon_col] = env[taxon_col].astype(str)

    env_cols = [c for c in env.columns if c != taxon_col]

    print("df taxa:", sorted(df[taxon_col].astype(str).unique())[:30])
    print("env taxa:", sorted(env[taxon_col].astype(str).unique())[:30])
    print("intersection:", set(df[taxon_col].astype(str)) & set(env[taxon_col].astype(str)))

    t2 = df[[taxon_col, omega_col, keep_col]].copy()
    t2 = t2.merge(env[[taxon_col] + env_cols], on=taxon_col, how="inner")
    analysis_df = t2[t2[keep_col]].copy()
    y = analysis_df[omega_col].copy()
    if log_transform:
        y = np.log1p(y)

    newick_text = _read_newick_text(tree_newick)
    vcv = _build_vcv_from_tree(newick_text)

    results = []
    for col in env_cols:
        x = pd.to_numeric(analysis_df[col], errors="coerce")
        valid = ~(y.isna() | x.isna())
        sub = analysis_df.loc[valid, [taxon_col]].copy()

        if len(sub) < min_n:
            continue

        sub["y"] = y[valid]
        sub["x"] = x[valid]

        taxa = [t for t in sub[taxon_col].tolist() if t in vcv.index]
        sub = sub[sub[taxon_col].isin(taxa)]

        if len(sub) < min_n:
            continue

        V = vcv.loc[sub[taxon_col], sub[taxon_col]].to_numpy()
        X = sm.add_constant(sub["x"].to_numpy())
        Y = sub["y"].to_numpy()

        try:
            model = sm.GLS(Y, X, sigma=V)
            res = model.fit()

            results.append({
                "env_factor": col,
                "n_taxa": len(sub),
                "beta": res.params[1],
                "std_error": res.bse[1],
                "t_value": res.tvalues[1],
                "p_value_raw": res.pvalues[1],
            })
        except Exception:
            continue

    pgls_df = pd.DataFrame(results)
    print(pgls_df)
    if not pgls_df.empty:
        reject, qvals, _, _ = multipletests(
            pgls_df["p_value_raw"], alpha=alpha, method="fdr_bh"
        )
        pgls_df["q_value_fdr"] = qvals
        pgls_df["significant_fdr"] = reject

        pgls_df = pgls_df.sort_values(
            by=["q_value_fdr", "p_value_raw", "env_factor"],
            ascending=[True, True, True],
            na_position="last",
        ).reset_index(drop=True)
    else:
        pgls_df = pd.DataFrame(columns=[
            "env_factor", "n_taxa", "beta", "std_error", "t_value",
            "p_value_raw", "q_value_fdr", "significant_fdr",
            "response_variable", "transformation", "model"
        ])

    pgls_export = pgls_df.rename(columns={
        "env_factor": "Environmental factor",
        "n_taxa": "N",
        "beta": "β",
        "std_error": "SE",
        "t_value": "t",
        "p_value_raw": "P-value",
        "q_value_fdr": "Q-value",
    })

    pgls_export["β"] = pgls_export["β"].map(_fmt_effect)
    pgls_export["SE"] = pgls_export["SE"].map(_fmt_effect)
    pgls_export["t"] = pgls_export["t"].map(_fmt_t)
    pgls_export["P-value"] = pgls_export["P-value"].map(_fmt_pq)
    pgls_export["Q-value"] = pgls_export["Q-value"].map(_fmt_pq)

    f_pgls = output_dir / "Table2_PGLS_environment_association.csv"
    pgls_export.to_csv(f_pgls, index=False)

    return {
        "TableS1": f_s1,
        "Table1": f_t1,
        "Table2": f_pgls,
    }


def run_absrel_env_pipeline(
        alignment: PathLike,
        tree: PathLike,
        env_csv: PathLike,
        outdir: PathLike,
        *,
        code: str = "Universal",
        branches: str = "All",
        multiple_hits: str = "None",
        srv: str = "No",
        timeout_sec: Optional[int] = None,
        log_transform: bool = True,
        alpha: float = 0.05,
        min_n: int = 5,
):
    
        # 先检查序列标签与环境表 taxon 是否一致
    validate_alignment_env_taxa(
        alignment=alignment,
        env_csv=env_csv,
        taxon_col="taxon",
    )

    df = run_absrel_and_parse(
        alignment=alignment,
        tree=tree,
        out_dir=outdir,
        code=code,
        branches=branches,
        multiple_hits=multiple_hits,
        srv=srv,
        timeout_sec=timeout_sec,
    )

    paths = run_pgls_reporting(
        absrel_df=df,
        env_csv=env_csv,
        tree_newick=tree,
        outdir=outdir,
        log_transform=log_transform,
        alpha=alpha,
        min_n=min_n,
    )
    print(paths)
    return {
        "absrel_df": df,
        "output_tables": paths,
    }



if __name__ == "__main__":
    out_dir = '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/scripts/integration'
    fasta_path = '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/data_creating/hhnew/4HPAAS3.fasta'
    tree_path = '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/data_creating/hhnewout/4HPAAS3/file_input/4HPAAS3.paml.tree'
    df = run_absrel_and_parse(alignment=fasta_path, tree=tree_path, out_dir=out_dir,)
    x = run_pgls_reporting(df,
                '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/phyloselect/integration/envs1.csv',
                tree_path,
                "/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/phyloselect/integration/out2"
                )
    print(x)
