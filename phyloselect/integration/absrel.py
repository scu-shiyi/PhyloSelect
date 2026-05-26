from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import math
import re
import pandas as pd
import shutil

class HyPhyError(RuntimeError):
    pass


@dataclass
class ABSRELBranchResult:
    branch: str
    omega_max: Optional[float] = None
    omega_weighted: Optional[float] = None
    baseline_omega:Optional[float] = None

    original_name: Optional[str] = None
    dn: Optional[float] = None
    ds: Optional[float] = None
    lrt: Optional[float] = None
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    notes: Optional[str] = None


class ABSRELRunner:
    """
    Minimal, non-interactive runner for `hyphy absrel`.

    Design goals:
    - Provide stable invocation for pipeline integration
    - Keep exposed parameters minimal (Nature-style CLI integration)
    - Parse JSON into branch-level metrics usable for environment association
    """

    def __init__(self):
        hyphy_path = shutil.which("hyphy")
        if not hyphy_path:
            raise FileNotFoundError(
                "HyPhy not found. Please install HyPhy and make sure it is on your PATH, "
                "or set the HYPHY environment variable to the absolute path of the executable.\n"
                "Examples:\n"
                "  conda install -c bioconda hyphy\n"
                "  export HYPHY=/path/to/hyphy"
            )
        self.hyphy_bin = hyphy_path

    def run(
        self,
        alignment: str | Path,
        tree: str | Path,
        output_json: str | Path,
        *,
        code: str = "Universal",
        branches: str = "All",
        multiple_hits: str = "None",
        srv: str = "No",
        extra_args: Optional[List[str]] = None,
        cwd: Optional[str | Path] = None,
        env: Optional[Dict[str, str]] = None,
        timeout_sec: Optional[int] = None,
    ) -> Path:
        """
        Execute HyPhy aBSREL and produce a JSON output file.
        """
        alignment = Path(alignment)
        tree = Path(tree)
        output_json = Path(output_json)

        if not alignment.exists():
            raise FileNotFoundError(f"Alignment not found: {alignment}")
        if not tree.exists():
            raise FileNotFoundError(f"Tree not found: {tree}")

        output_json.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.hyphy_bin,
            "absrel",
            "--alignment", str(alignment),
            "--tree", str(tree),
            "--output", str(output_json),
            "--code", code,
            "--branches", branches,
            "--multiple-hits", multiple_hits,
            "--srv", srv,
        ]

        if extra_args:
            cmd.extend(extra_args)

        try:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd is not None else None,
                env={**os.environ, **(env or {})},
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_sec,
                check=False,
            )
        except subprocess.TimeoutExpired as e:
            raise HyPhyError(f"HyPhy aBSREL timed out: {e}") from e

        if proc.returncode != 0:
            raise HyPhyError(
                "HyPhy aBSREL failed.\n"
                f"Command: {' '.join(cmd)}\n"
                f"STDOUT:\n{proc.stdout}\n"
                f"STDERR:\n{proc.stderr}\n"
            )

        if not output_json.exists():
            raise HyPhyError(
                "HyPhy completed but JSON output was not created.\n"
                f"Expected: {output_json}\n"
                f"STDOUT:\n{proc.stdout}\n"
                f"STDERR:\n{proc.stderr}\n"
            )

        return output_json

    def parse_branch_results(self, absrel_json: str | Path, *, tips_only: bool = True) -> List[ABSRELBranchResult]:
        absrel_json = Path(absrel_json)
        with absrel_json.open("r", encoding="utf-8") as f:
            data = json.load(f)

        branch_records = self._extract_branch_attribute_container(data)

        results: List[ABSRELBranchResult] = []
        for branch_name, payload in branch_records.items():
            if tips_only and not self._is_tip_name(str(branch_name)):
                continue

            omega_weighted, omega_max, baseline_omega, dn, ds, note = self._extract_branch_metrics(payload)
            lrt = self._extract_first_numeric(payload, ["LRT"])
            pval = self._extract_first_numeric(payload, ["Uncorrected P-value"])
            qval = self._extract_first_numeric(payload, ["Corrected P-value"])
            original_name = None
            if isinstance(payload, dict):
                original_name = payload.get("original name")
            results.append(
                ABSRELBranchResult(
                    branch=str(branch_name),
                    omega_weighted=omega_weighted,
                    omega_max=omega_max,
                    baseline_omega=baseline_omega,
                    dn=dn,
                    ds=ds,
                    lrt=lrt,
                    p_value=pval,
                    q_value=qval,
                    notes=note,
                    original_name=original_name,
                )
            )

        results.sort(key=lambda r: r.branch)
        return results


    def _is_tip_name(self, name):
        return re.match(r"^Node\d+$", name) is None

    def _extract_omega_metrics(self, payload: Any) -> Tuple[Optional[float], Optional[float], Optional[str]]:
        """
        Extract omega_max and omega_weighted from a branch payload.

        Supports HyPhy aBSREL JSON formats, including:
        - "Rate Distributions": [[omega, weight], ...]
        - other heuristic fallbacks (omegas/proportions, rate classes, etc.)
        """
        note_parts: List[str] = []

        if not isinstance(payload, dict):
            return None, None, "branch payload is not a dict"

        omega_candidates: List[float] = []
        weighted: Optional[float] = None

        # (A) Primary for your JSON: "Rate Distributions": [[omega, weight], ...]
        rd = payload.get("Rate Distributions")
        if isinstance(rd, list) and rd:
            omegas: List[float] = []
            weights: List[float] = []
            ok = True
            for item in rd:
                if (
                        isinstance(item, list)
                        and len(item) >= 2
                        and isinstance(item[0], (int, float))
                        and isinstance(item[1], (int, float))
                ):
                    om = float(item[0])
                    wt = float(item[1])
                    # 排除 NaN/inf
                    if not (math.isfinite(om) and math.isfinite(wt)):
                        continue
                    omegas.append(om)
                    weights.append(wt)
                else:
                    ok = False

            if omegas:
                omega_candidates.extend(omegas)
                s = sum(weights)
                if s > 0:
                    weighted = sum(o * w for o, w in zip(omegas, weights)) / s
                    note_parts.append("omega metrics from Rate Distributions")
            elif not ok:
                note_parts.append("Rate Distributions present but unparseable")

        # (B) Fallback 1) direct omega list
        if not omega_candidates:
            for k in ["omegas", "omega", "Omega", "dN/dS", "rates"]:
                v = payload.get(k)
                if isinstance(v, list) and v and all(isinstance(x, (int, float)) for x in v):
                    omega_candidates.extend([float(x) for x in v])

        # (C) Fallback 2) nested class dicts
        if weighted is None:
            for k in ["rate classes", "rate_classes", "classes", "mixture", "mixture distribution"]:
                v = payload.get(k)
                if isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
                    omegas = []
                    props = []
                    for cls in v:
                        om = self._extract_first_numeric(cls, ["omega", "Omega", "dN/dS"])
                        pr = self._extract_first_numeric(cls, ["proportion", "weight", "probability", "pi"])
                        if om is not None:
                            omegas.append(float(om))
                        if pr is not None:
                            props.append(float(pr))
                    if omegas:
                        omega_candidates.extend(omegas)
                    if omegas and props and len(omegas) == len(props):
                        s = sum(props)
                        if s > 0:
                            weighted = sum(o * p for o, p in zip(omegas, props)) / s
                            note_parts.append("omega_weighted from mixture proportions")

        # (D) Fallback 3) separate arrays
        if weighted is None:
            omegas = payload.get("omegas")
            props = payload.get("proportions") or payload.get("weights") or payload.get("probabilities")
            if (
                    isinstance(omegas, list) and isinstance(props, list)
                    and len(omegas) == len(props)
                    and all(isinstance(x, (int, float)) for x in omegas)
                    and all(isinstance(x, (int, float)) for x in props)
            ):
                s = float(sum(props))
                if s > 0:
                    weighted = sum(float(o) * float(p) for o, p in zip(omegas, props)) / s
                    omega_candidates.extend([float(x) for x in omegas])
                    note_parts.append("omega_weighted from omegas+proportions")

        omega_max = max(omega_candidates) if omega_candidates else None
        if omega_max is None and weighted is None:
            note_parts.append("no omega fields detected in payload")

        return omega_max, weighted, "; ".join(note_parts) if note_parts else None

    def _extract_first_numeric(self, d: Any, keys: List[str]) -> Optional[float]:
        if not isinstance(d, dict):
            return None
        for k in keys:
            v = d.get(k)
            if isinstance(v, (int, float)):
                return float(v)
        return None

    def _extract_branch_attribute_container(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Return a dict mapping branch_name -> attributes payload.
        Handles HyPhy aBSREL layouts like:
          "branch attributes": {"0": {branch->payload}, "attributes": {...}}
        """
        ba = data.get("branch attributes")

        if isinstance(ba, dict):
            # Most common in your file: {"0": {...branch map...}, "attributes": {...meta...}}
            if "0" in ba and isinstance(ba["0"], dict):
                # Ensure it's actually a branch map (values are dict payloads)
                if ba["0"] and all(isinstance(v, dict) for v in ba["0"].values()):
                    return ba["0"]

            # More general: sometimes numeric keys other than "0"
            for k, v in ba.items():
                if isinstance(k, str) and k.isdigit() and isinstance(v, dict) and v and all(
                        isinstance(vv, dict) for vv in v.values()):
                    return v

            # If it's already branch->payload (no wrapper layer)
            if ba and all(isinstance(v, dict) for v in ba.values()):
                # Guard: avoid returning wrapper dicts like {"attributes": {...}}
                # Require at least one branch-like key
                keys = [str(x) for x in ba.keys()]
                if any(re.match(r"^Node\d+$", kk) or re.search(r"[A-Za-z]", kk) for kk in keys):
                    return ba

        # Alternate: "branches"
        br = data.get("branches")
        if isinstance(br, dict) and br and all(isinstance(v, dict) for v in br.values()):
            return br

        raise HyPhyError(
            "Unrecognized aBSREL JSON structure: cannot locate branch attribute container."
        )

    import numpy as np

    def to_dataframe(
            self,
            absrel_json: Union[str, Path],
            *,
            tips_only: bool = True,
            ds_eps: float = 1e-6,
            dn_eps: float = 1e-6,
            omega_floor: float = 1e-6,
            omega_weighted_cap: float = 100,
            omega_max_cap: float = 1000,
            baseline_omega_cap: float = 100,
            omega_dn_ds_cap: float = 100,
            unstable_ds: float = 1e-4,
            unstable_omega: float = 10,
            dn_ds_ratio_cap: float = 10,
    ) -> pd.DataFrame:
        results = self.parse_branch_results(absrel_json, tips_only=tips_only)
        rows = []

        for r in results:
            dn = r.dn
            ds = r.ds

            omega_dn_ds = None
            if dn is not None and ds is not None and math.isfinite(dn) and math.isfinite(ds) and ds >= ds_eps:
                omega_dn_ds = dn / ds

            flag_ds_small = (ds is None) or (not math.isfinite(ds)) or (ds < ds_eps)
            flag_dn_small = (dn is None) or (not math.isfinite(dn)) or (dn < dn_eps)

            flag_low_signal = flag_dn_small and flag_ds_small

            flag_extreme_omega_weighted = (
                    r.omega_weighted is not None
                    and math.isfinite(r.omega_weighted)
                    and r.omega_weighted >= omega_weighted_cap
            )

            flag_extreme_omega_max = (
                    r.omega_max is not None
                    and math.isfinite(r.omega_max)
                    and r.omega_max >= omega_max_cap
            )

            flag_extreme_baseline_omega = (
                    r.baseline_omega is not None
                    and math.isfinite(r.baseline_omega)
                    and r.baseline_omega >= baseline_omega_cap
            )

            flag_extreme_omega_dn_ds = (
                    omega_dn_ds is not None
                    and math.isfinite(omega_dn_ds)
                    and omega_dn_ds >= omega_dn_ds_cap
            )

            flag_dn_ds_ratio_high = (
                    omega_dn_ds is not None
                    and math.isfinite(omega_dn_ds)
                    and omega_dn_ds > dn_ds_ratio_cap
            )

            flag_omega_cap = (
                    flag_extreme_omega_weighted
                    or flag_extreme_omega_max
                    or flag_extreme_baseline_omega
                    or flag_extreme_omega_dn_ds
            )

            flag_unstable_omega = (
                    ds is not None
                    and math.isfinite(ds)
                    and ds < unstable_ds
                    and (
                            (r.omega_weighted is not None and math.isfinite(
                                r.omega_weighted) and r.omega_weighted > unstable_omega)
                            or (omega_dn_ds is not None and math.isfinite(omega_dn_ds) and omega_dn_ds > unstable_omega)
                    )
            )

            keep = (
                    not flag_ds_small
                    and not flag_low_signal
                    and not flag_omega_cap
                    and not flag_unstable_omega
                    and not flag_dn_ds_ratio_high
            )

            exclusion_reasons = []
            if flag_ds_small:
                exclusion_reasons.append("small ds")
            if flag_low_signal:
                exclusion_reasons.append("low substitution signal")
            if flag_extreme_omega_weighted:
                exclusion_reasons.append("extreme omega_weighted")
            if flag_extreme_omega_max:
                exclusion_reasons.append("extreme omega_max")
            if flag_extreme_baseline_omega:
                exclusion_reasons.append("extreme baseline_omega")
            if flag_extreme_omega_dn_ds:
                exclusion_reasons.append("extreme omega_dn_ds")
            if flag_unstable_omega:
                exclusion_reasons.append("small ds; unstable omega")
            if flag_dn_ds_ratio_high:
                exclusion_reasons.append("dn/ds ratio > 10")

            rows.append({
                "taxon": r.original_name or r.branch,
                "hyphy_name": r.branch,
                "original_name": r.original_name,

                "omega_weighted": r.omega_weighted,
                "omega_max": r.omega_max,
                "baseline_omega": r.baseline_omega,

                "dn": dn,
                "ds": ds,
                "omega_dn_ds": omega_dn_ds,

                "lrt": r.lrt,
                "p_value": r.p_value,
                "q_value": r.q_value,
                "significant": (r.q_value is not None and r.q_value <= 0.05),

                "flag_ds_small": flag_ds_small,
                "flag_dn_small": flag_dn_small,
                "flag_low_signal": flag_low_signal,
                "flag_omega_cap": flag_omega_cap,
                "flag_extreme_omega_weighted": flag_extreme_omega_weighted,
                "flag_extreme_omega_max": flag_extreme_omega_max,
                "flag_extreme_baseline_omega": flag_extreme_baseline_omega,
                "flag_extreme_omega_dn_ds": flag_extreme_omega_dn_ds,
                "flag_unstable_omega": flag_unstable_omega,
                "flag_dn_ds_ratio_high": flag_dn_ds_ratio_high,

                "keep": keep,
                "exclusion_reason": "; ".join(exclusion_reasons),

                "notes": r.notes,
            })

        df = pd.DataFrame(rows)

        if not df.empty:
            df = df.sort_values("taxon").reset_index(drop=True)

        numeric_cols = [
            "omega_weighted", "omega_max", "baseline_omega",
            "dn", "ds", "omega_dn_ds",
            "lrt", "p_value", "q_value",
        ]

        for c in numeric_cols:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        return df

    def _extract_branch_metrics(self, payload: Any) -> Tuple[
        Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[str]
    ]:
        """
        Return (omega_weighted, omega_max, baseline_omega, dn, ds, note)
        """
        if not isinstance(payload, dict):
            return None, None, None, None, None, "branch payload is not a dict"

        note_parts: List[str] = []

        # 1) baseline omega
        baseline_omega = self._extract_first_numeric(payload, ["Baseline MG94xREV omega ratio"])

        # 2) dN / dS (adaptive model)
        dn = self._extract_first_numeric(payload, ["Full adaptive model (non-synonymous subs/site)"])
        ds = self._extract_first_numeric(payload, ["Full adaptive model (synonymous subs/site)"])

        # 3) omega from Rate Distributions: [[omega, weight], ...]
        omega_weighted: Optional[float] = None
        omega_max: Optional[float] = None

        lrt = self._extract_first_numeric(payload, ["LRT"])
        pval = self._extract_first_numeric(payload, ["Uncorrected P-value"])
        qval = self._extract_first_numeric(payload, ["Corrected P-value"])

        rd = payload.get("Rate Distributions")
        if isinstance(rd, list) and rd:
            omegas: List[float] = []
            weights: List[float] = []
            for item in rd:
                if (
                        isinstance(item, list)
                        and len(item) >= 2
                        and isinstance(item[0], (int, float))
                        and isinstance(item[1], (int, float))
                ):
                    om = float(item[0])
                    wt = float(item[1])
                    if not (math.isfinite(om) and math.isfinite(wt)):
                        continue
                    omegas.append(om)
                    weights.append(wt)

            if omegas:
                omega_max = max(omegas)
                s = sum(weights)
                if s > 0:
                    omega_weighted = sum(o * w for o, w in zip(omegas, weights)) / s
                note_parts.append("omega from Rate Distributions")
            else:
                note_parts.append("Rate Distributions present but empty/unparseable")

        # 4) fallback: if Rate Distributions missing, use baseline_omega as last resort
        if omega_weighted is None and baseline_omega is not None:
            omega_weighted = baseline_omega
            omega_max = baseline_omega
            note_parts.append("omega fallback to baseline_omega")

        return omega_weighted, omega_max, baseline_omega, dn, ds, "; ".join(note_parts) if note_parts else None


if __name__ == "__main__":
    import shutil
    hyphy_bin = shutil.which("hyphy")
    out_dir = '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/scripts/integration'
    fasta_path = '/Users/sy/PhyloSelect/准备删除/testing_data/模块1/demo1/rbcL/file_input/rbcL.phy'
    tree_path = '/Users/sy/PhyloSelect/准备删除/testing_data/模块1/demo1/rbcL/file_input/rbcL.paml.tree'
    out_name = Path(fasta_path).stem
    json_path = Path(out_dir) / f"{out_name}.ABSREL.json"

    runner = ABSRELRunner(hyphy_bin)
    # json_path = runner.run(fasta_path, tree_path,json_path)

    json_path = "/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/scripts/integration/rbcL.ABSREL.json"
    # runner = ABSRELRunner(hyphy_bin=hyphy_bin)
    df = runner.to_dataframe(json_path, tips_only=True)
    pd.set_option("display.float_format", lambda x: format(x, "g"))
    print(df)
    df.to_csv("/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/phyloselect/integration/rbcL.csv", index=False)
