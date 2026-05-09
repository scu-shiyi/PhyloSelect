from string import Template
from pathlib import Path
import subprocess
import shutil
import re
from scipy.stats import chi2
import os
import tempfile
import pandas as pd
import ast

TEMPLATES = {
    "M0": Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 3
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    "M3": Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 3
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 3
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    'M7': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 7
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0
method = 0
cleandata = 0'''),
    'M8': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 8
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    'FREERATIO': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 1
NSsites = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
method = 0
fix_blength = 0''')
}

def _run_codeml_async(codeml_bin, ctl_path, threads, work_dir):
    work_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(max(1, int(threads)))
    return subprocess.Popen([codeml_bin, str(ctl_path)], cwd=str(work_dir), env=env)

def _run_codeml_pair(codeml_bin, jobs, total_threads):
    """
    jobs: list of (ctl_path, work_dir)
    total_threads: 用户给的总线程预算
    """
    jobs = list(jobs)
    if not jobs:
        return

    total_threads = max(1, int(total_threads))

    # 只有一个任务，直接用全部线程
    if len(jobs) == 1:
        ctl_path, work_dir = jobs[0]
        p = _run_codeml_async(codeml_bin, ctl_path, total_threads, work_dir)
        ret = p.wait()
        if ret != 0:
            raise RuntimeError(f"codeml failed with exit code {ret}")
        return

    # 多个任务时：
    # 线程数太小就串行，避免超预算
    if total_threads <= 1:
        for ctl_path, work_dir in jobs:
            p = _run_codeml_async(codeml_bin, ctl_path, 1, work_dir)
            ret = p.wait()
            if ret != 0:
                raise RuntimeError(f"codeml failed with exit code {ret}")
        return

    # 否则并发平分线程
    t1 = max(1, total_threads // len(jobs))
    remain = total_threads - t1 * len(jobs)
    per_job_threads = [t1] * len(jobs)
    for i in range(remain):
        per_job_threads[i] += 1

    procs = []
    for (ctl_path, work_dir), t in zip(jobs, per_job_threads):
        procs.append(_run_codeml_async(codeml_bin, ctl_path, t, work_dir))

    for p in procs:
        ret = p.wait()
        if ret != 0:
            raise RuntimeError(f"codeml failed with exit code {ret}")
        
def _mlc_complete(p: Path) -> bool:
    try:
        if not p.exists() or p.stat().st_size < 200:
            return False
        txt = p.read_text(errors="ignore")
        return ("lnL(" in txt) and ("Time used" in txt)
    except Exception:
        return False


def _entropy_csv_complete(path):
    try:
        path = Path(path)
        if not path.exists() or path.stat().st_size == 0:
            return False
        df = pd.read_csv(path, index_col=0, header=0)
        if df.empty or df.shape[0] == 0 or df.shape[1] == 0:
            return False
        # 抽样检查前几个单元格是否可解析
        values = df.values.flatten()
        checked = 0
        valid = 0

        for x in values:
            if pd.isna(x):
                continue
            checked += 1
            try:
                if isinstance(x, str):
                    x = ast.literal_eval(x)
                float(x) if not isinstance(x, (list, tuple)) else float(x[0])
                valid += 1
            except Exception:
                pass

            if checked >= 10:
                break

        return valid > 0

    except Exception:
        return False
    
def write_ctl(model, seq_path, tree_path, mlc_path, ctl_path):
    model = model.upper()
    if model not in TEMPLATES:
        raise ValueError(f"Unknown model: {model}")
    text = TEMPLATES[model].substitute(
        seq_path=str(seq_path),
        tree_path=str(tree_path),
        output_mlc=str(mlc_path)
    )
    ctl_path.write_text(text)
    return ctl_path, mlc_path


def parse_mlc_np_lnl(path: Path):

    NUM = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
    PAT = re.compile(
    rf'lnL\s*\(\s*ntime:\s*\d+\s+np:\s*(\d+)\s*\)\s*:\s*({NUM})')

    text = Path(path).read_text(errors="ignore")
    matches = list(PAT.finditer(text))
    if not matches:
        raise ValueError(f"无法在 {path} 里找到 lnL(ntime: .. np: ..): 行。")
    m = matches[-1]
    return {
        "np": int(m.group(1)),
        "lnL": float(m.group(2)),
    }

def parse_parameters(path, model_type):
    path = Path(path)
    text = path.read_text(errors="ignore")

    base = parse_mlc_np_lnl(path)
    result = {
        "model": model_type.lower(),
        "np": base["np"],
        "lnL": base["lnL"],
    }

    NUM = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
    model_type = model_type.lower()

    if model_type == 'm0':
        omega_pat = re.compile(r'omega\s*\(dN/dS\)\s*=\s*(' + NUM + r')')
        m = omega_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 omega")
        result["omega"] = float(m.group(1))

        kappa_pat = re.compile(r'kappa\s*\(ts/tv\)\s*=\s*(' + NUM + r')')
        m = kappa_pat.search(text)
        result["kappa"] = float(m.group(1)) if m else None
        return result
    
    elif model_type == 'm3':
        # K 值
        k_pat = re.compile(r'MLEs of dN/dS \(w\) for site classes \(K\s*=\s*(\d+)\)')
        m = k_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 M3 的 K 值")
        k = int(m.group(1))
        result["K"] = k
        # p 行
        p_pat = re.compile(r'^\s*p:\s+(.+)$', re.MULTILINE)
        m = p_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 M3 的 p 值")
        p_values = [float(x) for x in re.findall(NUM, m.group(1))]

        # w 行
        w_pat = re.compile(r'^\s*w:\s+(.+)$', re.MULTILINE)
        m = w_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 M3 的 w 值")
        w_values = [float(x) for x in re.findall(NUM, m.group(1))]

        # 动态写入 p0/p1/p2... 和 w0/w1/w2...
        result["p_values"] = p_values[:k]
        result["w_values"] = w_values[:k]

        site_pat = re.compile(
            r'^\s*(\d+)\s+([A-Za-z\-\?])\s+(' + NUM + r')(\*{1,2})\s+(' + NUM + r')',
            re.MULTILINE
        )
        sig_sites = []
        for m in site_pat.finditer(text):
            pos = m.group(1)
            aa = m.group(2)
            stars = m.group(4)
            sig_sites.append(f"{pos}{aa}{stars}")

        result["significant_sites"] = sig_sites
        result["significant_sites_str"] = " ".join(sig_sites) if sig_sites else "None"

        return result

    elif model_type == 'm7':
        pq_pat = re.compile(
        r'Parameters in M7 \(beta\):\s*p\s*=\s*(' + NUM + r')\s+q\s*=\s*(' + NUM + r')',
        re.MULTILINE)
        m = pq_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 M7 的 p 和 q")

        result["p"] = float(m.group(1))
        result["q"] = float(m.group(2))
        return result

    elif model_type == 'm8':
        m8_pat = re.compile(
        r'Parameters in M8 \(beta&w>1\):\s*'
        r'p0\s*=\s*(' + NUM + r')\s+'
        r'p\s*=\s*(' + NUM + r')\s+'
        r'q\s*=\s*(' + NUM + r')\s*'
        r'\(\s*p1\s*=\s*(' + NUM + r')\s*\)\s*'
        r'w\s*=\s*(' + NUM + r')',
        re.DOTALL)
        m = m8_pat.search(text)
        if not m:
            raise ValueError(f"无法在 {path} 中解析 M8 的参数")

        result["p0"] = float(m.group(1))
        result["p"] = float(m.group(2))
        result["q"] = float(m.group(3))
        result["p1"] = float(m.group(4))
        result["omega"] = float(m.group(5))

        site_pat = re.compile(r'^\s*(\d+)\s+([A-Za-z\-\?])\s+(' + NUM + r')(\*{1,2})\s+(' + NUM + r')',re.MULTILINE)
        # 1. 提取 NEB 区块
        neb_block_pat = re.compile(
            r'Naive Empirical Bayes \(NEB\) analysis.*?'
            r'Pr\(w>1\).*?\n(.*?)(?:\n\s*Bayes Empirical Bayes \(BEB\) analysis|\Z)',
            re.DOTALL
        )
        neb_sites = []
        neb_match = neb_block_pat.search(text)
        if neb_match:
            neb_block = neb_match.group(1)
            for m in site_pat.finditer(neb_block):
                pos = m.group(1)
                aa = m.group(2)
                stars = m.group(4)
                neb_sites.append(f"{pos}{aa}{stars}")

        # 2. 提取 BEB 区块
        beb_block_pat = re.compile(
        r'Bayes Empirical Bayes \(BEB\) analysis.*?'
        r'Pr\(w>1\).*?\n(.*)',
        re.DOTALL
    )
        beb_sites = []
        beb_match = beb_block_pat.search(text)
        if beb_match:
            beb_block = beb_match.group(1)
            for m in site_pat.finditer(beb_block):
                pos = m.group(1)
                aa = m.group(2)
                stars = m.group(4)
                beb_sites.append(f"{pos}{aa}{stars}")

        result["neb_sites"] = neb_sites
        result["neb_sites_str"] = " ".join(neb_sites) if neb_sites else "None"

        result["beb_sites"] = beb_sites
        result["beb_sites_str"] = " ".join(beb_sites) if beb_sites else "None"

        return result

def lrt(null_lnl, null_np, alt_lnl, alt_np):
    """
    Likelihood Ratio Test: H0 = simpler model (null), H1 = more complex model (alt)
    Inputs: lnL and np from .mlc files
    Returns: example.csv statistic, degrees of freedom, p-value, significance decision
    """
    df = alt_np - null_np
    if df <= 0:
        raise ValueError(f"Degrees of freedom (df)={df} <= 0, please check input np values.")
    stat = 2.0 * (alt_lnl - null_lnl)
    if stat < 0:
        stat = 0.0
    p = chi2.sf(stat, df)
    crit_005 = chi2.ppf(0.95, df)  
    crit_001 = chi2.ppf(0.99, df)  

    return {"stat": stat, 
            "df": df, 
            "p": p,
            "critical": {"0.05": crit_005, "0.01": crit_001},
            "sig": {"0.05": stat > crit_005, "0.01": stat > crit_001},}


def _make_workdirs(base_dir, gene_name, model, tmp_root = None):
    if tmp_root is None:
        tmp_root = base_dir / "tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)

    run_root = Path(tempfile.mkdtemp(prefix=f"{gene_name}_{model}_", dir=str(tmp_root)))
    workdirs = {
        "M0": run_root / "M0",
        "M3": run_root / "M3",
        "M7": run_root / "M7",
        "M8": run_root / "M8",
        "FREERATIO": run_root / "FREERATIO",
    }

    for d in workdirs.values():
        d.mkdir(parents=True, exist_ok=True)

    return run_root, workdirs


def _cleanup_tmp(run_root: Path, keep_intermediate: bool):
    if (not keep_intermediate) and run_root.exists():
        shutil.rmtree(run_root, ignore_errors=True)

def run_pair_model(seq_path,
                   tree_path,
                   out_dir, 
                   model, 
                   gene_name=None, 
                   tmp_root = None,
                   threads=1,
                   keep_intermediate=False):


    codeml_bin = shutil.which('codeml')
    if not codeml_bin:
        raise FileNotFoundError(
            "PAML not found. Please install PAML (codeml) and ensure it is on your PATH, "
            "or set the CODEML environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda paml\n"
            "  export CODEML=/path/to/codeml"
        )
    seq_path  = Path(seq_path)
    tree_path = Path(tree_path)

    if gene_name is None:
        gene_name = seq_path.stem

    base = Path(out_dir)
    res = base / "result"
    base.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)

    model_upper = model.upper()
    run_root, workdirs = _make_workdirs(base, gene_name, model_upper, tmp_root)

    try:
        if model_upper == "M0M3":
            null_ctl = base / "M0.ctl"
            alt_ctl = base / "M3.ctl"
            null_mlc = res / "M0.mlc"
            alt_mlc = res / "M3.mlc"
            # 写 ctl
            write_ctl("M0", seq_path, tree_path, null_mlc, null_ctl)
            write_ctl("M3", seq_path, tree_path, alt_mlc, alt_ctl)
            # run codeml
            jobs = []
            if not _mlc_complete(null_mlc):
                jobs.append((null_ctl, workdirs["M0"]))
            if not _mlc_complete(alt_mlc):
                jobs.append((alt_ctl, workdirs["M3"]))

            _run_codeml_pair(codeml_bin, jobs, threads)


            null_res = parse_mlc_np_lnl(null_mlc)
            null_lnl = null_res["lnL"]
            null_np = null_res["np"]
            alt_res = parse_mlc_np_lnl(alt_mlc)
            alt_lnl = alt_res["lnL"]
            alt_np = alt_res["np"]
            lrt_result = lrt(null_lnl, null_np, alt_lnl, alt_np)



            return {
                    "test": "M0_vs_M3",
                    "stat": lrt_result["stat"],
                    "df": lrt_result["df"],
                    "p": lrt_result["p"],
                    "critical": lrt_result["critical"],
                    "sig": lrt_result["sig"],
                    "result_dir": str(res),
                    "null_mlc": str(null_mlc),
                    "alt_mlc": str(alt_mlc),
                    "rst_path": None,
                    "summary_tsv": str(base / "summary.tsv"),
                }


        elif model_upper == "M7M8":
            null_ctl = base / "M7.ctl"
            alt_ctl = base / "M8.ctl"
            null_mlc = res / "M7.mlc"
            alt_mlc = res / "M8.mlc"
            rst_path = res / "M8.rst"

            write_ctl("M7", seq_path, tree_path, null_mlc, null_ctl)
            write_ctl("M8", seq_path, tree_path, alt_mlc, alt_ctl)


            jobs = []
            if not _mlc_complete(null_mlc):
                jobs.append((null_ctl, workdirs["M7"]))
            if not _mlc_complete(alt_mlc):
                jobs.append((alt_ctl, workdirs["M8"]))

            _run_codeml_pair(codeml_bin, jobs, threads)
            
            tmp_rst = workdirs["M8"] / "rst"
            if tmp_rst.exists():
                shutil.copy2(tmp_rst, rst_path)

            null_res = parse_mlc_np_lnl(null_mlc)
            null_lnl = null_res["lnL"]
            null_np = null_res["np"]
            alt_res = parse_mlc_np_lnl(alt_mlc)
            alt_lnl = alt_res["lnL"]
            alt_np = alt_res["np"]
            lrt_result = lrt(null_lnl, null_np, alt_lnl, alt_np)




            return {
                "test": "M7_vs_M8",
                "stat": lrt_result["stat"],
                "df": lrt_result["df"],
                "p": lrt_result["p"],
                "critical": lrt_result["critical"],
                "sig": lrt_result["sig"],
                "result_dir": str(res),
                "null_mlc": str(null_mlc),
                "alt_mlc": str(alt_mlc),
                "rst_path": str(rst_path) if rst_path.exists() else None,
            }
        
        elif model_upper == "FREERATIO":
            alt_ctl = base / "FREERATIO.ctl"
            alt_mlc = res / "FREERATIO.mlc"
            null_ctl = base / "M0.ctl"
            null_mlc = res / "M0.mlc"

            write_ctl("M0", seq_path, tree_path, null_mlc, null_ctl)
            write_ctl("FREERATIO", seq_path, tree_path, alt_mlc, alt_ctl)
            jobs = []
            if not _mlc_complete(null_mlc):
                jobs.append((null_ctl, workdirs["M0"]))

            if not _mlc_complete(alt_mlc):
                jobs.append((alt_ctl, workdirs["FREERATIO"]))

            _run_codeml_pair(codeml_bin, jobs, threads)


            null_res = parse_mlc_np_lnl(null_mlc)
            null_lnl = null_res["lnL"]
            null_np = null_res["np"]
            alt_res = parse_mlc_np_lnl(alt_mlc)
            alt_lnl = alt_res["lnL"]
            alt_np = alt_res["np"]
            lrt_result = lrt(null_lnl, null_np, alt_lnl, alt_np)



            return {
                "test": "M0_vs_FREERATIO",
                "stat": lrt_result["stat"],
                "df": lrt_result["df"],
                "p": lrt_result["p"],
                "critical": lrt_result["critical"],
                "sig": lrt_result["sig"],
                "result_dir": str(res),
                "null_mlc": str(null_mlc),
                "alt_mlc": str(alt_mlc),
                "rst_path": None,
                }
        else:
            raise ValueError(f"Unsupported model specification: {model}")
    finally:
        _cleanup_tmp(run_root, keep_intermediate)



if __name__ == "__main__":
    path = '/Users/sy/PhyloSelect/PhyloSelect-v2.1.1/DEMO/DEMO1/test/rbcL/paml_output/M7M8/result/M8.mlc'
    print(parse_parameters(path, "M8"))