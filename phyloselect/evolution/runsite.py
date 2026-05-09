from pathlib import Path
from phyloselect.utils.Phylip_Prepare import prepare_paml_input
from phyloselect.utils.EvoScoring import position_entropy
import pandas as pd
import ast
from phyloselect.utils.site_model import run_pair_model, _entropy_csv_complete, parse_parameters
from phyloselect.evolution.plot import Visualization

def Run_Site(fasta_path, 
             output_dir, 
             *,
             tree_path=None,
             outgroups=None, 
             label_max_chars=20,
             conservation_threshold=1.4, 
             entropy_csv=None,
             run_paml=True,
             threads=1,
             force=False,
             keep_intermediate=False,
             infer_tree=False,
             codon_aligned=False,
             skip_plot=False):


    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    
    if outgroups is None:
        outgroups = []

    gene_name = fasta_path.stem
    gene_dir = output_dir / gene_name
    file_input = gene_dir / "file_input"
    evo_output = gene_dir / "evo_output"
    paml_output = gene_dir / "paml_output"
    logs_dir = gene_dir / "logs"
    tmp_dir = gene_dir / "tmp"

    if gene_dir.exists() and force:
        import shutil
        shutil.rmtree(gene_dir)

        
    for d in [file_input, evo_output, paml_output,logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    if  keep_intermediate:
        tmp_dir.mkdir(parents=True, exist_ok=True)


    phylip_file, paml_treefile, _ = prepare_paml_input(fasta_path, file_input, tree_path=tree_path, infer_tree=infer_tree, codon_aligned=codon_aligned, keep_intermediate=keep_intermediate)

    default_entropy_csv =  evo_output / f"{gene_name}_entropy.csv"
    if entropy_csv is not None:
        entropy_csv = Path(entropy_csv)
    elif _entropy_csv_complete(default_entropy_csv):
        entropy_csv = default_entropy_csv
    else:
        print("[Status] Running Evo2 entropy calculation...")
        entropy_csv = Path(position_entropy(phylip_file, str(evo_output)))


    heatmap_df = pd.read_csv(entropy_csv, index_col=0, header=0)
    heatmap_df = heatmap_df.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    rst_path = None
    paml_summary_rows = []

    if run_paml:
        m0m3_base_dir = paml_output / "M0M3"
        m0m3_result = run_pair_model(phylip_file, 
                                     paml_treefile, 
                                     m0m3_base_dir, 
                                     "M0M3",
                                     gene_name=gene_name,
                                     tmp_root=tmp_dir if keep_intermediate else None,
                                     threads=threads,
                                     keep_intermediate=keep_intermediate)

        m0_parsed = parse_parameters(m0m3_result["null_mlc"], "m0")
        m3_parsed = parse_parameters(m0m3_result["alt_mlc"], "m3")
        paml_summary_rows.append({
                "Model": "M0: One ratio",
                "np": m0_parsed["np"],
                "lnL": round(m0_parsed["lnL"], 2),
                "Parameters": f"ω = {m0_parsed['omega']:.3f}",
                "Positively selected sites": "None"
            })
        
        m3_params = []
        for i, (p_val, w_val) in enumerate(zip(m3_parsed["p_values"], m3_parsed["w_values"])):
            m3_params.append(f"p{i} = {p_val:.3f}, ω{i} = {w_val:.3f}")
        m3_param_str = "; ".join(m3_params)

        paml_summary_rows.append({
            "Model": "M3: Discrete",
            "np": m3_parsed["np"],
            "lnL": round(m3_parsed["lnL"], 2),
            "Parameters": m3_param_str,
            "Positively selected sites": ", ".join(m3_parsed["significant_sites"]) if m3_parsed["significant_sites"] else "None",
        })
        if m0m3_result["p"] < 0.05:
            print("LRT (M0 vs M3) is significant at the 0.05/0.01 level: evidence of site-specific ω heterogeneity.")

            m7m8_base_dir = paml_output / "M7M8"
            m7m8_result = run_pair_model(phylip_file,
                                          paml_treefile, 
                                          m7m8_base_dir, 
                                          "M7M8", 
                                          gene_name=gene_name,
                                          tmp_root=tmp_dir if keep_intermediate else None,
                                          threads=threads,
                                          keep_intermediate=keep_intermediate,)


            rst_path = m7m8_result["rst_path"]

            m7_parsed = parse_parameters(m7m8_result["null_mlc"], "m7")
            m8_parsed = parse_parameters(m7m8_result["alt_mlc"], "m8")

            # ---------- M7 ----------
            paml_summary_rows.append({
                "Model": "M7: β",
                "np": m7_parsed["np"],
                "lnL": round(m7_parsed["lnL"], 2),
                "Parameters": f"p = {m7_parsed['p']:.3f}, q = {m7_parsed['q']:.3f}",
                "Positively selected sites": "Not allowed"
            })

            # ---------- M8 ----------
            # 参数字符串
            m8_param_str = (
                f"p0 = {m8_parsed['p0']:.3f}; "
                f"p = {m8_parsed['p']:.3f}, q = {m8_parsed['q']:.3f}; "
                f"p1 = {m8_parsed['p1']:.3f}, ω = {m8_parsed['omega']:.3f}"
            )

            # DEMO2（优先用 BEB）
            sites = m8_parsed.get("beb_sites") or []
            sites_str = ", ".join(sites) if sites else "None"

            paml_summary_rows.append({
                "Model": "M8: β&ω",
                "np": m8_parsed["np"],
                "lnL": round(m8_parsed["lnL"], 2),
                "Parameters": m8_param_str,
                "Positively selected sites": sites_str
            })
            if m7m8_result["p"] < 0.05:
                print("LRT (M7 vs M8) is significant at the 0.05/0.01 level: positive selection sites detected.")
            else:
                print("LRT (M7 vs M8) is not significant: no positive selection sites detected.")
        else:
            print("LRT (M0 vs M3) is not significant: no strong evidence of site-specific ω heterogeneity.")
            rst_path = None

        cols = ["Model", "np", "lnL", "Parameters", "Positively selected sites"]
        if paml_summary_rows:
            df = pd.DataFrame(paml_summary_rows)[cols]
            df.to_csv(paml_output / "SiteTestSummary.csv", sep='\t',index=False)
    else:
        rst_path= None

    if skip_plot:
        png_path = None
    else:
        plots_dir = gene_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        outpng_path = plots_dir / f"{gene_name}EvolutionarySites.png"
        png_path = Visualization(paml_treefile, 
                                heatmap_df, 
                                outpng_path, 
                                rst_path=rst_path, 
                                label_max_chars=label_max_chars, 
                                conservation_threshold=conservation_threshold, 
                                run_paml=run_paml,
                                outgroups=outgroups)
    return {
        "gene_dir": str(gene_dir),
        "file_input": str(file_input),
        "evo_output": str(evo_output),
        "paml_output": str(paml_output),
        "plots_dir": str(plots_dir) if not skip_plot else None,
        "logs_dir": str(logs_dir),
        "png": str(png_path) if png_path else None,
        "entropy_csv": str(entropy_csv),
        "tree": str(paml_treefile),
        "rst": str(rst_path) if rst_path else None,
        "site_test_summary": str(paml_output / "site_test_summary.csv") if run_paml else None,    
    }

