import argparse
import os
import tarfile
import shutil
import subprocess
import shlex
import zlib
from urllib.parse import urljoin


def main():
    parser = argparse.ArgumentParser(
        prog="phyloselect",
        description="PhyloSelect: A comprehensive toolkit for phylogenetic analysis and evolutionary selection.",
        formatter_class=argparse.HelpFormatter
    )

    subparsers = parser.add_subparsers(dest="command", required=True)


    # GeneMiner
    COMMAND_HELP = '''
    filter    Reference-based filtering of raw reads
    refilter  Refinement of filtered reads
    assemble  Gene assembly using wDBG
    consensus Consensus generation on heterozygous sites
    trim      Flank sequence removal
    combine   Gene alignment, concatenation and cleanup choose from geneminer
    tree      Phylogenetic tree reconstruction
    '''
    parser_geneminer = subparsers.add_parser("geneminer", help="Extract phylogenetic marker genes(CDS) from genomic data for evolutionary studies.",
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_geneminer.add_argument('actions',
                        choices=('filter', 'assemble', 'consensus', 'trim', 'combine', 'tree'),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='action',
                        nargs='*')

    parser_geneminer.add_argument('-f', help='Sample list file', metavar='FILE', required=True)
    parser_geneminer.add_argument('-r', help='Reference directory', metavar='DIR', required=True)
    parser_geneminer.add_argument('-o', dest = "output_dir", help='Output directory', metavar='DIR', required=True)
    parser_geneminer.add_argument('-p', default=1, help='Number of parallel processes', metavar='INT', type=int)

    parser_geneminer.add_argument('-kf', default=31, help='Filter k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('-ka', default=0, help='Assembly k-mer size (default = auto)', metavar='INT', type=int)
    parser_geneminer.add_argument('-s', '--step-size', default=4, help='Filter step size', metavar='INT', type=int)
    parser_geneminer.add_argument('-e', '--error-threshold', default=2, help='Error threshold', metavar='INT', type=int)
    parser_geneminer.add_argument('-sb', '--soft-boundary', choices=('0', 'auto', 'unlimited'), default='auto', help='Soft boundary (default = auto)', type=str)
    parser_geneminer.add_argument('-i', '--iteration', default=4096, help='Search depth', metavar='INT', type=int)

    parser_geneminer.add_argument('-c', '--consensus-threshold', default='0.75', help='Consensus threshold (default = 0.75)', metavar='FLOAT', type=float)

    parser_geneminer.add_argument('-ts', '--trim-source', choices=('assembly', 'consensus'), default=None, help='Whether to trim the primary assembly or the consensus sequence (default = output of last step, assembly if no other command given)')
    parser_geneminer.add_argument('-tm', '--trim-mode', choices=('all', 'longest', 'terminal', 'isoform'), default='terminal', help='Trim mode (default = terminal)', type=str)
    parser_geneminer.add_argument('-tr', '--trim-retention', default=0, help='Retention length threshold (default = 0.0)', metavar='FLOAT', type=float)

    parser_geneminer.add_argument('-cs', '--combine-source', choices=('assembly', 'consensus', 'trimmed'), default=None, help='Whether to combine the primary assembly, the consensus sequences or the trimmed sequences (default = output of last step, assembly if no other command given)')
    parser_geneminer.add_argument('-cd', '--clean-difference', default=1, help='Maximum acceptable pairwise difference in an alignment (default = 1.0)', metavar='FLOAT', type=float)
    parser_geneminer.add_argument('-cn', '--clean-sequences', default=0, help='Number of sequences required in an alignment (default = 0)', metavar='INT', type=int)

    parser_geneminer.add_argument('-m', '--tree-method', choices=('coalescent', 'concatenation'), default='coalescent', help='Multi-phyloselect tree reconstruction method (default = coalescent)')
    parser_geneminer.add_argument('-b', '--bootstrap', default=1000, help='Number of bootstrap replicates', metavar='INT', type=int)

    parser_geneminer.add_argument('--max-reads', default=0, help='Maximum reads per file', metavar='INT', type=int)
    parser_geneminer.add_argument('--min-depth', default=50, help='Minimum acceptable depth during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-depth', default=768, help='Maximum acceptable depth during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-size', default=6, help='Maximum file size during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--min-ka', default=21, help='Minimum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-ka', default=51, help='Maximum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('--msa-program', choices=('clustalo', 'mafft', 'muscle'), default='mafft', help='Program for multiple sequence alignment', type=str)
    parser_geneminer.add_argument('--no-alignment', action='store_true', help='Do not perform multiple sequence alignment')
    parser_geneminer.add_argument('--no-trimal', action='store_true', default=False, help='Do not run trimAl on alignments')
    parser_geneminer.add_argument('--phylo-program', choices=('raxmlng', 'iqtree', 'fasttree', 'veryfasttree'), default='fasttree', help='Program for phylogenetic tree reconstruction', type=str)

    # siteview
    parser_site = subparsers.add_parser("siteview", 
                                        help="Perform phylogenetic analysis to study site-wise heterogeneity, positive selection, and entropy in evolutionary trees.",
                                        formatter_class=argparse.HelpFormatter)
    parser_site.add_argument("-s", "--seq", required=True, metavar="FILE", help="Input CDS FASTA file (use --codon-aligned if already aligned).")
    parser_site.add_argument("-o", "--output-dir", dest="output_dir", required=True, help="Output directory.")
    parser_site.add_argument("-t", "--tree", default=None, metavar="FILE", help="Phylogenetic tree in Newick format (use --infer-tree if not provided).")
    parser_site.add_argument("-g", "--outgroups", nargs="+", default=None, metavar="SPECIES", help="One or more outgroup taxa used for rooting.")
    
    parser_site.add_argument("-n", "--label-max-chars", type=int, default=20, metavar="INT", help="Maximum number of characters shown for each leaf label in the tree.")
    parser_site.add_argument("--conservation-threshold", dest="conservation_threshold",type=float, default=1.4, metavar="FLOAT", help="Entropy threshold used to distinguish conserved and poorly conserved sites.")
    
    parser_site.add_argument("--no-run-paml", action="store_false", dest="run_paml",help="Skip codeml site-model analyses and only perform Evo2-based scoring.")
    parser_site.add_argument("--entropy-file", dest="entropy_csv", default=None, metavar="FILE",help="Precomputed entropy matrix file, mainly for debugging or reproducing previous results.")
    parser_site.add_argument("-p","--threads",type=int,default=1,metavar="INT",help="Total number of threads used for codeml model runs and related analysis steps.")
    parser_site.add_argument("--force",action="store_true",help="Overwrite the output directory if it already exists.")
    parser_site.add_argument("--keep-intermediate",dest="keep_intermediate",action="store_true",help="Keep temporary and intermediate files generated during the run.")
    parser_site.add_argument("--infer-tree",action="store_true",help="Infer a phylogenetic tree automatically when no tree file is provided.")
    parser_site.add_argument("--codon-aligned",action="store_true",help="Treat the input FASTA file as an already codon-aligned sequence file and skip alignment-related preprocessing.")
    parser_site.add_argument("--skip-plot",action="store_true",help="Skip figure generation and output only tabular and model-based results.")


    # selection

    parser_dnds = subparsers.add_parser("selection", help="Analyze evolutionary selection on genetic sequences to identify positive selection and evolutionary trends.",
                                        formatter_class=argparse.HelpFormatter)
    parser_dnds.add_argument("-i", "--input", dest="fasta_input", required=True, metavar="PATH", help="Input FASTA directory or a semicolon-separated list of FASTA files.")
    parser_dnds.add_argument("-o", "--output-dir", dest="output_dir", required=True, metavar="DIR", help="Output directory.")
    parser_dnds.add_argument("-m", "--tree-file-map",dest="tree_map",metavar="FILE",default=None, help="Optional file mapping each FASTA file to its corresponding phylogenetic tree.")
    # other parameters
    parser_dnds.add_argument("-p","--threads",type=int,default=1,metavar="INT",help="Total number of threads used for codeml model runs and related analysis steps.")
    parser_dnds.add_argument("--force",action="store_true",help="Overwrite the output directory if it already exists.")
    parser_dnds.add_argument("--keep-intermediate",dest="keep_intermediate",action="store_true",help="Keep temporary and intermediate files generated during the run.")
    parser_dnds.add_argument("--skip-plot",action="store_true",help="Skip figure generation and output only tabular and model-based results.")



    # env——factor
    parser_envassoc = subparsers.add_parser("envassoc", help="Perform branch-level selection–environment association analysis using aBSREL and PGLS.",
                                        formatter_class=argparse.HelpFormatter)
    parser_envassoc.add_argument("-s", '--seq', required=True, metavar="FILE",
                                 help="Codon-aligned sequence file used for aBSREL analysis.")
    parser_envassoc.add_argument("-t", '--tree', required=True, metavar="FILE",
                                 help="Phylogenetic tree in Newick format matching the input alignment.")
    parser_envassoc.add_argument("-e", "--env", required=True, metavar="FILE",
                                 help="Environmental metadata table in CSV format.")
    parser_envassoc.add_argument("-o", "--output-dir", dest="output_dir", required=True, metavar="DIR",
                                 help="Output directory.")
    # absrel
    parser_envassoc.add_argument("--code", default="Universal",metavar="CODE",help="Genetic code used by HyPhy aBSREL (default: Universal).")
    parser_envassoc.add_argument("--branches", default="All", metavar="BRANCH_SET",
                                 help="Branch set to be tested in aBSREL (default: All).")
    parser_envassoc.add_argument("--multiple-hits", dest="multiple_hits", default="None", metavar="MODEL",
                                 help="Multiple-hit substitution model used by aBSREL (default: None).")
    parser_envassoc.add_argument("--srv", default="No", metavar="YES/NO",
                                 help="Enable synonymous rate variation in aBSREL (default: No).")
    parser_envassoc.add_argument("--timeout-sec", dest="timeout_sec", type=int, default=None, metavar="SECONDS",
                                 help="Maximum runtime allowed for aBSREL.")
    # pgls
    parser_envassoc.add_argument("--no-log-transform", dest="no_log_transform", action="store_true",
                                 help="Disable log1p transformation of omega_weighted before PGLS analysis.")
    parser_envassoc.add_argument("--alpha", type=float, default=0.05, metavar="FLOAT",
                                 help="Significance threshold for FDR-adjusted Q-values (default: 0.05).")
    parser_envassoc.add_argument("--min-n", dest="min_n", type=int, default=5, metavar="INT",
                                 help="Minimum number of taxa required for each environmental variable.")
    parser_envassoc.add_argument( "-p", "--threads", type=int, default=1, metavar="INT", help="Number of threads used for analysis.")
    parser_envassoc.add_argument("--force", action="store_true", help="Overwrite the output directory if it already exists.")
 
    # Docking 
    parser_dock = subparsers.add_parser("docking", help="Simulate molecular interactions and evaluate ligand binding affinities in enzyme-substrate systems.",
                                        formatter_class=argparse.HelpFormatter)
    parser_dock.add_argument("-c", "--dock-config",dest = 'mapping_csv', required=True, metavar="FILE", help="Docking configuration file in CSV format.")
    parser_dock.add_argument("-t", "--tree", dest='tree_path', required=True, metavar="FILE", help="Phylogenetic tree in Newick format.")
    parser_dock.add_argument("-o", "--output-dir", dest = 'output_dir', required=True, metavar="DIR", help="Output directory.")

    parser_dock.add_argument("-g", "--outgroups", nargs="+", default=None, metavar="SPECIES", help="One or more outgroup taxa used for rooting.")
    parser_dock.add_argument("--force", action="store_true", help="Overwrite the output directory if it already exists.")
    parser_dock.add_argument("--keep-intermediate", dest="keep_intermediate", action="store_true", help="Keep temporary and intermediate files.")

    parser_check = subparsers.add_parser(
        "check",
        help="Check whether required third-party tools are available.",
        formatter_class=argparse.HelpFormatter
    )


    args = parser.parse_args()

    if args.command == "siteview":
        os.makedirs(args.output_dir, exist_ok=True)
        from phyloselect.evolution.core import RunSite
        RunSite(args)
    elif args.command == "selection":
        os.makedirs(args.output_dir, exist_ok=True)
        from phyloselect.evolution.core import RunBranch
        RunBranch(args)
    elif args.command == "docking":
        os.makedirs(args.output_dir, exist_ok=True)
        from phyloselect.docking.core import run as run_docking
        run_docking(args)
    elif args.command == "envassoc":
        os.makedirs(args.output_dir, exist_ok=True)
        from phyloselect.integration.core import RunAbsrelEnv
        RunAbsrelEnv(args)
    elif args.command == "geneminer":
        os.makedirs(args.output_dir, exist_ok=True)
        from phyloselect.seq_mining.core import run as run_geneminer

        geneminer_argv = [
            "-f", args.f,
            "-r", args.r,
            "-o", args.output_dir,
            "-p", str(args.p),
            "-kf", str(args.kf),
            "-ka", str(args.ka),
            "-s", str(args.step_size),
            "-e", str(args.error_threshold),
            "-sb", str(args.soft_boundary),
            "-i", str(args.iteration),
            "-c", str(args.consensus_threshold),
            "-tm", args.trim_mode,
            "-tr", str(args.trim_retention),
            "-cd", str(args.clean_difference),
            "-cn", str(args.clean_sequences),
            "-m", args.tree_method,
            "-b", str(args.bootstrap),
            "--max-reads", str(args.max_reads),
            "--min-depth", str(args.min_depth),
            "--max-depth", str(args.max_depth),
            "--max-size", str(args.max_size),
            "--min-ka", str(args.min_ka),
            "--max-ka", str(args.max_ka),
            "--msa-program", args.msa_program,
            "--phylo-program", args.phylo_program,
        ]

        if args.trim_source is not None:
            geneminer_argv += ["-ts", args.trim_source]

        if args.combine_source is not None:
            geneminer_argv += ["-cs", args.combine_source]

        if args.no_alignment:
            geneminer_argv.append("--no-alignment")

        if args.no_trimal:
            geneminer_argv.append("--no-trimal")

        if args.actions:
            geneminer_argv += args.actions

        run_geneminer(geneminer_argv)
    elif args.command == "check":
        import shutil
        import sys
        def check_tool(name, aliases=None):
            if aliases is None:
                aliases = [name]

            for cmd in aliases:
                path = shutil.which(cmd)
                if path:
                    return path

            return None

        tools = {
                "codeml (PAML)": ["codeml"],
                "hyphy": ["hyphy"],
                "iqtree": ["iqtree", "iqtree2"],
                "muscle": ["muscle"],
                "trimal": ["trimal"],
                "pal2nal": ["pal2nal.pl"],
                "perl": ["perl"],
                "fpocket": ["fpocket"],
                "openbabel": ["obabel", "openbabel"],
                "autodock-vina": ["vina", "autodock-vina"],
            }

        all_ok = True
        false_tool = []
        for tool_name, aliases in tools.items():
            path = check_tool(tool_name, aliases)
            if not path:
                false_tool.append(tool_name)
                all_ok = False
        if all_ok:
            print("[Ready] All dependencies are satisfied!")
        else:
            print("[Warning] Some dependencies are missing. Missing dependencies:")
            for t in false_tool:
                print(f" - {t}")
            print("\nPlease install them before running the pipeline.")
            sys.exit(1)

if __name__ == "__main__":
    main()
