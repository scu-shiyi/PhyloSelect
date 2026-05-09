from phyloselect.evolution.runsite import Run_Site
from phyloselect.evolution.RunEvoDnDs import RunEvoDnDs
from pathlib import Path


def RunSite(args):
    fasta_path = str(Path(args.seq).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    tree_path = str(Path(args.tree).resolve()) if args.tree else None

    outgroups = args.outgroups if args.outgroups else []
    label_max_chars = args.label_max_chars if args.label_max_chars is not None else 20
    conservation_threshold = (args.conservation_threshold if args.conservation_threshold is not None else 1.4)
    entropy_csv = str(Path(args.entropy_csv).resolve()) if args.entropy_csv else None
    run_paml = args.run_paml
    threads = args.threads if args.threads else 1
    force = args.force if hasattr(args, "force") else False
    keep_intermediate = args.keep_intermediate

    infer_tree = args.infer_tree
    codon_aligned = args.codon_aligned
    skip_plot = args.skip_plot

    output_path = Run_Site(
        fasta_path=fasta_path,
        output_dir=output_dir,
        tree_path=tree_path,
        outgroups=outgroups,
        label_max_chars=label_max_chars,
        conservation_threshold=conservation_threshold,
        entropy_csv=entropy_csv,
        run_paml=run_paml,
        threads=threads,
        force=force,
        keep_intermediate=keep_intermediate,
        infer_tree = infer_tree,
        codon_aligned = codon_aligned,
        skip_plot = skip_plot,
        )

    return output_path

def RunBranch(args):
    fasta_input = str(Path(args.fasta_input).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    fasta_tree_map = str(Path(args.tree_map).resolve()) if args.tree_map else None

    threads = args.threads if args.threads else 1
    force = args.force if hasattr(args, "force") else False
    keep_intermediate = args.keep_intermediate
    skip_plot = args.skip_plot
    
    output_png = RunEvoDnDs(
        fasta_input=fasta_input,
        output_dir=output_dir,
        fasta_tree_map=fasta_tree_map,
        threads=threads,
        force=force,
        keep_intermediate=keep_intermediate,
        skip_plot=skip_plot,
    )
    return output_png