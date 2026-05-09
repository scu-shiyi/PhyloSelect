from pathlib import Path
from phyloselect.docking.AutoDocking import AutoDocking

def run(args):
    """Thin wrapper around pipeline.AutoDocking so CLI / API share the same entry."""
    mapping_csv       = str(Path(args.mapping_csv).resolve())
    tree_path         = str(Path(args.tree_path).resolve())
    output_dir        = str(Path(args.output_dir).resolve())
    outgroups = args.outgroups if args.outgroups else []
    AutoDocking(
        mapping_csv=mapping_csv,
        tree_path=tree_path,
        output_dir=output_dir,
        outgroups=outgroups)
    return None


if __name__ == "__main__":
    run()