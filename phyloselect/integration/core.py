from phyloselect.integration.integration import run_absrel_env_pipeline
from pathlib import Path

def RunAbsrelEnv(args):
	alignment = str(Path(args.seq).resolve())
	tree = str(Path(args.tree).resolve())
	env_csv = str(Path(args.env).resolve())

	gene_name = Path(alignment).stem
	outdir = Path(args.output_dir).resolve() / gene_name
	outdir.mkdir(parents=True, exist_ok=True)


	result = run_absrel_env_pipeline(
		alignment=alignment,
		tree=tree,
		env_csv=env_csv,
		outdir=outdir,
		code=args.code,
		branches=args.branches,
		multiple_hits=args.multiple_hits,
		srv=args.srv,
		timeout_sec=args.timeout_sec,
		log_transform=not args.no_log_transform,
		alpha=args.alpha,
		min_n=args.min_n,
	)

	return result
