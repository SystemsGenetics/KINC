"""
Command-line interface for KINC.
"""

import argparse
import sys
from .analytics import import_emx, export_emx, Similarity, extract_network


def create_parser():
    """Create argument parser for KINC CLI."""
    parser = argparse.ArgumentParser(
        description='KINC - Knowledge Independent Network Construction',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # import-emx
    import_parser = subparsers.add_parser('import-emx',
                                          help='Import gene expression matrix from text file')
    import_parser.add_argument('--input', required=True, help='Input text file')
    import_parser.add_argument('--output', required=True, help='Output HDF5 file')
    import_parser.add_argument('--nan', default='NA', help='Token for missing values (default: NA)')
    import_parser.add_argument('--no-header', action='store_true', help='Input has no header row')
    import_parser.add_argument('--no-rownames', action='store_true', help='Input has no row names')

    # export-emx
    export_parser = subparsers.add_parser('export-emx',
                                          help='Export gene expression matrix to text file')
    export_parser.add_argument('--input', required=True, help='Input HDF5 file')
    export_parser.add_argument('--output', required=True, help='Output text file')
    export_parser.add_argument('--nan', default='NA', help='Token for missing values (default: NA)')

    # similarity
    sim_parser = subparsers.add_parser('similarity',
                                       help='Compute pairwise gene correlations')
    sim_parser.add_argument('--input', required=True, help='Input expression matrix (HDF5)')
    sim_parser.add_argument('--cmx', required=True, help='Output correlation matrix (HDF5)')
    sim_parser.add_argument('--ccm', required=True, help='Output cluster composition matrix (HDF5)')
    sim_parser.add_argument('--corrmethod', default='pearson', choices=['pearson', 'spearman'],
                           help='Correlation method (default: pearson)')
    sim_parser.add_argument('--clusmethod', default='none', choices=['none', 'gmm'],
                           help='Clustering method (default: none)')
    sim_parser.add_argument('--minexpr', type=float, default=float('-inf'),
                           help='Minimum expression threshold')
    sim_parser.add_argument('--maxexpr', type=float, default=float('inf'),
                           help='Maximum expression threshold')
    sim_parser.add_argument('--minsamp', type=int, default=30,
                           help='Minimum samples for correlation (default: 30)')
    sim_parser.add_argument('--minclus', type=int, default=1,
                           help='Minimum clusters to test (default: 1)')
    sim_parser.add_argument('--maxclus', type=int, default=5,
                           help='Maximum clusters to test (default: 5)')
    sim_parser.add_argument('--crit', default='ICL', choices=['AIC', 'BIC', 'ICL'],
                           help='Model selection criterion (default: ICL)')
    sim_parser.add_argument('--preout', type=lambda x: x.lower() == 'true', default=True,
                           help='Remove pre-clustering outliers (default: true)')
    sim_parser.add_argument('--postout', type=lambda x: x.lower() == 'true', default=True,
                           help='Remove post-clustering outliers (default: true)')
    sim_parser.add_argument('--mincorr', type=float, default=0.5,
                           help='Minimum |correlation| to save (default: 0.5)')
    sim_parser.add_argument('--maxcorr', type=float, default=1.0,
                           help='Maximum |correlation| to save (default: 1.0)')
    sim_parser.add_argument('--jobs', type=int, default=1,
                           help='Number of parallel jobs (default: 1)')

    # extract
    extract_parser = subparsers.add_parser('extract',
                                          help='Extract network edges to text file')
    extract_parser.add_argument('--input', required=True, help='Input correlation matrix (HDF5)')
    extract_parser.add_argument('--output', required=True, help='Output network file')
    extract_parser.add_argument('--ccm', help='Cluster composition matrix (HDF5, optional)')
    extract_parser.add_argument('--mincorr', type=float, default=0.0,
                               help='Minimum |correlation| threshold (default: 0.0)')
    extract_parser.add_argument('--maxcorr', type=float, default=1.0,
                               help='Maximum |correlation| threshold (default: 1.0)')
    extract_parser.add_argument('--format', default='text', choices=['text', 'tidy', 'minimal'],
                               help='Output format (default: text)')

    return parser


def main():
    """Main entry point for KINC CLI."""
    parser = create_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        if args.command == 'import-emx':
            import_emx(
                input_file=args.input,
                output_file=args.output,
                nan_token=args.nan,
                has_header=not args.no_header,
                has_row_names=not args.no_rownames
            )

        elif args.command == 'export-emx':
            export_emx(
                input_file=args.input,
                output_file=args.output,
                nan_token=args.nan
            )

        elif args.command == 'similarity':
            sim = Similarity(
                input_file=args.input,
                cmx_file=args.cmx,
                ccm_file=args.ccm,
                correlation_method=args.corrmethod,
                clustering_method=args.clusmethod,
                min_expression=args.minexpr,
                max_expression=args.maxexpr,
                min_samples=args.minsamp,
                min_clusters=args.minclus,
                max_clusters=args.maxclus,
                criterion=args.crit,
                preout=args.preout,
                postout=args.postout,
                min_correlation=args.mincorr,
                max_correlation=args.maxcorr
            )
            sim.run(n_jobs=args.jobs)

        elif args.command == 'extract':
            extract_network(
                cmx_file=args.input,
                output_file=args.output,
                ccm_file=args.ccm,
                min_correlation=args.mincorr,
                max_correlation=args.maxcorr,
                output_format=args.format
            )

        else:
            print(f"Unknown command: {args.command}")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
