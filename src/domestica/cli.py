import argparse
import logging
from pathlib import Path

# Assuming flat local directory imports structure.
from domestica.core import run_pipeline


def setup_logging(verbose: bool):
    """Sets up terminal output. Verbose mode shows DEBUG info."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S"
    )


def build_parser() -> argparse.ArgumentParser:
    """Configures the command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Protein Designer Pipeline: Excel/FASTA -> Params -> DNA Opt -> Excel"
    )

    # Required File IO
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input fasta (.fasta) or Excel file (.xlsx)")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output Excel file (.xlsx)")

    # Feature toggles
    parser.add_argument(
        "--params",
        nargs="+",
        default=[],
        choices=["mw", "pi", "net_charge", "gravy", "length", "all"],
        help="Select which protein parameters to calculate. Use 'all' for all parameters."
    )
    parser.add_argument("--optimize", action="store_true",
                        help="Perform codon optimization to generate nucleotide sequences")

    # Optional parameters
    parser.add_argument("--name-col", default="Name", help="Column header for protein names (Default: Name)")
    parser.add_argument("--seq-col", default="Sequence", help="Column header for protein sequences (Default: Sequence)")
    parser.add_argument("-n", "--nstruct", type=int, default=10,
                        help="Number of optimized DNA structures to generate per input (Default: 1)")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for net charge calculation (Default: 7.4)")
    parser.add_argument("--n-tag", type=str, default="",
                        help="Amino acid sequence to append to the N-terminus (used for analysis ONLY)")
    parser.add_argument("--c-tag", type=str, default="",
                        help="Amino acid sequence to append to the C-terminus (used for analysis ONLY)")


    # Optional DNA optimization and IDT arguments
    parser.add_argument("--idt_credentials_dir", type=str, default="~/.idt_credentials",
                        help="A path to the place to search for you stored IDT API credentials. If no info.json file is found, then you will be prompted to enter new ones and they will be stored there")
    parser.add_argument("-v", "--vector", type=Path, help="Genbank vector file for insertion")
    parser.add_argument("--skip-idt", action="store_true", help="Skip IDT complexity checking API calls")
    parser.add_argument('--idt_type', type=str, help='type of sequence to query', default='gene',
                        choices=['gene', 'gblock', 'gblock_hifi', 'eblock', 'old'])
    parser.add_argument("--idt_threshold", type=float, default=7,
                        help="automatically accept the first solution with IDT score under this threshold")
    # Utilities
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)

    logging.info(f"Starting pipeline...")
    logging.info(f"Reading sequences from {args.input}...")

    if not args.params and not args.optimize:
        logging.warning("Neither --params nor --optimize was selected. The output will just mirror the input.")

    # Pass configuration to core pipeline
    run_pipeline(
        input_path=args.input,
        output_path=args.output,
        params=args.params,
        optimize=args.optimize,
        vector_path=args.vector,
        nstruct=args.nstruct,
        skip_idt=args.skip_idt,
        ph=args.ph,
        name_col=args.name_col,
        seq_col=args.seq_col,
        idt_type=args.idt_type,
        idt_credentials_dir=args.idt_credentials_dir,
        idt_threshold=args.idt_threshold,
        n_tag=args.n_tag,
        c_tag=args.c_tag
    )


if __name__ == "__main__":
    main()