
import argparse
import typing
import sys


class Args(typing.NamedTuple):
    filename: str
    alignment_mode: typing.Literal["unaligned", "aligned"]
    gap_open: float
    gap_extension: float
    molecule: str


def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("-f","--filename", type=str, required=True, help="Input file with sequences to align")
    parser.add_argument("-a","--alignment_mode", type=str, choices=("unaligned", "aligned"), required=True, 
                        help='Choose the format of sequences: {unaligned, aligned}')
    parser.add_argument("--gap-open", type=float, required=True, help="Penalty for gap opening")
    parser.add_argument("--gap-extension", type=float, required=True, help="Penalty for gap extension")
    parser.add_argument("-m","--molecule", type=str, choices=('DNA', 'protein'), required=True, 
                        help='Choose the type of sequences: {DNA, protein}')
    return parser


def parse_args() -> Args:
    parser = create_parser()
    args = parser.parse_args(sys.argv[1:])

    return Args(
        filename=args.filename,
        alignment_mode=args.alignment_mode,
        gap_open=args.gap_open,
        gap_extension=args.gap_extension,
        molecule=args.molecule
    )

