"""Command-line entry point."""

import argparse

from .study import run


def main():
    parser = argparse.ArgumentParser(description="Reproduce the HeatStart diffusion study")
    parser.add_argument("--output", default="examples/results", help="output directory")
    args = parser.parse_args()
    result = run(args.output)
    print(
        f"Results written to {args.output}; "
        f"CN first minimum = {result['pulse']['methods']['cn']['first_min']:.6f}"
    )


if __name__ == "__main__":
    main()
