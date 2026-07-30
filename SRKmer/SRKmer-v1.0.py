#!/usr/bin/env python3

import argparse
import sys
import time
from pathlib import Path

import pysam

try:
    import mainKmer
except ImportError as exc:
    raise SystemExit(
        "Could not import the compiled 'mainKmer' module.\n"
        "From the SRKmer directory, compile it first with:\n"
        "    python setup.py build_ext --inplace"
    ) from exc


KMER_SIZE = 24


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a CRAM or BAM file for Yq12, Yq12-subunit, "
            "Y-centromere, and DYZ3 24-mers."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="Input CRAM or BAM file.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        help=(
            "Reference FASTA used to create the input CRAM. "
            "Required for CRAM input; optional for BAM input."
        ),
    )
    parser.add_argument(
        "-q",
        "--yq12",
        required=True,
        type=Path,
        help="Text file containing the main Yq12 24-mers, one per line.",
    )
    parser.add_argument(
        "-s",
        "--yq12s",
        dest="yq12_subunit",
        required=True,
        type=Path,
        help="Text file containing Yq12 subunit-specific 24-mers, one per line.",
    )
    parser.add_argument(
        "-c",
        "--cent",
        dest="centromere",
        required=True,
        type=Path,
        help="Text file containing Y-centromere 24-mers, one per line.",
    )
    parser.add_argument(
        "-d",
        "--cents",
        dest="dyz3",
        required=True,
        type=Path,
        help="Text file containing DYZ3-specific 24-mers, one per line.",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=Path,
        help="Output directory. It will be created if it does not exist.",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help=(
            "Append to an existing summary TSV. By default, existing output "
            "files are overwritten."
        ),
    )

    return parser


def require_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {path}")
    return path


def load_kmer_file(path: Path, label: str) -> dict:
    """Load one 24-mer per line and initialize the counters expected by Cython."""
    kmers = {}

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            kmer = line.strip().upper()

            if not kmer:
                continue

            if len(kmer) != KMER_SIZE:
                raise ValueError(
                    f"{label}: line {line_number} has length {len(kmer)}, "
                    f"but SRKmer requires {KMER_SIZE}-mers."
                )

            invalid = set(kmer) - set("ACGT")
            if invalid:
                raise ValueError(
                    f"{label}: line {line_number} contains non-ACGT "
                    f"characters: {''.join(sorted(invalid))}"
                )

            kmers[kmer] = []

    if not kmers:
        raise ValueError(f"{label} contains no valid {KMER_SIZE}-mers: {path}")

    kmers["lengths"] = 0
    kmers["total_reads"] = 0
    return kmers


def alignment_mode(path: Path) -> str:
    suffix = path.suffix.lower()

    if suffix == ".cram":
        return "rc"
    if suffix == ".bam":
        return "rb"

    raise ValueError(
        f"Unsupported alignment format: {path.name}. "
        "The input filename must end in .cram or .bam."
    )


def sample_name_from_path(path: Path) -> str:
    return path.stem


def write_summary(
    output_path: Path,
    append: bool,
    total_reads: int,
    yq12_reads: dict,
    cent_reads: dict,
    yq12_sub_reads: dict,
    dyz3_reads: dict,
    runtime_seconds: float,
) -> None:
    columns = [
        "Total_Reads",
        "Total_Yq12_Reads",
        "Total_Yq12_Read_Length",
        "Total_Yq12_Subunit_Reads",
        "Total_Yq12_Subunit_Read_Lengths",
        "Total_Centromere_Reads",
        "Total_Centromere_Read_Lengths",
        "Total_DYZ3_Reads",
        "Total_DYZ3_Length",
        "RunTime",
    ]

    values = [
        total_reads,
        yq12_reads["total_reads"],
        yq12_reads["lengths"],
        yq12_sub_reads["total_reads"],
        yq12_sub_reads["lengths"],
        cent_reads["total_reads"],
        cent_reads["lengths"],
        dyz3_reads["total_reads"],
        dyz3_reads["lengths"],
        runtime_seconds,
    ]

    file_exists_with_content = output_path.exists() and output_path.stat().st_size > 0
    mode = "a" if append else "w"

    with output_path.open(mode, encoding="utf-8") as handle:
        if not append or not file_exists_with_content:
            handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(map(str, values)) + "\n")


def write_yq12_subunit_reads(output_path: Path, yq12_sub_reads: dict) -> None:
    if int(yq12_sub_reads["total_reads"]) == 0:
        return

    with output_path.open("w", encoding="utf-8") as handle:
        for kmer, reads in yq12_sub_reads.items():
            if kmer in {"total_reads", "lengths"}:
                continue
            for read in reads:
                handle.write(f"{kmer}\t{read}\n")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    start_time = time.time()

    try:
        input_path = require_file(args.input, "Input alignment file")
        yq12_path = require_file(args.yq12, "Yq12 k-mer file")
        yq12_subunit_path = require_file(
            args.yq12_subunit, "Yq12 subunit k-mer file"
        )
        centromere_path = require_file(args.centromere, "Centromere k-mer file")
        dyz3_path = require_file(args.dyz3, "DYZ3 k-mer file")

        mode = alignment_mode(input_path)

        reference_path = None
        if mode == "rc":
            if args.reference is None:
                parser.error("--reference is required when --input is a CRAM file.")
            reference_path = require_file(args.reference, "Reference FASTA")
        elif args.reference is not None:
            reference_path = require_file(args.reference, "Reference FASTA")

        output_dir = args.out.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        yq12 = load_kmer_file(yq12_path, "Yq12 k-mer file")
        yq12_subunit = load_kmer_file(
            yq12_subunit_path, "Yq12 subunit k-mer file"
        )
        centromere = load_kmer_file(centromere_path, "Centromere k-mer file")
        dyz3 = load_kmer_file(dyz3_path, "DYZ3 k-mer file")

        open_kwargs = {}
        if reference_path is not None:
            open_kwargs["reference_filename"] = str(reference_path)

        with pysam.AlignmentFile(str(input_path), mode, **open_kwargs) as samfile:
            (
                total_reads,
                yq12_reads,
                cent_reads,
                yq12_sub_reads,
                dyz3_reads,
            ) = mainKmer.myFunction(
                samfile,
                yq12,
                centromere,
                yq12_subunit,
                dyz3,
            )

        sample_name = sample_name_from_path(input_path)
        runtime_seconds = time.time() - start_time

        summary_path = output_dir / f"{sample_name}_kmerSampleResults.tsv"
        subunit_reads_path = output_dir / f"{sample_name}_yq12Sub_reads.txt"

        write_summary(
            output_path=summary_path,
            append=args.append,
            total_reads=total_reads,
            yq12_reads=yq12_reads,
            cent_reads=cent_reads,
            yq12_sub_reads=yq12_sub_reads,
            dyz3_reads=dyz3_reads,
            runtime_seconds=runtime_seconds,
        )
        write_yq12_subunit_reads(subunit_reads_path, yq12_sub_reads)

        print(f"Summary written to: {summary_path}")
        if int(yq12_sub_reads["total_reads"]) > 0:
            print(f"Yq12 subunit reads written to: {subunit_reads_path}")
        print(f"Completed in {runtime_seconds:.2f} seconds.")

        return 0

    except (FileNotFoundError, ValueError, OSError, pysam.utils.SamtoolsError) as exc:
        parser.exit(status=1, message=f"Error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
