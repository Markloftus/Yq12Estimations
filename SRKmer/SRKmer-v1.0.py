#!/usr/bin/env python3

import argparse
import json
import sys
import time
from collections import defaultdict
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
DEFAULT_DYZ1_UNIT_LENGTH = 3569.0
DEFAULT_DYZ2_UNIT_LENGTH = 2420.0


def positive_float(value: str) -> float:
    """Argparse type requiring a finite value greater than zero."""
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected a number greater than zero, received: {value}"
        ) from exc

    if not (number > 0.0) or number == float("inf"):
        raise argparse.ArgumentTypeError(
            f"Expected a finite number greater than zero, received: {value}"
        )

    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a CRAM or BAM file for Yq12, Yq12-subunit, "
            "Y-centromere, and DYZ3 24-mers, then calculate an "
            "SRKmer Yq12 length estimate."
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
        "-j",
        "--classification-json",
        required=True,
        type=Path,
        help=(
            "JSON dictionary mapping Yq12 subunit 24-mers to classifications "
            "such as DYZ1 or DYZ2. For example: "
            "my_95Percent_Identity_SubunitKmers.json"
        ),
    )
    parser.add_argument(
        "--sample-depth",
        required=True,
        type=positive_float,
        help=(
            "Sequencing-depth value used to normalize the summed read lengths. "
            "This must be greater than zero."
        ),
    )
    parser.add_argument(
        "--dyz1-unit-length",
        type=positive_float,
        default=DEFAULT_DYZ1_UNIT_LENGTH,
        help=(
            "DYZ1 repeat-unit length used to estimate the DYZ1 subunit count "
            f"(default: {DEFAULT_DYZ1_UNIT_LENGTH:g} bp)."
        ),
    )
    parser.add_argument(
        "--dyz2-unit-length",
        type=positive_float,
        default=DEFAULT_DYZ2_UNIT_LENGTH,
        help=(
            "DYZ2 repeat-unit length used to estimate the DYZ2 subunit count "
            f"(default: {DEFAULT_DYZ2_UNIT_LENGTH:g} bp)."
        ),
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
    """Load one 24-mer per line and initialize counters expected by Cython."""
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


def load_classifications(path: Path) -> dict[str, str]:
    """
    Load the k-mer classification JSON.

    Expected format:
        {
            "AAAAAAAA...": "DYZ1",
            "CCCCCCCC...": "DYZ2"
        }
    """
    with path.open("r", encoding="utf-8") as handle:
        classifications = json.load(handle)

    if not isinstance(classifications, dict):
        raise ValueError(
            "Classification JSON must contain a dictionary mapping "
            "k-mers to classifications."
        )

    cleaned = {}
    for raw_kmer, raw_classification in classifications.items():
        kmer = str(raw_kmer).strip().upper()
        classification = str(raw_classification).strip().upper()

        if len(kmer) != KMER_SIZE:
            raise ValueError(
                "Classification JSON contains a key with length "
                f"{len(kmer)} rather than {KMER_SIZE}: {raw_kmer!r}"
            )

        invalid = set(kmer) - set("ACGT")
        if invalid:
            raise ValueError(
                "Classification JSON contains a k-mer with non-ACGT "
                f"characters: {raw_kmer!r}"
            )

        if not classification:
            raise ValueError(
                f"Classification JSON contains an empty class for k-mer {kmer}."
            )

        cleaned[kmer] = classification

    if not cleaned:
        raise ValueError(f"Classification JSON contains no entries: {path}")

    return cleaned


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


def parse_read_length(read_identifier: str) -> int:
    """
    Parse the length from the Cython output format:
        read_name_<length>

    rsplit is used because read names may themselves contain underscores.
    """
    try:
        length_text = str(read_identifier).rsplit("_", 1)[1]
        length = int(length_text)
    except (IndexError, ValueError) as exc:
        raise ValueError(
            "Could not parse a read length from the Yq12-subunit record "
            f"{read_identifier!r}. Expected a value ending in '_<length>'."
        ) from exc

    if length < 0:
        raise ValueError(
            f"Read length cannot be negative in record: {read_identifier!r}"
        )

    return length


def calculate_sample_lengths(
    sample_name: str,
    yq12_sub_reads: dict,
    classifications: dict[str, str],
    sample_depth: float,
    dyz1_unit_length: float,
    dyz2_unit_length: float,
) -> tuple[dict, dict]:
    """
    Reproduce the downstream notebook calculation directly from the
    in-memory Yq12-subunit read dictionary.

    Unpolished:
        Sum all read lengths classified as DYZ1 or DYZ2, then divide by depth.

    Polished:
        Deduplicate exact (k-mer, read_identifier) rows, sum their read
        lengths, then divide by depth. The original notebook used weight 1.0.

    DYZ1_Counts and DYZ2_Counts:
        Divide each depth-normalized class total by its repeat-unit length.
    """
    per_class_bp = defaultdict(int)
    unique_rows = set()
    polished_total_bp = 0
    total_records = 0
    unclassified_records = 0

    for raw_kmer, read_identifiers in yq12_sub_reads.items():
        if raw_kmer in {"total_reads", "lengths"}:
            continue

        kmer = str(raw_kmer).upper()
        classification = classifications.get(kmer, "NO")

        for read_identifier in read_identifiers:
            read_length = parse_read_length(read_identifier)
            total_records += 1
            per_class_bp[classification] += read_length

            if classification == "NO":
                unclassified_records += 1

            row = (kmer, str(read_identifier))
            if row not in unique_rows:
                unique_rows.add(row)
                polished_total_bp += read_length

    dyz1_total_bp = int(per_class_bp.get("DYZ1", 0))
    dyz2_total_bp = int(per_class_bp.get("DYZ2", 0))

    sample_lengths = {
        sample_name: {
            "Polished": polished_total_bp / sample_depth,
            "Unpolished": (dyz1_total_bp + dyz2_total_bp) / sample_depth,
            "DYZ1_Counts": (dyz1_total_bp / sample_depth) / dyz1_unit_length,
            "DYZ2_Counts": (dyz2_total_bp / sample_depth) / dyz2_unit_length,
        }
    }

    calculation_details = {
        "Sample": sample_name,
        "Sample_Depth": sample_depth,
        "DYZ1_Unit_Length_bp": dyz1_unit_length,
        "DYZ2_Unit_Length_bp": dyz2_unit_length,
        "Total_Yq12_Subunit_Records": total_records,
        "Unique_Yq12_Subunit_Records": len(unique_rows),
        "Unclassified_Yq12_Subunit_Records": unclassified_records,
        "DYZ1_Total_Read_Length_bp": dyz1_total_bp,
        "DYZ2_Total_Read_Length_bp": dyz2_total_bp,
        "Polished_Total_Read_Length_bp": polished_total_bp,
        "Per_Class_Read_Length_bp": dict(per_class_bp),
    }

    return sample_lengths, calculation_details


def write_json(output_path: Path, data: dict) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=4, sort_keys=True)
        handle.write("\n")


def write_summary(
    output_path: Path,
    append: bool,
    total_reads: int,
    yq12_reads: dict,
    cent_reads: dict,
    yq12_sub_reads: dict,
    dyz3_reads: dict,
    sample_depth: float,
    estimates: dict,
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
        "Sample_Depth",
        "Polished",
        "Unpolished",
        "DYZ1_Counts",
        "DYZ2_Counts",
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
        sample_depth,
        estimates["Polished"],
        estimates["Unpolished"],
        estimates["DYZ1_Counts"],
        estimates["DYZ2_Counts"],
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
        classification_path = require_file(
            args.classification_json, "Classification JSON"
        )

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
        classifications = load_classifications(classification_path)

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

        sample_lengths, calculation_details = calculate_sample_lengths(
            sample_name=sample_name,
            yq12_sub_reads=yq12_sub_reads,
            classifications=classifications,
            sample_depth=args.sample_depth,
            dyz1_unit_length=args.dyz1_unit_length,
            dyz2_unit_length=args.dyz2_unit_length,
        )

        runtime_seconds = time.time() - start_time
        estimates = sample_lengths[sample_name]

        summary_path = output_dir / f"{sample_name}_kmerSampleResults.tsv"
        subunit_reads_path = output_dir / f"{sample_name}_yq12Sub_reads.txt"
        sample_lengths_path = output_dir / f"{sample_name}_SampleLengths.json"
        details_path = output_dir / f"{sample_name}_SRKmerCalculationDetails.json"

        write_summary(
            output_path=summary_path,
            append=args.append,
            total_reads=total_reads,
            yq12_reads=yq12_reads,
            cent_reads=cent_reads,
            yq12_sub_reads=yq12_sub_reads,
            dyz3_reads=dyz3_reads,
            sample_depth=args.sample_depth,
            estimates=estimates,
            runtime_seconds=runtime_seconds,
        )
        write_yq12_subunit_reads(subunit_reads_path, yq12_sub_reads)
        write_json(sample_lengths_path, sample_lengths)
        write_json(details_path, calculation_details)

        print(f"Summary written to: {summary_path}")
        if int(yq12_sub_reads["total_reads"]) > 0:
            print(f"Yq12 subunit reads written to: {subunit_reads_path}")
        print(f"SampleLengths dictionary written to: {sample_lengths_path}")
        print(f"Calculation details written to: {details_path}")
        print(f"Completed in {runtime_seconds:.2f} seconds.")

        return 0

    except (FileNotFoundError, ValueError, OSError, pysam.utils.SamtoolsError) as exc:
        parser.exit(status=1, message=f"Error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
