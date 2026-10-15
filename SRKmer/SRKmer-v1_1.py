#!/usr/bin/env python3
import argparse
import json
import multiprocessing
import os
import signal
from concurrent.futures import ProcessPoolExecutor
import sqlite3
import tempfile
from contextlib import closing, contextmanager
import sys
import time
from collections import defaultdict, deque
from pathlib import Path
from typing import Dict

__version__ = "1.2"

MIN_PYTHON_VERSION = (3, 7)
if sys.version_info < MIN_PYTHON_VERSION:
    raise SystemExit(
        "SRKmer requires Python {}.{} or newer, but this interpreter is "
        "{}.{}. Load a newer Python module or environment and try again.".format(
            MIN_PYTHON_VERSION[0],
            MIN_PYTHON_VERSION[1],
            sys.version_info[0],
            sys.version_info[1],
        )
    )

try:
    import pysam
except ImportError as exc:
    raise SystemExit(
        "Could not import 'pysam', which SRKmer needs to read CRAM/BAM files.\n"
        "Install it with:\n"
        "    pip install pysam\n"
        "or:\n"
        "    conda install -c bioconda pysam"
    ) from exc



KMER_SIZE = 24
DEFAULT_DYZ1_UNIT_LENGTH = 3569.0
DEFAULT_DYZ2_UNIT_LENGTH = 2420.0
COUNTER_KEYS = ("total_reads", "lengths")

MAX_BATCH_PAYLOAD = 4 * 1024 * 1024
_WORKER_TARGETS = None


def log(message: str) -> None:
    """Progress and warning output goes to stderr so stdout stays parseable."""
    print(message, file=sys.stderr, flush=True)


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


def positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Expected a positive integer.") from exc
    if number < 1:
        raise argparse.ArgumentTypeError("Expected an integer >= 1.")
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
        "--version",
        action="version",
        version=f"SRKmer {__version__}",
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
            "Add a row to this sample's existing summary TSV instead of "
            "replacing it. The summary file is named after the input file, so "
            "this accumulates repeated runs of the same sample; it does not "
            "combine different samples into one file. By default the summary "
            "is overwritten."
        ),
    )
    parser.add_argument(
        "--allow-disjoint-kmers",
        action="store_true",
        help=(
            "Continue even when a subunit k-mer file shares no k-mers with its "
            "parent file (Yq12 subunit with Yq12, DYZ3 with centromere). "
            "Without this flag SRKmer stops early, because such k-mers can "
            "never be counted."
        ),
    )
    parser.add_argument(
        "--allow-unclassified",
        action="store_true",
        help=(
            "Continue even when no k-mer in the subunit file appears in the "
            "classification JSON. Without this flag SRKmer stops early, "
            "because every estimate would be zero."
        ),
    )

    parser.add_argument(
        "--temp-dir", type=Path,
        help=("Directory for the temporary on-disk read database (default: --out). "
              "Use a disk-backed directory with adequate free space. The "
              "database is removed on normal completion or handled errors."),
    )
    parser.add_argument(
        "--workers", type=positive_int, default=1,
        help=("K-mer scanning processes (default: 1, runs inline). For an "
              "8-CPU Slurm allocation, start with 7 workers and the default "
              "--io-threads 1, leaving one CPU for reading and merging."),
    )
    parser.add_argument(
        "--io-threads", type=positive_int, default=1,
        help=("pysam BAM/CRAM decoding thread setting (default: 1). Additional "
              "decoding threads need CPU and memory beyond the scan workers."),
    )
    parser.add_argument(
        "--batch-size", type=positive_int, default=2000,
        help=("Maximum alignment records per parallel task (default: 2000). "
              "Payload also targets <=4 MiB; one long record may exceed this. "
              "At most twice --workers tasks are submitted at once."),
    )
    parser.add_argument(
        "--start-method", choices=["auto", "spawn"] + (
            ["fork"] if "fork" in multiprocessing.get_all_start_methods() else []
        ),
        default="auto",
        help=("Worker process start method. Auto uses fork on Linux and spawn "
              "elsewhere. Fork initially shares lookup pages using copy-on-write; "
              "spawn copies the lookup to each worker."),
    )
    return parser


def require_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")
    try:
        with path.open("rb"):
            pass
    except OSError as exc:
        raise OSError(f"{label} cannot be read: {path} ({exc})") from exc
    return path


def load_kmer_file(path: Path, label: str) -> dict:
    """Load validated 24-mers in file order without allocating read lists."""
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

            kmers[kmer] = None

    if not kmers:
        raise ValueError(f"{label} contains no valid {KMER_SIZE}-mers: {path}")

    kmers["lengths"] = 0
    kmers["total_reads"] = 0
    return kmers


def load_classifications(path: Path) -> Dict[str, str]:
    """
    Load the k-mer classification JSON.

    Expected format:
        {
            "AAAAAAAA...": "DYZ1",
            "CCCCCCCC...": "DYZ2"
        }
    """
    with path.open("r", encoding="utf-8") as handle:
        try:
            classifications = json.load(handle)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Classification JSON is not valid JSON ({path}): {exc}"
            ) from exc

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


def check_classification_overlap(
    yq12_subunit: dict,
    classifications: Dict[str, str],
    allow_unclassified: bool,
) -> None:

    subunit_total = len(yq12_subunit) - len(COUNTER_KEYS)
    matched = sum(1 for kmer in classifications if kmer in yq12_subunit)

    if matched == 0:
        message = (
            "None of the k-mers in the classification JSON appear in the Yq12 "
            "subunit k-mer file. Every DYZ1/DYZ2 estimate would be zero. Check "
            "that the two files come from the same k-mer set."
        )
        if not allow_unclassified:
            raise ValueError(message + " Use --allow-unclassified to run anyway.")
        log(f"Warning: {message}")
        return

    if matched < subunit_total:
        percent = 100.0 * matched / subunit_total
        log(
            f"Note: {matched} of {subunit_total} subunit k-mers "
            f"({percent:.1f}%) have a classification; the rest will be "
            "counted as unclassified."
        )


def check_kmer_subset(
    child: dict,
    parent: dict,
    child_label: str,
    parent_label: str,
    allow_disjoint: bool,
) -> None:

    child_total = len(child) - len(COUNTER_KEYS)
    missing = sum(
        1 for kmer in child if kmer not in COUNTER_KEYS and kmer not in parent
    )

    if missing == 0:
        return

    if missing == child_total:
        message = (
            f"No k-mer in the {child_label} appears in the {parent_label}. "
            "SRKmer only counts a subunit k-mer when that same k-mer is also "
            f"present in the {parent_label}, so every count would be zero. "
            "Check that these two files belong to the same k-mer set."
        )
        if not allow_disjoint:
            raise ValueError(message + " Use --allow-disjoint-kmers to run anyway.")
        log(f"Warning: {message}")
        return

    percent = 100.0 * missing / child_total
    log(
        f"Warning: {missing} of {child_total} k-mers in the {child_label} "
        f"({percent:.1f}%) are absent from the {parent_label} and can never be "
        "counted."
    )


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


def check_alignment_readable(path: Path, mode: str, open_kwargs: dict) -> None:
    """Open and close the alignment before the slower setup work begins."""
    try:
        with pysam.AlignmentFile(str(path), mode, **open_kwargs):
            pass
    except (ValueError, OSError) as exc:
        hint = ""
        if mode == "rc":
            hint = (
                " If this CRAM was made with a different reference, pass the "
                "matching FASTA with --reference."
            )
        raise ValueError(f"Could not open {path}: {exc}.{hint}") from exc


def sample_name_from_path(path: Path) -> str:
    return path.stem


def prepare_targets(yq12, centromere, yq12_subunit, dyz3):

    subunit_kmers = [k for k in yq12_subunit if k not in COUNTER_KEYS]
    subunit_indices = {k: i for i, k in enumerate(subunit_kmers)}
    # (region index, Yq12-subunit index or None, is DYZ3)
    targets = {
        k: (0, subunit_indices.get(k), False)
        for k in yq12 if k not in COUNTER_KEYS
    }
    for kmer in centromere:
        if kmer not in COUNTER_KEYS and kmer not in targets:
            targets[kmer] = (1, None, kmer in dyz3)
    return targets, subunit_kmers


class DiskReadStore:


    def __init__(self, connection):
        self.connection = connection
        connection.execute("PRAGMA page_size=4096")
        connection.execute("PRAGMA cache_size=-8192")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA mmap_size=0")
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute(
            "CREATE TABLE hits (kmer_index INTEGER NOT NULL, "
            "identifier TEXT NOT NULL, read_length INTEGER NOT NULL)"
        )
        self.cursor = connection.cursor()

    def add(self, kmer_index, identifier, read_length):
        self.cursor.execute(
            "INSERT INTO hits VALUES (?, ?, ?)",
            (kmer_index, identifier, read_length),
        )

    def polished_totals(self):
        self.connection.commit()
        unique_count = 0
        total_bp = 0

        rows = self.connection.execute(
            "SELECT read_length FROM hits GROUP BY kmer_index, identifier"
        )
        for (read_length,) in rows:
            unique_count += 1
            total_bp += read_length
        return unique_count, total_bp

    def write_reads(self, path, subunit_kmers):
        written = 0

        rows = self.connection.execute(
            "SELECT kmer_index, identifier FROM hits ORDER BY kmer_index, rowid"
        )
        with path.open("w", encoding="utf-8") as handle:
            for kmer_index, identifier in rows:
                handle.write(f"{subunit_kmers[kmer_index]}\t{identifier}\n")
                written += 1
        return written


def scan_alignment(samfile, targets, subunit_kmers, classifications, store):

    counts = [{"lengths": 0, "total_reads": 0} for _ in range(4)]
    subunit_counts, dyz3_counts = counts[2:]
    per_class_bp = defaultdict(int)
    unclassified_records = 0
    total_reads = 0
    lookup = targets.get
    add_record = store.add

    for record in samfile.fetch(until_eof=True):
        total_reads += 1
        sequence = record.query_sequence
        qualities = record.query_qualities
        if sequence is None or qualities is None:
            continue
        read_length = len(sequence)
        for start in range(read_length - KMER_SIZE + 1):
            target = lookup(sequence[start:start + KMER_SIZE])
            if target is None:
                continue
            if min(qualities[start:start + KMER_SIZE]) < 20:
                continue
            region, subunit_index, is_dyz3 = target
            counts[region]["total_reads"] += 1
            counts[region]["lengths"] += read_length
            if subunit_index is not None:
                subunit_counts["total_reads"] += 1
                subunit_counts["lengths"] += read_length
                kmer = subunit_kmers[subunit_index]
                classification = classifications.get(kmer, "NO")
                per_class_bp[classification] += read_length
                if classification == "NO":
                    unclassified_records += 1
                # str(None) also preserves the original absent-QNAME behavior.
                identifier = str(record.query_name) + "_" + str(read_length)
                add_record(subunit_index, identifier, read_length)
            elif is_dyz3:
                dyz3_counts["total_reads"] += 1
                dyz3_counts["lengths"] += read_length
            break
    unique_count, polished_total_bp = store.polished_totals()
    return (total_reads, counts, dict(per_class_bp), unclassified_records,
            unique_count, polished_total_bp)


def initialize_scan_worker(targets):
    global _WORKER_TARGETS
    _WORKER_TARGETS = targets
    # The coordinator handles Ctrl-C and closes the bounded worker pool.
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def worker_ready():
    return os.getpid()


@contextmanager
def scan_executor(workers, start_method, targets):

    if workers == 1:
        yield None
        return
    if start_method == "auto":
        start_method = "fork" if sys.platform.startswith("linux") else "spawn"
    context = multiprocessing.get_context(start_method)
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=context,
        initializer=initialize_scan_worker, initargs=(targets,),
    ) as executor:

        ready = [executor.submit(worker_ready) for _ in range(workers)]
        for future in ready:
            future.result()
        yield executor


def alignment_batches(samfile, batch_size):

    batch = []
    payload = 0
    for record in samfile.fetch(until_eof=True):
        sequence = record.query_sequence
        qualities = record.query_qualities
        if sequence is None or qualities is None:
            item = (None, None, None)
            size = 64
        else:
            name = str(record.query_name)
            qualities = bytes(qualities)
            item = (sequence, qualities, name)
            size = len(sequence) + len(qualities) + len(name) + 64
        if batch and (len(batch) >= batch_size or
                      payload + size > MAX_BATCH_PAYLOAD):
            yield batch
            batch = []
            payload = 0
        batch.append(item)
        payload += size
    if batch:
        yield batch


def scan_batch(batch):

    counts = [{"lengths": 0, "total_reads": 0} for _ in range(4)]
    subunit_hits = []
    lookup = _WORKER_TARGETS.get
    for sequence, qualities, name in batch:
        if sequence is None:
            continue
        read_length = len(sequence)
        for start in range(read_length - KMER_SIZE + 1):
            target = lookup(sequence[start:start + KMER_SIZE])
            if target is None:
                continue
            if min(qualities[start:start + KMER_SIZE]) < 20:
                continue
            region, subunit_index, is_dyz3 = target
            counts[region]["total_reads"] += 1
            counts[region]["lengths"] += read_length
            if subunit_index is not None:
                counts[2]["total_reads"] += 1
                counts[2]["lengths"] += read_length
                subunit_hits.append(
                    (subunit_index, name + "_" + str(read_length), read_length)
                )
            elif is_dyz3:
                counts[3]["total_reads"] += 1
                counts[3]["lengths"] += read_length
            break
    return len(batch), counts, subunit_hits


def scan_alignment_parallel(samfile, subunit_kmers, classifications, store,
                            executor, workers, batch_size):

    counts = [{"lengths": 0, "total_reads": 0} for _ in range(4)]
    per_class_bp = defaultdict(int)
    unclassified_records = 0
    total_reads = 0
    pending = deque()

    def merge(result):
        nonlocal total_reads, unclassified_records
        batch_total, batch_counts, subunit_hits = result
        total_reads += batch_total
        for total, partial in zip(counts, batch_counts):
            total["total_reads"] += partial["total_reads"]
            total["lengths"] += partial["lengths"]
        for subunit_index, identifier, read_length in subunit_hits:
            kmer = subunit_kmers[subunit_index]
            classification = classifications.get(kmer, "NO")
            per_class_bp[classification] += read_length
            if classification == "NO":
                unclassified_records += 1
            store.add(subunit_index, identifier, read_length)

    try:
        for batch in alignment_batches(samfile, batch_size):
            if len(pending) >= 2 * workers:
                merge(pending.popleft().result())
            pending.append(executor.submit(scan_batch, batch))
        while pending:
            merge(pending.popleft().result())
    finally:
        # Cancel tasks not yet running if decoding, merging or a worker fails.
        for future in pending:
            future.cancel()
    unique_count, polished_total_bp = store.polished_totals()
    return (total_reads, counts, dict(per_class_bp), unclassified_records,
            unique_count, polished_total_bp)


def calculate_streaming_lengths(sample_name, scan_result, sample_depth,
                                dyz1_unit_length, dyz2_unit_length):
    (_, counts, per_class_bp, unclassified_records,
     unique_count, polished_total_bp) = scan_result
    dyz1_total_bp = int(per_class_bp.get("DYZ1", 0))
    dyz2_total_bp = int(per_class_bp.get("DYZ2", 0))
    # Preserve v1.0's arithmetic order for identical floating-point results.
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
        "Total_Yq12_Subunit_Records": counts[2]["total_reads"],
        "Unique_Yq12_Subunit_Records": unique_count,
        "Unclassified_Yq12_Subunit_Records": unclassified_records,
        "DYZ1_Total_Read_Length_bp": dyz1_total_bp,
        "DYZ2_Total_Read_Length_bp": dyz2_total_bp,
        "Polished_Total_Read_Length_bp": polished_total_bp,
        "Per_Class_Read_Length_bp": per_class_bp,
    }
    return sample_lengths, calculation_details


def warn_on_empty_results(calculation_details: dict) -> None:
    total = calculation_details["Total_Yq12_Subunit_Records"]
    unclassified = calculation_details["Unclassified_Yq12_Subunit_Records"]
    dyz1 = calculation_details["DYZ1_Total_Read_Length_bp"]
    dyz2 = calculation_details["DYZ2_Total_Read_Length_bp"]

    if total == 0:
        log(
            "Warning: no Yq12 subunit reads were found, so all length "
            "estimates are zero. Confirm the input alignment contains Y "
            "chromosome reads and that the k-mer files are the intended ones."
        )
        return

    if dyz1 + dyz2 == 0:
        log(
            f"Warning: all {total} Yq12 subunit records were unclassified, so "
            "the DYZ1 and DYZ2 estimates are zero. Confirm the classification "
            "JSON matches the subunit k-mer file."
        )
    elif unclassified:
        percent = 100.0 * unclassified / total
        log(
            f"Note: {unclassified} of {total} Yq12 subunit records "
            f"({percent:.1f}%) had no classification and are excluded from the "
            "DYZ1/DYZ2 totals."
        )


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


def run(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    start_time = time.time()

    input_path = require_file(args.input, "Input alignment file")
    yq12_path = require_file(args.yq12, "Yq12 k-mer file")
    yq12_subunit_path = require_file(args.yq12_subunit, "Yq12 subunit k-mer file")
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
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(
            f"--out must be a directory, but this path is a file: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    open_kwargs = {"threads": args.io_threads}
    if reference_path is not None:
        open_kwargs["reference_filename"] = str(reference_path)

    check_alignment_readable(input_path, mode, open_kwargs)

    log("Loading k-mer files ...")
    yq12 = load_kmer_file(yq12_path, "Yq12 k-mer file")
    yq12_subunit = load_kmer_file(yq12_subunit_path, "Yq12 subunit k-mer file")
    centromere = load_kmer_file(centromere_path, "Centromere k-mer file")
    dyz3 = load_kmer_file(dyz3_path, "DYZ3 k-mer file")
    classifications = load_classifications(classification_path)

    check_kmer_subset(
        yq12_subunit,
        yq12,
        "Yq12 subunit k-mer file",
        "Yq12 k-mer file",
        args.allow_disjoint_kmers,
    )
    check_kmer_subset(
        dyz3,
        centromere,
        "DYZ3 k-mer file",
        "centromere k-mer file",
        args.allow_disjoint_kmers,
    )
    check_classification_overlap(
        yq12_subunit, classifications, args.allow_unclassified
    )

    targets, subunit_kmers = prepare_targets(yq12, centromere, yq12_subunit, dyz3)

    del yq12, centromere, yq12_subunit, dyz3
    temp_parent = (args.temp_dir.expanduser().resolve()
                   if args.temp_dir is not None else output_dir)
    if not temp_parent.is_dir():
        raise ValueError(f"Temporary directory does not exist: {temp_parent}")

    sample_name = sample_name_from_path(input_path)
    summary_path = output_dir / f"{sample_name}_kmerSampleResults.tsv"
    subunit_reads_path = output_dir / f"{sample_name}_yq12Sub_reads.txt"
    sample_lengths_path = output_dir / f"{sample_name}_SampleLengths.json"
    details_path = output_dir / f"{sample_name}_SRKmerCalculationDetails.json"

    log(f"Streaming {input_path.name}; temporary records in {temp_parent} ...")
    log(f"Scan workers: {args.workers}; decoder thread setting: {args.io_threads}.")
    if args.workers > 1:
        available_text = os.environ.get("SLURM_CPUS_PER_TASK")
        if available_text and available_text.isdigit():
            # A planning heuristic, not an exact OS-thread count.
            planned = args.workers + 1 + max(0, args.io_threads - 1)
            if planned > int(available_text):
                log(f"Warning: workers plus reading/decoding may need about "
                    f"{planned} CPUs; SLURM_CPUS_PER_TASK={available_text}.")
    with scan_executor(args.workers, args.start_method, targets) as executor:
        with tempfile.TemporaryDirectory(prefix=".srkmer-", dir=str(temp_parent)) as tmp:
            with closing(sqlite3.connect(str(Path(tmp) / "reads.sqlite3"))) as connection:
                store = DiskReadStore(connection)
                with pysam.AlignmentFile(str(input_path), mode, **open_kwargs) as samfile:
                    if executor is None:
                        result = scan_alignment(
                            samfile, targets, subunit_kmers, classifications, store
                        )
                    else:
                        result = scan_alignment_parallel(
                            samfile, subunit_kmers, classifications, store,
                            executor, args.workers, args.batch_size,
                        )
                sample_lengths, calculation_details = calculate_streaming_lengths(
                    sample_name, result, args.sample_depth,
                    args.dyz1_unit_length, args.dyz2_unit_length,
                )
                warn_on_empty_results(calculation_details)
                log("Writing results in the original k-mer/record order ...")
                reads_written = store.write_reads(subunit_reads_path, subunit_kmers)


    total_reads, counts = result[:2]
    yq12_reads, cent_reads, yq12_sub_reads, dyz3_reads = counts
    runtime_seconds = time.time() - start_time
    estimates = sample_lengths[sample_name]
    write_summary(
        output_path=summary_path, append=args.append, total_reads=total_reads,
        yq12_reads=yq12_reads, cent_reads=cent_reads,
        yq12_sub_reads=yq12_sub_reads, dyz3_reads=dyz3_reads,
        sample_depth=args.sample_depth, estimates=estimates,
        runtime_seconds=runtime_seconds,
    )
    write_json(sample_lengths_path, sample_lengths)
    write_json(details_path, calculation_details)

    print(f"Summary written to: {summary_path}")
    if reads_written:
        print(f"Yq12 subunit reads written to: {subunit_reads_path}")
    else:
        print(f"No Yq12 subunit reads; wrote an empty {subunit_reads_path}")
    print(f"SampleLengths dictionary written to: {sample_lengths_path}")
    print(f"Calculation details written to: {details_path}")
    print(f"Completed in {runtime_seconds:.2f} seconds.")

    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        return run(args, parser)
    except KeyboardInterrupt:
        parser.exit(status=130, message="Interrupted before results were written.\n")
    except (
        FileNotFoundError,
        ValueError,
        OSError,
        KeyError,
        TypeError,
        pysam.utils.SamtoolsError,
        sqlite3.Error,
        RuntimeError,
    ) as exc:
        parser.exit(status=1, message=f"Error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
