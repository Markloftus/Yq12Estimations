import gzip
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
import sqlite3
import tempfile
from collections import Counter
import pandas as pd
import os, json, numpy as np
from pathlib import Path
from numba.typed import Dict as NDict
from numba import njit, types
from numba.typed import Dict
import argparse as ap


def parse_args():
    parser = ap.ArgumentParser(
        description=(
            "Estimate Yq12 read content from a gzip-compressed FASTQ using "
            "precomputed 24-mer lookup tables."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input FASTQ file compressed with gzip (.fastq.gz or .fq.gz).",
    )
    parser.add_argument(
        "-o", "--out", dest="output", required=True,
        help="Path for the output CSV file.",
    )
    parser.add_argument(
        "--msy-kmers", required=True,
        help=(
            "Path to XDR_Kmers_filled_Filtered.json. This file supplies "
            "the MSY normalization k-mers and their initial counts."
        ),
    )
    parser.add_argument(
        "--cache-dir", required=True,
        help=(
            "Directory containing known_k{K}.npy, msy_k{K}.npy, and "
            "cmap_k{K}.npy."
        ),
    )
    parser.add_argument(
        "--known-kmers", default=None,
        help=(
            "Path to Yq12_24Mers_filled.json. Required only when the cache "
            "files are missing or --rebuild-cache is used."
        ),
    )
    parser.add_argument(
        "--rebuild-cache", action="store_true",
        help="Rebuild the .npy lookup tables from the JSON k-mer files.",
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=24,
        help="K-mer length used by the lookup files (default: 24).",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4,
        help="Number of worker threads (default: 4).",
    )
    parser.add_argument(
        "--temp-dir", default=None,
        help="Parent directory for an automatically cleaned temporary SQLite database "
             "(default: system temporary directory). Use fast local storage if available.",
    )
    parser.add_argument(
        "--batch-size", type=int, default=512,
        help="Reads per worker task (default: 512). Smaller uses less memory; "
             "larger may reduce task scheduling overhead.",
    )
    return parser.parse_args()


def require_file(path_value, label):
    path = Path(path_value).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} was not found: {path}")
    return path


def load_json(path):
    with path.open("r") as handle:
        return json.load(handle)

# 2-bit map (Python-side)
_B2I_py = {ord('A'):0, ord('C'):1, ord('G'):2, ord('T'):3,
           ord('a'):0, ord('c'):1, ord('g'):2, ord('t'):3}

def _enc_kmer_skipN(s: str) -> int | None:
    v = 0
    for ch in s.encode('ascii'):
        b = _B2I_py.get(ch, -1)
        if b < 0:   # N or invalid
            return None
        v = (v << 2) | b
    return v

def build_kmer_cache(k: int,
                     known_set: set[str],
                     MSY_lookup: set[str],
                     MSY_canon_map: dict[str, str],
                     cache_dir: str):
    """
    Encode once and save to disk:
      - known_k{k}.npy       : int64 array of known k-mers
      - msy_k{k}.npy         : int64 array of MSY k-mers (both strands if provided)
      - cmap_k{k}.npy        : int64 Nx2 array of (kmer_int, canon_int)
      - meta_k{k}.json       : small metadata (k, sizes)
    Skips kmers containing 'N'.
    """
    os.makedirs(cache_dir, exist_ok=True)

    known_int = []
    for s in known_set:
        v = _enc_kmer_skipN(s)
        if v is not None:
            known_int.append(v)

    msy_int = []
    for s in MSY_lookup:
        v = _enc_kmer_skipN(s)
        if v is not None:
            msy_int.append(v)

    pairs = []
    for km, canon in MSY_canon_map.items():
        vk = _enc_kmer_skipN(km)
        vc = _enc_kmer_skipN(canon)
        if vk is not None and vc is not None:
            pairs.append((vk, vc))

    # Deduplicate & sort 
    known_arr = np.unique(np.asarray(known_int, dtype=np.int64))
    msy_arr   = np.unique(np.asarray(msy_int,   dtype=np.int64))
    if pairs:
        cmap_arr = np.unique(np.asarray(pairs, dtype=np.int64), axis=0)
    else:
        cmap_arr = np.zeros((0,2), dtype=np.int64)

    np.save(os.path.join(cache_dir, f"known_k{k}.npy"), known_arr)
    np.save(os.path.join(cache_dir, f"msy_k{k}.npy"),   msy_arr)
    np.save(os.path.join(cache_dir, f"cmap_k{k}.npy"),  cmap_arr)

    meta = {
        "k": int(k),
        "known_size": int(known_arr.size),
        "msy_size": int(msy_arr.size),
        "cmap_size": int(cmap_arr.shape[0]),
    }
    with open(os.path.join(cache_dir, f"meta_k{k}.json"), "w") as fh:
        json.dump(meta, fh)



def load_kmer_cache(k: int, cache_dir: str):
    """Load .npy tables and populate global lookups ."""
    global K_VALUE, K_MASK, KNOWN_D, MSY_D, CMAP_D

    known_arr = np.load(os.path.join(cache_dir, f"known_k{k}.npy"), mmap_mode="r")
    msy_arr   = np.load(os.path.join(cache_dir, f"msy_k{k}.npy"),   mmap_mode="r")
    cmap_arr  = np.load(os.path.join(cache_dir, f"cmap_k{k}.npy"),  mmap_mode="r")

    K_VALUE = int(k)
    K_MASK  = (1 << (2*k)) - 1

    # Numba typed dicts for the hot JIT path
    KNOWN_D = NDict.empty(key_type=types.int64, value_type=types.boolean)
    for v in known_arr:
        KNOWN_D[int(v)] = True

    MSY_D = NDict.empty(key_type=types.int64, value_type=types.boolean)
    for v in msy_arr:
        MSY_D[int(v)] = True

    CMAP_D = NDict.empty(key_type=types.int64, value_type=types.int64)
    # Expect shape (N,2); safe if empty
    for i in range(cmap_arr.shape[0]):
        CMAP_D[int(cmap_arr[i,0])] = int(cmap_arr[i,1])

def cache_exists(cache_dir, k):
    paths = [f"known_k{k}.npy", f"msy_k{k}.npy", f"cmap_k{k}.npy"]
    return all(os.path.exists(os.path.join(cache_dir, p)) for p in paths)


# Reverse complement 
_trans = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def _revcomp_str(s: str) -> str:
    return s.translate(_trans)[::-1]

# 2-bit helpers
_INT2B = "ACGT"
def _dec_kmer(code: int, k: int) -> str:
    out=[]; v=code
    for _ in range(k):
        out.append(_INT2B[v & 3]); v >>= 2
    return "".join(reversed(out))

# Strict encoder that SKIPS kmers containing non-ACGT
_B2I_py = {ord('A'):0, ord('C'):1, ord('G'):2, ord('T'):3,
           ord('a'):0, ord('c'):1, ord('g'):2, ord('t'):3}

def _enc_kmer_skipN(s: str):
    v = 0
    for ch in s.encode('ascii'):
        b = _B2I_py.get(ch, -1)
        if b < 0:
            return None
        v = (v << 2) | b
    return v

# Fast byte->base map as a NumPy array 
_B2I_ARR = np.full(256, -1, dtype=np.int8)
for c,v in [(b'A',0),(b'C',1),(b'G',2),(b'T',3),
            (b'a',0),(b'c',1),(b'g',2),(b't',3)]:
    _B2I_ARR[c[0]] = v

def prepare_kmer_tables(k, known_set, MSY_lookup, MSY_canon_map):
    """Build Numba-typed dicts & masks. Call ONCE before threading."""
    global K_VALUE, K_MASK, KNOWN_D, MSY_D, CMAP_D
    K_VALUE = int(k)
    K_MASK  = (1 << (2*k)) - 1

    # Encode to ints (skip kmers with N)
    known_int = [v for s in known_set if (v := _enc_kmer_skipN(s)) is not None]
    msy_int   = [v for s in MSY_lookup if (v := _enc_kmer_skipN(s)) is not None]

    KNOWN_D = Dict.empty(key_type=types.int64, value_type=types.boolean)
    for v in known_int:
        KNOWN_D[v] = True

    MSY_D = Dict.empty(key_type=types.int64, value_type=types.boolean)
    for v in msy_int:
        MSY_D[v] = True

    CMAP_D = Dict.empty(key_type=types.int64, value_type=types.int64)
    for kmer, canon in MSY_canon_map.items():
        vk = _enc_kmer_skipN(kmer)
        vc = _enc_kmer_skipN(canon)
        if vk is not None and vc is not None:
            CMAP_D[vk] = vc

@njit(cache=True, nogil=True, fastmath=True)
def _scan_windows(bs, k, mask, b2i, KNOWN_D, MSY_D, CMAP_D):
    n = bs.size
    if n < k:
        return 0, 0, np.empty(0, np.int64), np.empty(0, np.int32)

    code = 0
    ok = 0
    hit = 0
    bad = 0

    # small dynamic table for per-read MSY counts
    cap = 32
    keys = np.empty(cap, np.int64)
    cnts = np.zeros(cap, np.int32)
    sz = 0

    for i in range(n):
        b = b2i[bs[i]]
        if b < 0:
            ok = 0
        else:
            code = ((code << 2) & mask) | b
            if ok < k:
                ok += 1

        if i >= k - 1:
            if ok < k:
                bad += 1
                continue

            # known hit
            if KNOWN_D.get(code, False):
                hit += 1
                continue

            # MSY hit
            if MSY_D.get(code, False):
                canon = CMAP_D.get(code, code)

                # linear probe in tiny per-read table
                found = -1
                for j in range(sz):
                    if keys[j] == canon:
                        found = j
                        break

                if found >= 0:
                    cnts[found] += 1
                else:
                    if sz == cap:
                        # grow by 2x
                        cap2 = cap * 2
                        k2 = np.empty(cap2, np.int64)
                        c2 = np.zeros(cap2, np.int32)
                        for t in range(cap):
                            k2[t] = keys[t]
                            c2[t] = cnts[t]
                        keys = k2
                        cnts = c2
                        cap = cap2
                    keys[sz] = canon
                    cnts[sz] = 1
                    sz += 1
            else:
                bad += 1

    return hit, bad, keys[:sz].copy(), cnts[:sz].copy()

def process_read(header, seq):
    if KNOWN_D is None or MSY_D is None or CMAP_D is None:
        raise RuntimeError("Call prepare_kmer_tables(...) before process_read().")

    n = len(seq)
    if n < K_VALUE:
        return header, 0.0, n, {}

    bs = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    hit, _bad, msy_keys, msy_cnts = _scan_windows(bs, K_VALUE, K_MASK, _B2I_ARR, KNOWN_D, MSY_D, CMAP_D)

    allCount = n - K_VALUE + 1
    proportion = hit / allCount if allCount else 0.0

    # build MSY counts dict (int k-mers)
    local_msy_int = {int(msy_keys[i]): int(msy_cnts[i]) for i in range(msy_keys.size)}

    return header, proportion, n, local_msy_int


def stream_fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline()
            fh.readline()  # plus
            fh.readline()  # qual
            yield h.strip(), s.strip()

def _process_read_batch(reads):
    """Batch Python -> Numba calls, retaining the original per-read algorithm."""
    return [process_read(*read) for read in reads]


def _make_temporary_database(path):
    """Store one final result per read header (same duplicate-header behavior as dict)."""
    conn = sqlite3.connect(str(path))
    # Safe for disposable scratch data; avoids journal overhead.
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA cache_size=-8192")  # Approx. 8 MiB SQLite cache
    conn.execute("""
        CREATE TABLE reads (
            header TEXT PRIMARY KEY,
            proportion REAL NOT NULL,
            size INTEGER NOT NULL
        )
    """)
    return conn


def read_fastq_multithreaded(path, conn, num_threads=16, batch_size=128):
    """Bound outstanding reads, stream results to SQLite, accumulate MSY counts."""
    total_msy_int = Counter()
    _ = _scan_windows(np.frombuffer(b"ACGTACGTACGTACGT", dtype=np.uint8),
                      max(4, K_VALUE), (1 << (2*max(4, K_VALUE))) - 1,
                      _B2I_ARR, KNOWN_D, MSY_D, CMAP_D)

    sql = """INSERT INTO reads (header, proportion, size)
             VALUES (?, ?, ?)
             ON CONFLICT(header) DO UPDATE SET
                proportion=excluded.proportion,
                size=excluded.size"""
    # Like the original dict, the first occurrence sets output row order,
    # while a later duplicate header replaces its values. All occurrences
    # still contribute MSY k-mer counts, just as before.
    pending = deque()
    outstanding_limit = max(2, num_threads * 2)
    uncommitted = 0

    def consume(future):
        nonlocal uncommitted
        results = future.result()
        conn.executemany(sql, (
            (header, prop, size)
            for header, prop, size, _ in results
        ))
        for _, _, _, local_msy_int in results:
            if local_msy_int:
                total_msy_int.update(local_msy_int)
        uncommitted += len(results)
        if uncommitted >= 4096:
            conn.commit()
            uncommitted = 0

    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        stream = stream_fastq(path)
        while True:
            reads = list(islice(stream, batch_size))
            if not reads:
                break
            pending.append(ex.submit(_process_read_batch, reads))
            if len(pending) >= outstanding_limit:
                consume(pending.popleft())
        while pending:
            consume(pending.popleft())

    conn.commit()
    return total_msy_int


def write_output_csv(conn, output_path, depth, rows_per_chunk=10000):
    """Serialize the same pandas columns and read order without a giant frame."""
    cur = conn.execute(
        "SELECT header, proportion, size FROM reads ORDER BY rowid"
    )
    first = True
    while True:
        rows = cur.fetchmany(rows_per_chunk)
        if not rows:
            break
        df2 = pd.DataFrame.from_records(
            rows, columns=["header", "proportion", "size"]
        )
        df2.set_index("header", inplace=True)
        df2.index.name = None  # replicate DataFrame.from_dict(..., orient='index')
        df2["Depth"] = depth
        df2["normalizedSize"] = (
            df2["size"] / depth if depth != 0 else np.nan
        )
        df2.to_csv(output_path, mode="w" if first else "a", header=first)
        first = False

    if first:
        # No input reads: preserve a valid output schema in this corner case.
        pd.DataFrame(columns=[
            "proportion", "size", "Depth", "normalizedSize"
        ]).to_csv(output_path)


if __name__ == "__main__":
    args = parse_args()

    if args.kmer_size <= 0:
        raise ValueError("--kmer-size must be greater than zero.")
    if args.threads <= 0:
        raise ValueError("--threads must be greater than zero.")
    if args.batch_size <= 0:
        raise ValueError("--batch-size must be greater than zero.")

    k = args.kmer_size
    num_threads = args.threads
    fastq_path = require_file(args.input, "Input FASTQ")
    msy_kmers_path = require_file(args.msy_kmers, "MSY k-mer JSON")
    cache_dir = Path(args.cache_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    MSYuniqueKmers = load_json(msy_kmers_path)
    if not isinstance(MSYuniqueKmers, dict):
        raise TypeError(
            "--msy-kmers must contain a JSON object mapping k-mers to counts."
        )

    cache_available = cache_exists(str(cache_dir), k)
    if args.rebuild_cache or not cache_available:
        if args.known_kmers is None:
            missing = ", ".join(
                str(cache_dir / name)
                for name in (
                    f"known_k{k}.npy",
                    f"msy_k{k}.npy",
                    f"cmap_k{k}.npy",
                )
                if not (cache_dir / name).is_file()
            )
            raise FileNotFoundError(
                "The required k-mer cache is incomplete and --known-kmers "
                "was not supplied. Missing cache file(s): " + missing
            )

        known_kmers_path = require_file(
            args.known_kmers, "Yq12 known-kmer JSON"
        )
        known_kmers = load_json(known_kmers_path)
        if isinstance(known_kmers, dict):
            known_set = set(known_kmers.keys())
        elif isinstance(known_kmers, (list, set)):
            known_set = set(known_kmers)
        else:
            raise TypeError(
                "--known-kmers must contain either a JSON object or a JSON list."
            )

        # The strand maps are only needed when generating lookup cache files.
        MSY_lookup = set()
        MSY_canon_map = {}
        for kmer in MSYuniqueKmers:
            rc = _revcomp_str(kmer)
            canon = kmer if kmer <= rc else rc
            MSY_lookup.add(kmer)
            MSY_canon_map[kmer] = canon
            MSY_lookup.add(rc)
            MSY_canon_map[rc] = canon

        cache_dir.mkdir(parents=True, exist_ok=True)
        build_kmer_cache(k, known_set, MSY_lookup, MSY_canon_map, str(cache_dir))
        # Drop large cache-build-only JSON objects before creating Numba dicts.
        del known_kmers, known_set, MSY_lookup, MSY_canon_map

    load_kmer_cache(k, str(cache_dir))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_parent = None
    if args.temp_dir is not None:
        temp_parent = Path(args.temp_dir).expanduser().resolve()
        temp_parent.mkdir(parents=True, exist_ok=True)

    # Depth depends on all reads, so disk-backed per-read results enable
    # one FASTQ pass without keeping all reads in memory.
    with tempfile.TemporaryDirectory(prefix="lrkmer_", dir=temp_parent) as scratch:
        conn = _make_temporary_database(Path(scratch) / "reads.sqlite3")
        try:
            msy_counts_int = read_fastq_multithreaded(
                str(fastq_path), conn, num_threads=num_threads,
                batch_size=args.batch_size,
            )

            # Identical normalization calculation to v1.0.
            for code_int, count in msy_counts_int.items():
                kmer = _dec_kmer(code_int, K_VALUE)
                MSYuniqueKmers[kmer] = MSYuniqueKmers.get(kmer, 0) + count

            depth_vals = [value for value in MSYuniqueKmers.values() if value != 0]
            depth = np.median(depth_vals) if depth_vals else 0.0
            write_output_csv(conn, output_path, depth)
        finally:
            conn.close()
