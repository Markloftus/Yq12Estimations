import gzip
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
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
        "-t", "--threads", type=int, default=8,
        help="Number of worker threads (default: 8).",
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
    global K_VALUE, K_MASK, KNOWN_D, MSY_D, CMAP_D, KNOWN_SET_PY, MSY_SET_PY

    known_arr = np.load(os.path.join(cache_dir, f"known_k{k}.npy"), mmap_mode="r")
    msy_arr   = np.load(os.path.join(cache_dir, f"msy_k{k}.npy"),   mmap_mode="r")
    cmap_arr  = np.load(os.path.join(cache_dir, f"cmap_k{k}.npy"),  mmap_mode="r")

    K_VALUE = int(k)
    K_MASK  = (1 << (2*k)) - 1

    # Python sets for quick membership in the rare second pass
    KNOWN_SET_PY = set(np.ndarray.tolist(known_arr))
    MSY_SET_PY   = set(np.ndarray.tolist(msy_arr))

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
        return header, 0.0, n, {}, 'NONE'

    bs = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    hit, bad, msy_keys, msy_cnts = _scan_windows(bs, K_VALUE, K_MASK, _B2I_ARR, KNOWN_D, MSY_D, CMAP_D)

    allCount = n - K_VALUE + 1
    proportion = hit / allCount if allCount else 0.0

    # build MSY counts dict (int k-mers)
    local_msy_int = {int(msy_keys[i]): int(msy_cnts[i]) for i in range(msy_keys.size)}

    # Only materialize badmers when needed
    if proportion >= 0.90 and bad > 0:
        badmers = []
        code = 0; ok = 0; mask = K_MASK
        for i, ch in enumerate(bs):
            b = int(_B2I_ARR[ch])   # <— cast to Python int to avoid int8 overflow
            if b < 0:
                ok = 0
            else:
                code = ((code << 2) & mask) | b
                if ok < K_VALUE:
                    ok += 1
            if i >= K_VALUE - 1:
                start = i - (K_VALUE - 1)
                if ok < K_VALUE:
                    badmers.append(seq[start:start+K_VALUE]); continue
                # membership checks...
                if (code in KNOWN_D) or (code in MSY_D):
                    continue
                badmers.append(seq[start:start+K_VALUE])
        return header, proportion, n, local_msy_int, badmers

    return header, proportion, n, local_msy_int, 'NONE'


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

def read_fastq_multithreaded(path, num_threads=16):
    readDict = {}
    total_msy_int = Counter()
    _ = _scan_windows(np.frombuffer(b"ACGTACGTACGTACGT", dtype=np.uint8),
                      max(4, K_VALUE), (1 << (2*max(4, K_VALUE))) - 1,
                      _B2I_ARR, KNOWN_D, MSY_D, CMAP_D)

    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        for header, prop, size, local_msy_int, badmers in ex.map(lambda r: process_read(*r), stream_fastq(path)):
            readDict[header] = {"proportion": prop, "size": size, "badmers": badmers}
            if local_msy_int:
                total_msy_int.update(local_msy_int)

    return readDict, total_msy_int  # NOTE: int-coded k-mers


if __name__ == "__main__":
    args = parse_args()

    if args.kmer_size <= 0:
        raise ValueError("--kmer-size must be greater than zero.")
    if args.threads <= 0:
        raise ValueError("--threads must be greater than zero.")

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

    MSY_lookup = set()
    MSY_canon_map = {}
    for kmer in list(MSYuniqueKmers.keys()):
        rc = _revcomp_str(kmer)
        canon = kmer if kmer <= rc else rc
        MSY_lookup.add(kmer)
        MSY_canon_map[kmer] = canon
        MSY_lookup.add(rc)
        MSY_canon_map[rc] = canon

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

        cache_dir.mkdir(parents=True, exist_ok=True)
        build_kmer_cache(k, known_set, MSY_lookup, MSY_canon_map, str(cache_dir))

    load_kmer_cache(k, str(cache_dir))

    readDict, msy_counts_int = read_fastq_multithreaded(
        str(fastq_path), num_threads=num_threads
    )

    # Fold observed MSY k-mer counts back into the normalization dictionary.
    for code_int, count in msy_counts_int.items():
        kmer = _dec_kmer(code_int, K_VALUE)
        MSYuniqueKmers[kmer] = MSYuniqueKmers.get(kmer, 0) + count

    df2 = pd.DataFrame.from_dict(readDict, orient="index")
    depth_vals = [value for value in MSYuniqueKmers.values() if value != 0]
    depth = np.median(depth_vals) if depth_vals else 0.0
    df2["Depth"] = depth
    df2["normalizedSize"] = (
        df2["size"] / depth if depth != 0 else np.nan
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df2.to_csv(output_path)

