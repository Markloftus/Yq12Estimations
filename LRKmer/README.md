## Current version: v1.1

# STEP 1: Unzip the Cache Directory zip file
You will need these files for the LRKmer to work.</br>
What files need to be/can be found in the cache-dir:

- known_k24.npy </br>
- msy_k24.npy </br>
- cmap_k24.npy </br>
- XDR_Kmers_filled_Filtered.json 


# STEP 2: Run the LRKmer
How to run:</br>

python LRKmer-v1.1.py \
    --input sample.fastq.gz \
    --out sample_LRKmer.csv \
    --msy-kmers /YourDirectory/path/XDR_Kmers_filled_Filtered.json \
    --cache-dir /YourDirectory/path/cache/ \
    --threads 8 \
    --batch-size 512 \
    --temp-dir /YourDirectory/tmp/

| Argument | Description |
|---|---|
| `--input` | Input gzip-compressed FASTQ file. |
| `--out` | Path to the output CSV file. |
| `--msy-kmers` | Path to the MSY normalization k-mer JSON file. |
| `--cache-dir` | Directory containing the precomputed k-mer lookup arrays (`known_k24.npy`, `msy_k24.npy`, and `cmap_k24.npy`). |
| `--threads` | Number of worker threads (default: 8). |
| `--batch-size` | Number of reads processed per worker task (default: 128; recommended: 512). |
| `--temp-dir` | Directory for temporary SQLite files (optional). Temporary files are automatically removed after processing. |

# STEP 3: Analyze
After running the LRKmer you will receive a CSV file with all the read Yq12 kmer proportions. I suggest filtering by keeping all reads >=.50 in the 'proportion' column then summing the 'normalizedSize' column.
<br>
## Log:
November 19th, 2026 v1.1 was created. This uses around 75% less memory and is much faster. Added a temp folder parameter as well.
