STEP 1: Unzip the Cache Directory zip file
# You will need these files for the LRKmer to work.
# What files need to be in the cache-dir:

known_k24.npy
msy_k24.npy
cmap_k24.npy


STEP2: Run the LRKmer
# How to run:

python LRKmer-v1.0.py \
    --input sample.fastq.gz \
    --out sample_LRKmer.csv \
    --msy-kmers /YourDirectory/path/XDR_Kmers_filled_Filtered.json \
    --cache-dir /YourDirectory/path/ \
    --threads 8



