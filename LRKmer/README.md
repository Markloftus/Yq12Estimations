# STEP 1: Unzip the Cache Directory zip file
You will need these files for the LRKmer to work.</br>
What files need to be/can be found in the cache-dir:

- known_k24.npy </br>
- msy_k24.npy </br>
- cmap_k24.npy </br>
- XDR_Kmers_filled_Filtered.json (through this in the zip but you can place it anywhere just point to it obviously)


# STEP 2: Run the LRKmer
How to run:</br>

python LRKmer-v1.0.py \
    --input sample.fastq.gz \
    --out sample_LRKmer.csv \
    --msy-kmers /YourDirectory/path/XDR_Kmers_filled_Filtered.json \
    --cache-dir /YourDirectory/path/ \
    --threads 8



