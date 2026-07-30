# STEP 1:
-NOTE: You will need to calculate the average XDR sample depth before running the SRKmer as it takes in the depth as a paramter (I am using 15 as an example below; GRCh38.win.bed is supplied in supplemental tables but I will supply it again in this directory).

1) samtools bedcov -Q 20 GRCh38.win.bed CRAMFILE >> $SAMPLENAME.Cov.out
2) Then I read in the cov.out in python:
    - keyDF = pd.read_csv(directory+sample, sep='\t', header=None,on_bad_lines='skip')
    Then I take the mean:
    - float(np.mean(keyDF[3])/1000)
    This number gets passed as the --sample-depth float/integer below


# STEP 2: Compile the Cython extension

1) python setup.py build_ext --inplace

This produces a compiled module with a name similar to:
    -mainKmer.cpython-311-x86_64-linux-gnu.so
The exact filename will depend on the operating system and Python version.

# STEP 3: Run the SRKmer

How to run SRKmer: </br>
-(note: I suppled kmerENV.yml this is the conda environment I utilized during the runs on the compute cluster in case you just want to remake it - not necessary)

python SRKmer-v1.0.py \
    --input HG00477.cram \
    --reference GRCh38.primary_assembly.genome.fa \
    --yq12 Yq12_kmers.txt \
    --yq12s Yq12_subunit_kmers.txt \
    --cent Y_centromere_kmers.txt \
    --cents DYZ3_kmers.txt \
    --classification-json my_95Percent_Identity_SubunitKmers.json \
    --sample-depth 15 \
    --out results/

## The program will write these files:
results/HG00477_kmerSampleResults.tsv </br>
results/HG00477_yq12Sub_reads.txt </br>
results/HG00477_SampleLengths.json </br>
results/HG00477_SRKmerCalculationDetails.json </br>


You will want to utilize the _sampleLengths.json information (I use the 'Unpolished' also named 'Refined' in the manuscript text. The Polished tends to over estimate but that is useful as a max length if you want a min-max window.):

{
    "HG00477": {
        "DYZ1_Counts": 5554.010268,
        "DYZ2_Counts": 5725.581489,
        "Polished": 39131890.0,
        "Unpolished": 33678170.0
    }
}
