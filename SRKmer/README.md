# STEP 1:
-NOTE: You will need to calculate the average XDR (Y-chromosome X degenerate region) sample depth before running the SRKmer as it takes in this depth as a parameter (I am using 15 as an example below; GRCh38.win.bed is supplied in supplemental tables but I will supply it again in this directory).

1) samtools bedcov -Q 20 GRCh38.win.bed CRAMFILE >> $SAMPLENAME.Cov.out
2) Then I read in the cov.out in python:
    - keyDF = pd.read_csv(directory+sample, sep='\t', header=None,on_bad_lines='skip')
    Then I take the mean:
    - float(np.mean(keyDF[3])/1000)
    This number gets passed as the --sample-depth float/integer below


# STEP 2: Compile the Cython extension (ONLY IF USING V1.0. The newer versions do NOT need this anymore)

1) python setup.py build_ext --inplace

This produces a compiled module with a name similar to:
    -mainKmer.cpython-311-x86_64-linux-gnu.so
The exact filename will depend on the operating system and Python version.

# STEP 3: Run the SRKmer

How to run SRKmer: </br>
-(note: I suppled kmerENV.yml this is the conda environment I utilized during the runs on the compute cluster in case you just want to remake it - not necessary). <br><br>

python SRKmer-v1.0.py \
    --input HG00477.cram \
    --reference GRCh38.primary_assembly.genome.fa \
    --yq12 yq12_CombinedUnion_1.0.txt \
    --yq12s Yq12_Union_1.0.txt \
    --cent centromere_CombinedUnion_1.0.txt \
    --cents Centromere_Union_1.0.txt \
    --classification-json my_95Percent_Identity_SubunitKmers.json \
    --sample-depth 15 \
    --out results/

If you are using the updated SRKmer-v1_1.py (or eventually newer) the same commands need to be utilized/files be fed to the SRKmer. Though you no longer need to do step 2. You can also specify: <br><br>
--temp-dir /Your/Path/srkmer_tmp/ <br><br>
--workers 7 <br><br>
This (--temp-dir) will tell the SRKmer where to write temporary files and --workers is how many different threads you want reading the file at a time. If you do not pass it --temp-dir it will build the temporary files in --out. You should set the number of workers to one less than the total CPUs you ask for. So if you are running a slurm job (#SBATCH --cpus-per-task=8) you should specify 7 workers. This gives seven scanning workers plus one CPU for reading and merging. For each additional worker you should provide it about 0.5Gb of memory. <br>
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
