# EXTRA FILES:
The sampleDepths.json file contains average sample read depths across the XDR regions of 1500+ samples from the 1000 genomes. In case anyone is using those files and wanted that info or wanted to know the depths I used per sample in the manuscript. 

# Analysis Notebooks

This directory contains the Jupyter notebooks used to combine Yq12 estimates, benchmark the SRKmer and LRKmer methods, evaluate population-level patterns, and perform father–son statistical analyses.

These notebooks represent the analysis workflow used during manuscript preparation. They contain exploratory code, intermediate checks, figure-generation code, and final statistical analyses. Some file paths are specific to the original analysis environment and would need to be updated before the notebooks could be rerun elsewhere.

## `Part1_AllTests.ipynb`

This notebook contains the primary data integration and statistical analysis workflow for the Yq12 estimation study.

The notebook combines results from multiple sources, including:

* SRKmer Yq12 estimates
* LRKmer Yq12 estimates
* Completely assembled Yq12 regions
* Sample sequencing-depth information
* Y-chromosome haplogroup annotations
* Population and pedigree information
* Father–son relationships from the 1000 Genomes Project and related datasets

### Major analyses

#### SRKmer and LRKmer benchmarking

The notebook compares SRKmer and LRKmer estimates against samples with completely assembled Yq12 regions.

These analyses include:

* Pearson correlations between estimated and assembled Yq12 lengths
* Rank-based comparisons using Spearman-style rank correlations and Kendall’s tau
* Comparisons between SRKmer and LRKmer estimates
* Scatterplots comparing estimated and assembled Yq12 lengths
* Examination of estimation residuals and rank differences

#### Yq12 repeat composition

The notebook calculates and evaluates the estimated numbers of DYZ1 and DYZ2 repeat units.

This includes:

* Calculation of the DYZ1:DYZ2 subunit ratio
* Evaluation of the relationship between DYZ1 and DYZ2 counts
* Visualization of the distribution of DYZ1:DYZ2 ratios
* Comparison of SRKmer-derived repeat composition with assembled Yq12 repeat composition

#### Population-level Yq12 variation

The notebook summarizes the distribution of estimated Yq12 lengths across the population-scale dataset.

Analyses include:

* Descriptive summaries of Yq12 size
* Histograms of estimated Yq12 lengths
* Comparison of Yq12 length across Y-chromosome haplogroups
* Boxplots and individual-sample plots organized by broad and detailed haplogroups
* Kruskal–Wallis tests for differences among haplogroups
* Dunn post-hoc tests with multiple-testing correction
* Hodges–Lehmann pairwise differences with bootstrap confidence intervals
* Additional exploratory effect-size and regression analyses

#### Father–son Yq12 length analysis

The notebook identifies father–son pairs and constructs paired datasets containing:

* Estimated Yq12 length
* DYZ1:DYZ2 ratio
* Y-chromosome haplogroup
* Population assignment
* Percent difference between fathers and sons

The primary father–son analysis excludes pairs with an absolute estimated Yq12 length difference of 20% or greater, while retaining the complete dataset for sensitivity analyses.

Statistical analyses include:

* Wilcoxon signed-rank tests
* Exact binomial sign tests
* Hodges–Lehmann paired pseudomedian differences
* Bootstrap confidence intervals
* Comparison of paired differences among haplogroups
* Fisher’s exact tests relating the direction of father–son differences to the father’s position relative to global or haplogroup-specific Yq12 medians

The notebook also generates paired father–son plots, percent-difference distributions, haplogroup-stratified figures, and intermediate CSV files used by subsequent analyses.

### Intermediate outputs

Among its intermediate outputs, this notebook creates the father–son datasets used by the Part 2 notebook:

```text
fatherSon_Yq12.csv
fatherSon_RatioDF_Yq12.csv
```

---

## `Part2_More_FatherVsSonStatisticalTesting.ipynb`

This notebook performs a focused follow-up analysis of father–son DYZ1:DYZ2 subunit ratios.

It uses the paired father–son datasets generated in Part 3 and asks whether sons tend to have DYZ1:DYZ2 ratios that are closer to an equal ratio of 1 than those of their fathers.

### Main analysis

Father and son ratios are log-transformed so that:

```text
log(ratio) = 0
```

corresponds to a DYZ1:DYZ2 ratio of 1.

For each pair, the notebook calculates the difference between the son’s and father’s absolute log-distance from 1. A negative value indicates that the son’s ratio is closer to 1 than the father’s ratio.

The statistical analyses include:

* The proportion of sons whose ratios are closer to 1
* An exact binomial confidence interval for this proportion
* A one-sided exact binomial sign test
* A one-sided Wilcoxon signed-rank test of the paired distance differences
* A regression-through-the-origin analysis of log-transformed father and son ratios
* A test of whether the father-to-son slope is less than 1
* The proportion of pairs remaining on the same side of a ratio of 1

### Visualizations

The notebook generates paired dumbbell plots showing father and son ratios using:

* Raw DYZ1:DYZ2 ratios
* Log-transformed ratios
* Absolute log-distance from a ratio of 1

These figures illustrate whether each son’s ratio moved toward or away from an equal DYZ1:DYZ2 ratio relative to his father.

---

## Relationship between the notebooks

`Part3_Publication-TheBigKahuna_AllTests.ipynb` is the main analysis notebook. It integrates the datasets, benchmarks the estimators, evaluates population and haplogroup patterns, and performs the primary father–son Yq12 length analyses.

`Part4_FatherVsSonStatisticalTesting.ipynb` is a narrower follow-up notebook focused specifically on intergenerational changes in the DYZ1:DYZ2 subunit ratio.

The notebooks are provided primarily for transparency and reproducibility of the analyses reported in the associated study.
