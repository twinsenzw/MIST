# MIST
Bayesian reconstruction of strain haplotypes from shotgun metagenomics data

MIST (Metagenomic Inference of Strain Types) reconstructs strain haplotypes from shotgun metagenomic data using Bayesian inference by Markov chain Monte Carlo (MCMC). It is derived from the software ["LINEAGE"](http://www.genetics.org/content/197/3/925), with the following modifications:

* MIST implements Metropolis coupled Markov chain Monte Carlo [(MC)3], which runs 4 MCMC chains (1 main/cold chain and 3 heated chains) simultaneously on different threads. Briefly, (MC)3 prevents the MCMC chain from getting stuck in local optima. Please see [Altekar et al. 2004](https://academic.oup.com/bioinformatics/article/20/3/407/186341) for a detailed explanation of (MC)3.

* MIST has improved sampling strategies (see the figure below) that accelerate the convergence of MCMC runs. Briefly, at the beginning of the MCMC run (step < 0.2 x total steps), MIST tends to introduce big alterations into the state space, including changing all strain relative abundances in a sample, changing all SNP alleles within a range of a strain haplotype, and re-assignment of alleles across strains at a SNP locus. Then during the later steps of the MCMC run, MIST tends to introduce fine adjustments to the state space, including bitflipping one SNP locus in one strain (i.e. changing one locus from the reference base to the alternative base, or vice versa) and redistributing strain frequencies between two strains in a single sample. 

* MIST only reconstructs the haplotypes but does not infer the phylogeny. This is because 1) co-estimation of haplotypes and phylogeny significantly increases the parameter space and therefore runtime, 2) the constant-population-size coalescent prior of the phylogeny, as implemented in LINEAGE, are based on assumptions (i.e. constant population size) that are often unreasonable in a microbiota, due to its dynamic temporal fluctuations. 

## To run

MIST is implemented in c++. It can be compiled without anything other than the standard library:
```
g++ -std=c++11 -pthread MIST.cpp -o MIST
```
MIST takes four arguments, two of which are input files:
```
./MIST {Number of strain haplotypes to infer} {Number of MCMC steps (MCMC chain length)} {input file for counts of the "reference" base} {input file for counts of the "alternative" base}
```
Examples of the input files are referenceCount_input.txt and alternativeCount_input.txt.
In these files, each row represent a SNP site, and each column represent a sample in the dataset. For example, according to the example files, there are 2 reference alleles and 0 alternative alleles in the first SNP site in the first sample.
 
Definitions of the reference base and alternative base are the same as LINEAGE: for each biallelic SNP loci, one allele is arbitrarily denoted "reference" and the other "alternative". Because both LINEAGE and MIST models nucleotide substitution using the biallelic Jukes-Cantor model, the actual nucleotide of the reference and alternative alleles does not influence the inference (e.g. reference=A, alternative=T is the same as reference=G, alternative=C). Therefore, MIST denotes reference alleles by "T" and alternative alleles by "A".

MIST output multiple files, most of them are for development and trouble-shooting. The relevant output file is "log_cold.txt". Each line in log_cold.txt represent one step (sample of state) in the main (cold) MCMC chain. 
* The first tab-separated column in each line represents the step number. 
* The second column represents the log posterior of the present parameter state. 
* The third column represents the haplotype of the estimated strains. Each strain haplotype is contained in "{}", where A represents the alternative base and T represents the reference base. As is mentioned before, the actual type of the bases does not influence phylogenetic relationship of the strains; only the difference matters.
* The fourth column represents the relative abundances of the strains in each sample.

## Benchmark

MIST was benchmarked on an in silico synthetic community composed of 10 Staphylococcus aureus strains (NC_xxx) mixed in 10 communities, with each community containing a random subset of strains with random relative abundances (experimental_profile.txt). 

Reads were simulated from the synthetic communities using MetaSim at 5x coverage and mapped back to the Staphylococcus-aureus-specific MetaPhlAn marker gene sequence (NC_xxx.metaphlan.marker.fa), a strategy used in previous strain reconstruction tools such as Strainphlan and ConStrains. To avoid sequencing error, we generated allele count only on sites that are 1) biallelic, and 2) has at least 2x coverage of each allele in at least one of the samples. The allele counts were organized into the input files referenceCount_input.txt and alternativeCount_input.txt. 

We ran MIST on the dataset to infer 2 to 20 strains, and determine the most likely number of strains based on the Bayes factor. The Bayes factor was estimated using an inportance sampling strategy as described in [O'Brian et al. 2014](http://www.genetics.org/content/197/3/925). Finally, the estimated strain haplotypes and their relative abundances were compared to the true haplotypes and relative abundances used to generate the synthetic communities:
![alt text](https://github.com/twinsenzw/MIST/blob/master/Benchmark/bargraph_10in10_5x_6strains.svg)
* Left panel: strain 1~6, the 6 strains inferred by MIST. NC_xx, the 10 strains used to generate the synthetic communities. The phylogeny is generated using the SNP alleles estimated by MIST and the SNP alleles observed in the Staphylococcus-aureus specific MetaPhlAn marker gene sequences (NC_xxx.metaphlan.marker.fa).
* Upper right panel: the estimated relative abundances of the strains. The colors match the color codes in the phylogeny.
* Lower right panel: the true relative abundances of the strains that constitutes the synthetic communities. The colors match the color codes in the phylogeny.

