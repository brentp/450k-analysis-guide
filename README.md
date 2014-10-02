A Practical (And Opinionated) Guide To Analyzing 450K Data
==========================================================

There are a number of options for normalizing and modeling data from the Illumina Infinium HumanMethylation450 BeadChip, this document describes what I believe are good, practical choices for performing analyses. It is not a substitute for individual package documentation, but a broad overview of a workflow. Please fork and fix/add changes.

Design
------

If you have a choice, don't do studies on mixed cell types, like whole-blood. There are smart people [figuring out ways](http://www.biomedcentral.com/1471-2105/13/86) to infer the individual components of a mixture of cell types, but these will remove some of the signal and you'll be less sure of the results.

If you have a choice, don't do studies with, e.g. 3 samples per treatment. This may be fine for cancer, or something that's causing huge changes in cell lines, but it's going to leave you (again) unsure of the results for treatments or conditions with subtle changes.

There are batch and position effects for 450K so randomize your samples across and within plates.

QC and Normalization
--------------------

The `minfi` R package will produce an extensive quality control output from input samples. Examine this, looking out for samples with an odd beta distribution.

Normalization consists of "background correction, color bias (or dye bias) adjustment and Infinium I/II-type bias correction." [cite](http://bib.oxfordjournals.org/content/early/2013/08/27/bib.bbt054.full)

There are a number of methods for normalizing the 450K data. Schalkwyk's group [does a thorough analysis](http://www.ncbi.nlm.nih.gov/pubmed/23631413) of a number of different methods and finds a `dasen` method to be the best performing by their metrics. They provide a plethora of normalization methods in their [wateRmelon package](http://www.bioconductor.org/packages/release/bioc/html/wateRmelon.html).

For my analyses, to date, I use the [SWAN method](http://www.ncbi.nlm.nih.gov/pubmed/22703947&refdoi=10.1186/1471-2164-14-293) as implemented in the [minfi package](http://www.bioconductor.org/packages/release/bioc/html/minfi.html). As noted in the [wateRmelon paper](TODO), I have found that SWAN does do weird things at the extremes so it may be necessary to truncate the Beta values to (0 + 𝛿, 1 - 𝛿). Otherwise values very close to 1 may become outliers after logit-transforming to M-values.

The [minfi paper](http://www.ncbi.nlm.nih.gov/pubmed/24478339) introduces stratified quantile normalization, available as `preprocessQuantile` which may also be a good alternative. But, that paper shows how small the differences in normalization methods really are. The same group has also developed [`preprocessFunnorm`](http://biorxiv.org/content/biorxiv/early/2014/02/23/002956.full.pdf) which uses control probes in a way that can also mediate batch effects (but that function is currently only available in the [devel version](https://bioconductor.org/packages/devel/bioc/html/minfi.html) of `minfi`)....

The authors of [this paper](http://bib.oxfordjournals.org/content/early/2013/08/27/bib.bbt054.full) suggest that: 1) it is better to remove probes with high detection p-values before normalization; 2) high-intensity (M + U signals) type I probes are often unreliable; 3) researchers "do not recommend applying any between-array normalization method to Infinium HumanMethylation450 data". They also suggest a delta-Beta cutoff in addition to p-value cut-offs in downstream analyses. [The paper](http://bib.oxfordjournals.org/content/early/2013/08/27/bib.bbt054.full) is a great overview if you are interested in normalization methods.

PCA/MDS
=======
We often see correlation with some PC's for Plate and for SentrixPosition. These are poster-child batch (and position) effects. You should run PCA (or MDS), plot the first few dimensions and look for outliers. You can also test the principal components for correlation or separation among all of your known laboratory and clinical variables. Usually, gender will be separately perfectly by the 1st or 2nd principal component. This can be used to find sample mixups. See [this software](https://github.com/brentp/clinical-components) and the plot therein. There is also [shiny-methyl](http://f1000research.com/articles/3-175/v1) which allows interactive visualization and is good for QC.

Batch-Effect Removal/Mitigation
===============================
Batch effects *will* be present in microarray data. Randomization across and within plates can mitigate those batch effects, but there may be others that depend somehow on the day of sample preparation, or the technicican who prepared the samples. It is best if the batch effects correspond to known covariates (these can often be found by correlation with principal components as described above). For known batch effects [ComBat](http://biostatistics.oxfordjournals.org/content/8/1/118.abstract) has been [shown to be](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0017238) one of the most effective methods of removing batch effects.

[sva](http://www.bioconductor.org/packages/release/bioc/html/sva.html) is the most common choice for unknown batch effects. It essentially uses the principal components of the residuals of a specified model as potential covariates. Those new covariates can be used as *surrogate variables* in downstream analyses. [PEER](http://www.ncbi.nlm.nih.gov/pubmed/22343431) does something similar with a Bayesian bent. 

Beta vs M
=========
Beta values fall between 0 and 1 and as such are easily interpretable to percent methylation at a given site. To date, most modeling is done on the logit-transformed M values. Currently, this (using M-values) is the most acceptable way to do the modeling, though I believe there are some papers in the pipeline that show no difference using beta-values. 

Modeling
========
Once the data is normalized, and batch effects are handled, the modeling can be quite simple.
Sending the M-values (and the covariates matrix) to [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) is (to my knowledge) the most common way to get a p-value for each probe. The moderated t-statistic should improve power and `limma` is very easy to use. In part, because of this (and the simplified matrix operations), it is most common to put methylation as the dependent variable with disease/treatment as an independent variable even though for most hypotheses (you have a hypothesis, right?) it would make more sense to have disease status and the dependent variable.

There has been recent work showing that a beta-binomial distributed also works well for modelling the beta values [find citation].

In any case, association can be driven by outliers so carefully examine probes with suspiciously low p-values and decide if you want to perform **robust** regression. There is a `robust` option for *limma* that may be appropriate.

Examine the p-value distribution with a [qq-plot](http://en.wikipedia.org/wiki/Q%E2%80%93Q_plot) to look for genomic inflation. It is actually [common](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003678#pgen.1003678.s001) to see inflation but it's good to understand the source. [This paper](http://www.nature.com/ejhg/journal/v19/n7/full/ejhg201139a.html) shows that polygenic inheritance can lead to apparent inflation; by my reading, if we believe that methylation at many sites affects the phenotype, then this may apply. Regardless, [this thread](https://groups.google.com/forum/#!topic/epigenomicsforum/11MXvDqJpFA) has several suggestions for things to look into.

*Remember, you're testing 450K+ sites. So do your multiple-testing correction.*

Differentially Methylated Regions (DMRs)
========================================
We know that methylation has a regional nature--that is, that a CpG site is likely to have similar methylation status as its neighbor due to the machinery that maintains methylation. As such, it makes sense to find regions of methylation that are associated with disease. This can improve power and pertinence over single-probe tests. There are a huge number of methods for finding DMRs. See [this paper](http://biorxiv.org/content/early/2014/07/15/007120) for an enumeration of some of the methods. 

It would be folly to ignore anything from the Irizarry (or Hansen) group on methylation. They propose [bumphunting](http://www.bioconductor.org/packages/release/bioc/html/bumphunter.html) (and see their many papers) which, in short, compares the sum of estimated coefficients in a region to the coefficient sums from data simulated by adding the shuffled residuals of a reduced model (without the covariate of interest) to the fitted values of the reduced model. **NOTE**, that the current implementation actually doesn't permute the residuals and therefore suggests not to use a design matrix with > 1 column. They also loess-smooth and apply an simulate-derived cutoff for the coefficients. In my experience, `bumphunting` can be conservative and the smoothing doesn't always work as I intend but regions that it does report usually look very good.

Although I find the implementation and use less than stellar, the idea from the [A-clustering paper](http://www.ncbi.nlm.nih.gov/pubmed/23990415) is excellent. Rather than modelling all the probes separately and then combining somehow as do most methods including `bumphunting`, they first create clusters of data based on the correlation of the raw-methylation data. (This clustering is unbiased -- without knowing the model to be applied.) They then model the cluster as a group using [Generalized Estimating Equations (GEE's)](http://en.wikipedia.org/wiki/Generalized_estimating_equation). The GEE handles the correlation within samples and the clustering addresses a common problem--how to do multiple-testing correction and/or p-value cutoffs. With tools like bumphunting and other peak-finding tools, we have to choose a cutoff and then report only sites that meet those cutoffs. By modeling clusters, as `Aclust` does, the user can get every tested site and then perform multiple-testing correction as needed. In addition, the user can, for example decide she only wants to test clusters that have more than 8 probes. This can greatly reduce the multiple-testing burden.

Many of the available methods for DMR-finding force one to use a certain class of model (GEE for A-clust) or to make assumptions (Gaussian errors for most methods). I have [developed a method](https://github.com/brentp/combined-pvalues/) that operates only on the p-values. It works by combining the p-values in the manner described in [our paper](http://www.ncbi.nlm.nih.gov/pubmed/22954632) and in the [original work of Katerina Kechris](http://www.ncbi.nlm.nih.gov/pubmed/20812907). With this, we can have more complex study-designs involving **twins, families, or repeated-measures** and use mixed-models to get the p-values and then we can combine the p-values with this method. Given a bed file (columns of chrom, start, end, p-value) of p-values from any study, comb-p can be run with a single command that results in annotated DMRs. **end shameless plug**

Again, there are dozens of methods for finding DMRs; I haven't tried them all.

Mixture of Cells
================
As stated above, if you have a choice (I usually do not), don't work with mixed cell-types like whole-blood or PBMCs. If you must work with mixed-cell types then try to get cell-sorted data as well.
If you can't get cell-sorted data, then you'll have to use [Houseman's method](http://www.biomedcentral.com/1471-2105/13/86) available in `minfi` as `estimateCellCounts` to infer the cell types (from blood). We use this and then plot and test the values for case-vs-control to see which (if any) cell-types are different. Those differences could be driving any and all changes seen in methylation. They could also be spurious differences because the method is not perfect. Any sane reviewer/reader will at least question your whether your results are driven by differential mixtures of cells.

Houseman also has a recent method for [reference-free adjustment for cell mixture](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4016702/) which means you can use it even without a reference of sorted measurements (as is used in `minfi::estimateCellCounts()`). It doesn't give an estimate of different cell-types, but gives you some surrogate variables you can add to your model. 

SNPs
====
The 50-base probes are designed according to a reference sequence. Variants in a sample anywhere in that probe can affect binding and therefore the estimate of methylation. It has been shown that the greatest effect occurs within 2 bases of the CpG [CITE?] so it is common to use the bioconductor package [Illumina450ProbeVariants](http://www.bioconductor.org/packages/release/data/experiment/html/Illumina450ProbeVariants.db.html) to find and remove probes with a variant within 2 bases of the CpG with an allele frequency > 0.05 in your study population. Depending on the population and the stringency (AF, distance to CpG), this can remove about 10% of the 450K probes.

For DMR analyses, our group has left in those sites, but annotated the DMRs by the number of common variants in the probe and in the 2 bases near the measured CpG with the idea that the regional analysis mitigates the effect of a single variant. In fact, if there is a genetic variant driving an association, it's probably valuable to know that. Discarding probes *a priori* may discard information, therefore it may be better to annotate results with variant information rather than discard them out-right. In our studies, I have never seen a change in genomic inflation after removing probes with known SNPs.

For studies with sequence data and 450K data, a more sophisticated analysis can be done to see which exact individuals have a variant within the probe.


Relation to Expression
======================
Your collaborators will be happy if every methylation change that you find occurs in the gene promoter and the expression of that gene goes up as methylation goes down and vice-versa. *That will not happen.*

The most common way to relate methylation to expression globally is to take all methylation changes and find all expression probes (or transcripts for RNA-Seq) within XX KB (where XX is 2, 5, 10, 20 depending on the researcher) and to do a plot of methylation % change vs expression (log2) fold-change. In general, we see an enrichment for pairs in the off-diagonal--that is, sites where methylation has decreased and expression has increased (or vice-versa). But we also see sites where both measurements increase in cases relative to controls. It is becoming less common to expect this inverse relationship without considering other factors, but you will likely deal with this expectation from collaborators.

In addition to examining the correspondence of change, one can also do a correlation of methylation with expression for all samples at a given methylation site and a given transcript. This gives more direct evidence of their association.

Summary
=======
+ don't do mixed cell-types
+ choose a defendable normalization method
+ find and adjust for batch effects
+ perform appropriate modeling
+ find methylation regions
+ relate to expression

References
==========

see refs.bib
