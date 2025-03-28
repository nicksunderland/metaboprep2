# metaboprep data preparation summary report

### Date: 27 March, 2025

metaboprep report relates to:

  * Project: TEST_QC  
  * Platform: Metabolon

The `metaboprep` `R` package performs three operations: 

1. Provides an assessment and summary statistics of the raw metabolomics data.
2. Performs data filtering on the metabolomics data.
3. Provides an assessment and summary statistics of the filtered metabolomics data, particularly in the context of batch variables when available.

This report provides descriptive information for raw and filtered metabolomics data for the project TEST_QC. 








## The data preparation workflow is as follows:



1. Issues can be raised on [GitHub](https://github.com/MRCIEU/metaboprep/issues). 
2. Questions relating to the `metaboprep` pipeline can be directed to [David Hughes: d.a.hughes@bristol.ac.uk](mailto:d.a.hughes@bristol.ac.uk).
3. `metaboprep` is published in [Journal to be determined]() and can be cited as:


# 1. Summary of raw data

## Sample size of TEST_QC data set

<div class="figure" style="text-align: center">
<img src="figure/Sample_size-1.png" alt="plot of chunk Sample_size"  />
<p class="caption">plot of chunk Sample_size</p>
</div>


## Missingness
Missingness is evaluated across samples and features using the original/raw data set. 

### Visual structure of missingness in your raw data set.

<div class="figure" style="text-align: center">
<img src="figure/MissingnessMatrix-1.png" alt="plot of chunk MissingnessMatrix"  />
<p class="caption">plot of chunk MissingnessMatrix</p>
</div>
Figure Legend: Missingness structure across the raw data table. White cells depict missing data. Individuals are in rows, metabolites are in columns. 


### Summary of sample and feature missingness 

<div class="figure" style="text-align: center">
<img src="figure/run.missingness-1.png" alt="plot of chunk run.missingness"  />
<p class="caption">plot of chunk run.missingness</p>
</div>
Figure Legend: Raw data - (a) Distribution of sample missingness with sample mean illustrated by the red vertical line. (b) Distribution of feature missingness sample mean illustrated by the red vertical line. (c) Table of sample and feature missingness percentiles. A tabled version of plot a and b.  (d) Estimates of study samples sizes under various levels of feature missingness.

# 2. Data Filtering 

## Exclusion summary

<div class="figure" style="text-align: center">
<img src="figure/exclusion_table-1.png" alt="plot of chunk exclusion_table"  />
<p class="caption">plot of chunk exclusion_table</p>
</div>
Table Legend: Six primary data filtering exclusion steps were made during the preparation of the data.
(1) Samples with missingness >=80% were excluded. (2) features with missingness >=80% were excluded (xenobiotics are not included in this step). (3) sample exclusions based on the user defined threshold were excluded. (4) feature exclusions based on user defined threshold were excluded (xenobiotics are not included in this step). (5) samples with a total-peak-area or total-sum-abundance that is >= N standard deviations from the mean, where N was defined by the user, were excluded. (6) samples that are >= N standard deviations from the mean on principal component axis 1 and 2, where N was defined by the user, were excluded.


## Metabolite or feature reduction and principal components

A data reduction was carried out to identify a list of representative features for generating a sample principal component analysis. This step reduces the level of inter-correlation in the data to ensure that the principal components are not driven by groups of correlated features.

<div class="figure" style="text-align: center">
<img src="figure/PCA_1_2-1.png" alt="plot of chunk PCA_1_2"  />
<p class="caption">plot of chunk PCA_1_2</p>
</div>
Figure Legend: The data reduction table on the left presents the number of metabolites at each phase of the data reduction (Spearman's correlation distance tree cutting) analysis. On the right principal components 1 and 2 are plotted for all individuals, using the representitve features identified in the data reduction analysis. The red vertical and horizontal lines indicate the standard deviation (SD) cutoffs for identifying individual outliers, which are plotted in red. The standard deviations cuttoff were defined by the user.



# 3. Summary of filtered data

## Sample size (N) 
  * The number of samples in data = 574  
  * The number of features in data = 1062  

## Relative to the raw data
  * 0 samples were filtered out, given the user's criteria.
  * 214 features were filtered out, given the user's criteria.
  * Please review details above and your log file for the number of features and samples excluded and why. 
  

<div class="figure" style="text-align: center">
<img src="figure/sample_overview-1.png" alt="plot of chunk sample_overview"  />
<p class="caption">plot of chunk sample_overview</p>
</div>

Figure Legend: Filtered data summary. Distributions for sample missingness, feature missingness, and total abundance of samples. Row two of the figure provides a Spearman's correlation distance clustering dendrogram highlighting the metabolites used as representative features in blue, the clustering tree cut height is denoted by the horizontal line. Row three provides a summary of the metabolite data reduction in the table, a Scree plot of the variance explained by each PC and a plot of principal component 1 and 2, as derived from the representative metabolites. The Scree plot also identifies the number of PCs estimated to be informative (vertical lines) by the Cattel's Scree Test acceleration factor (red,  n = 2) and Parallel Analysis (green, n = 42). Individuals in the PC plot were clustered into 4 kmeans (k) clusters, using data from PC1 and PC2. The kmeans clustering and color coding is strictly there to help provide some visualization of the major axes of variation in the sample population(s).


## Structure among samples

<div class="figure" style="text-align: center">
<img src="figure/PCApairsplot-1.png" alt="plot of chunk PCApairsplot"  />
<p class="caption">plot of chunk PCApairsplot</p>
</div>

Figure Legend: A matrix (pairs) plot of the top five principal components including demarcations of the 3rd (grey), 4th (orange), and 5th (red) standard deviations from the mean. Samples are color coded as in the summary PC plot above using a kmeans analysis of PC1 and PC2 with a k (number of clusters) set at 4. The choice of k = 4 was not robustly chosen it was a choice of simplicity to help aid visualize variation and sample mobility across the PCs.

## Feature Distributions

### Estimates of normality: W-statistics for raw and log transformed data




<div class="figure" style="text-align: center">
<img src="figure/shapiroW-1.png" alt="plot of chunk shapiroW"  />
<p class="caption">plot of chunk shapiroW</p>
</div>

Figure Legend: Histogram plots of Shapiro W-statistics for raw (left) and log transformed (right) data distributions. A W-statistic value of 1 indicates the sample distribution is perfectly normal and value of 0 indicates it is perfectly uniform. Please note that log transformation of the data *may not* improve the normality of your data.

Analysis details: Of the 1062 features in the data 17 features were excluded from this analysis because of no variation or too few observations (n < 40). Of the remaining 1045 metabolite features, a total of 222 may be considered normally distributed given a Shapiro W-statistic >= 0.95.


## Outliers

Evaluation of the number of samples and features that are outliers across the data.


<div class="figure" style="text-align: center">
<img src="figure/outlier_summary-1.png" alt="plot of chunk outlier_summary"  />
<p class="caption">plot of chunk outlier_summary</p>
</div>
Table Legend: The table reports the average number of outlier values for samples and features in the data set. The minimum, max and other quantiles of the sample and feature distributions are reported as well.

### Notes on outlying samples at each metabolite|feature

There may be extreme outlying observations at individual metabolites|features that have not been accounted for. You may want to:

1. Turn these observations into NAs.
2. Winsorize the data to some maximum value.
3. Rank normalize the data which will place those outliers into the top of the ranked standard normal distribution.
4. Turn these observations into NAs and then impute them along with other missing data in your data set. 


# 4. Variation in filtered data by available variables

### feature missingness

Feature missingness may be influenced by the metabolites' (or features') biology or pathway classification, or your technologies methodology. The figure(s) below provides an illustrative evaluation of the proportion of *feature missigness* as a product of the variable(s) available in the raw data files. 







<div class="figure" style="text-align: center">
<img src="figure/featuremis_1or2-1.png" alt="plot of chunk featuremis_1or2"  />
<p class="caption">plot of chunk featuremis_1or2</p>
</div>



Figure Legend: Box plot illustration(s) of the relationship that feature variables have with feature missingness.


### sample missingness

The figure provides an illustrative evaluation of the proportion of *sample missigness* as a product of sample batch variables provided by your supplier. This is the univariate influence of batch effects on *sample missingness*.


```
 -- A total of 5 possible feature level batch variables were identified in the sample data table. -- 
```

```
 -- No feature level batch variables identified or all were invariable -- 
```

















Figure Legend: Box plot illustration(s) of the relationship that available batch variables have with sample missingness.

## Multivariate evaluation: batch variables


```
[1] " -- No sample level batch variables were provided or all were invariable -- "
```

Table Legend: TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation of the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *sample missingness*. Presence of NA's would indicate that the model is inappropriate.


# 5. Total peak or abundance area (TA) of samples:

The total peak or abundance area (TA) is simply the sum of the abundances measured across all features. TA is one measure that can be used to identify unusual samples given their entire profile. However, the level of missingness in a sample may influence TA. To account for this we:  

1. Evaluate the correlation between TA estimates across all features with PA measured using only those features with complete data (no missingness).  
2. Determine if the batch effects have a measurable impact on TA.


## Relationship with missingness

Correlation between total abundance (TA; at complete features) and missingness

<div class="figure" style="text-align: center">
<img src="figure/TPA_missingness-1.png" alt="plot of chunk TPA_missingness"  />
<p class="caption">plot of chunk TPA_missingness</p>
</div>

Figure Legend: Relationship between total peak area at complete features (x-axis) and sample missingness (y-axis).



## Univariate evaluation: batch effects

The figure below provides an illustrative evaluation of the  *total abundance* (at complete features) as a product of sample batch variables provided by your supplier. 


```
[1] " -- No sample level batch variables were provided or all were invariable -- "
```









Figure Legend: Violin plot illustration(s) of the relationship between total abundance (TA; at complete features) and sample batch variables that are available in your data.


### Multivariate evaluation: batch variables


```
[1] " -- No sample level batch variables were provided or all were invariable -- "
```

Table Legend: TypeII ANOVA: the eta-squared (eta-sq) estimates are an estimation on the percent of variation explained by each independent variable, after accounting for all other variables, as derived from the sum of squares. This is a multivariate evaluation of batch variables on *total peak|abundance area* at complete features.


# 6. Power analysis

Exploration for case/control and continuous outcome data using the filtered data set

Analytical power analysis for both continuous and imbalanced presence/absence correlation analysis.

<div class="figure" style="text-align: center">
<img src="figure/power_exploration-1.png" alt="plot of chunk power_exploration"  />
<p class="caption">plot of chunk power_exploration</p>
</div>

Figure Legend: Simulated  effect sizes (standardized by trait SD) are illustrated by their color in each figure. Figure (A) provides estimates of power for continuous traits with the total sample size on the x-axis and the estimated power on the y-axis. Figure (B) provides estimates of power for presence/absence (or binary) traits in an imbalanced design. The estimated power is on the y-axis. The total sample size is set to 574 and the x-axis depicts the number of individuals present (or absent) for the trait. The effects sizes illustrated here were chosen by running an initial set of simulations which identified effects sizes that would span a broad range of power estimates given the sample population's sample size.

