# this data is not typical

All DNA is stored in a DNA Extraction Plate.   the location of any one sample is defined by the DNA plate ID and the well position.  DNA is transferred from the DNA plate to a quant plate.  
@Sharon Magnuson
 is mainly concerned with the DNA plate and well address, as opposed to the fluor quant plate and well address.
for a given project, there may be 1 to several quant runs due to the limitations of the machine (384 wells).  Consequently, we may have data from either 1 or several quant runs. We generally run duplicate or triplicate quants for each DNA sample.
Each quant run has a standard curve, consequently, each quant run needs the DNA concentrations estimated independently of the other quant runs.  There is substantial variation from quant run to quant run even with identical kit and recipe and samples. Therefore, the concentrations need to be modeled independently.  Because of this variation, when multiple plates need to be used, it's best to spread replicates of the same sample across quant runs.
To process the data, we use 2 scripts.  The first script processes each quant run independently of the other, figures out the best model for rfu vs ng/ul  then estimates the ng/ul for all samples.  It outputs images of the std curves, model ic, and concentrations for each DNA plate id and well address.
The second script eats the output of the first, calculates mean concentration and 95CI, and identifies samples with undesirable qualities for requant.
Samples can be flagged for requanting due to unusually high concentration relative to the other samples, or unusually wide confidence intervals.  the script should output visualaizations of the concentration and 95CI as well as a tabular list of samples (DNA plate id, well address) that require requanting. For those flagged for unusually high confidence intervals, a dilution factor other than 0 : 1 should be suggested which would bring the concentration down in to the range of the other samples.  Any other info that might be useful in identifying why samples were flagged would be useful here. (edited) 


Jason Selwyn
  Friday at 1:09 PM
@Kevin Labrador
 I've made a couple of changes to tests/prj_bird_rota_blue_damselfly_quant_analysis/scripts/summarize_quant_results.R to hopefully better estimate the mean ng/uL of DNA and flag samples of concern.
Estimating DNA Quantity
First the amount of DNA in each sample is now estimated with a model based on however many separate quant plate predictions are made. Specifically:
glmmTMB(ng_per_ul ~ (1 | ID), 
                     dispformula = ~ID,
                     family = Gamma(link = 'log'),
                     dna_amounts)
This models the mean and separate variances for each sample based on however many replicate measurements were made of that sample. I also use a gamma distribution to ensure that all output CIs and means are >0.
Flagging Problem Samples
There are two ways for a sample to be problematic.
Too much DNA is present and it needs to be diluted to not take up too much of the sequencing lane.
Too much uncertainty in our estimate of the amount of DNA.
We make two models to identify these samples. The first is:
lm(log(ng_per_ul_mean) ~ 1, quant_files_summarized) where we estimate the average amount of DNA across all samples and any samples which are outside the 95% prediction interval (i.e. we expect 95% of samples to fall within that interval) This identifies samples which have too much DNA and need to be diluted.
To identify the second class of problem samples we approach it similarly but instead of modelling the mean we model the (upr95%CI - lwr95%CI) / mean. Basically that lets us identify samples which we just don't have a particularly precise measurement of the amount of DNA present relative to the estimated mean. Specifically the model used is: lm(log(ng_per_ul_normspread) ~ 1, data = quant_files_summarized)and then we use the 95% prediction interval to find upper cutoffs which identify samples too variable given their mean.
Here is an updated version of the plot (log10 transformed on the x-axis) with the colors showing how samples were flagged.
@ChrisBird
 
@Kevin Labrador
 
@Sharon Magnuson
 Thoughts? Is this identifying what you would want flagged?
@Sharon Magnuson
 For the dilution information for samples with too much DNA is it useful to say how much water to add to 1uL of sample to get DNA at a concentration equal to the average sample (e.g. if the sample had 77.7 ng/uL the dilution column would say 22.9 indicating that to get to the average of 3.26 ng/uL you should add 22.9 uL of water to 1 uL of sample)? Or is there a different value you would like to target?
@Kevin Labrador
 
@ChrisBird
 I have also modified the jackknife exclusion of standards to 1) include cubic polynomial models and 2) only to exclude a standard if all 4 models identify it as an outlier when fitting the model without the standard vs with it.
@Kevin Labrador
 All the changes I made are in the scripts in tests/prj_bird_rota_blue_damselfly_quant_analysis/scripts  not the main scripts. When its decided how we want to go you can incorporate the changes into the main functions/scripts.
