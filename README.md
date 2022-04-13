# HHstudy_FluMatsumoto2014-15
Replication code and data accompanying: Endo et al. "Fine-scale family structure shapes influenza transmission risk in households: Insights from primary schools in Matsumoto city, 2014/15" [https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007589](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007589)

## Files
* Flu_Matsumoto_HHtransmission.R

R code to produce posterior samples using MCMC.

* HHdata_FluMatsumoto2014-15.csv

Dataset of influenza in households of primary school students in Matsumoto city, Japan

## Dataset description
Guardians of students in all 29 public primary schools in Matsumoto city were asked to respond to a questionnaire on influenza episodes in the 2014/15 season. The dataset includes 10,486 responses.

Each row represents specific patterns of household composition and the existence of influenza episodes. The column "Frequency" shows the number of observation of such household composition and influenza acquisition patterns.

Columns:
* [***.num]: the number of household members (categorised as ***) in the household.
* [***.inf]: the number of household members (categorised as ***) who had influenza during the season.
* [Frequency]: the number of observations (***.num and ***.inf) in the dataset

Household members are categorised into:
* Student: primary school students targeted by the survey
* Sibling: those indicated as brothers or sisters of the students in the response
* Father, Mother: fathers and mothers of the students
* Other: those indicated as other than Student, Sibling, Father or Mother. Breakdowns: grandparents (90%), uncles/aunts (7%), others (3%).

For further details of data collection procedures, see Uchida et al., 2016 (doi: 10.1016/j.pmedr.2016.12.002).

## Dependencies
* R ver. 3.4.1
* {Rcpp} ver. 1.0.0
* {RcppArmadillo} ver. 0.8.300.1.0
* {LaplacesDemon} ver. 16.0.1

## Licence

[MIT](https://github.com/akira-endo/HHstudy_FluMatsumoto2014-15/blob/master/LICENSE)

## Author

[Akira Endo](https://github.com/akira-endo), Mitsuo Uchida, Adam Kucharski, Sebastian Funk
