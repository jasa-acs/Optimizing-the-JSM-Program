## Data

**Abstract (Mandatory)**
The file “inputdata_JSM.RData” includes three datasets used in the analysis: “abstract_info”
concerns the detailed information about the sessions of JSM 2020, “seeded_words” is the table
of the seeded words we used in seeded LDA, and finally “trayvon_martin” is made up of all the
words of Trayvon Martin Corpus we used in data cleaning. The file
“allocation_matrix_list.RData”, instead, is a list of 100 allocation matrices we generated through
our algorithm. We attached it to give the opportunity to avoid running the final loop (commented)
and, consequently, save time.

**Availability (Mandatory)**
No restrictions for the data, they are available.

**Description (Mandatory if data available)**
All data we used were freely available at the official JSM webpage.

## Code

**Abstract (Mandatory)**
The file “jsm_script.R” includes the code. It allows to reproduce in R the analysis we carried out.

**Description (Mandatory)**
R version 3.6.1.
To run the code the R packages dependencies are: text2vec 0.5.1, data.table 1.13.4, stringr
1.4.0, stopwords 1 .0, corpus 0.10.0, slam 0.1- 45 , topicmodels 0.2- 8 and ggplot2 3.3.2.

## Instructions for Use

**Reproducibility (Mandatory)**
Figure 1 and Tables 1,2, 3 ,5 are reproducible (workflow information) from the code.
