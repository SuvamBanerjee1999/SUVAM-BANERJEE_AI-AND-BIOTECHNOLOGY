

> getwd()
[1] "C:/Users/Suvam/Documents/Module_I"
> setwd("C:/Users/Suvam/Downloads/AI_Omics_Internship_2025")
> dir.create("raw_data")
> dir.create("clean_data")
> dir.create("scripts")
> dir.create("results or Tasks")
> dir.create("plots")
> data<-read.csv(file.choose())
> view(data)
Error in view(data) : could not find function "view"

> library("tidyverse")
── Attaching core tidyverse packages ──────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.5
✔ forcats   1.0.0     ✔ stringr   1.5.1
✔ ggplot2   3.5.2     ✔ tibble    3.3.0
✔ lubridate 1.9.4     ✔ tidyr     1.3.1
✔ purrr     1.1.0     
── Conflicts ────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package to force all conflicts to become errors
> view(data)
> str(data)
'data.frame':	20 obs. of  6 variables:
  $ patient_id: chr  "P001" "P002" "P003" "P004" ...
$ age       : int  34 28 45 39 50 30 41 36 55 29 ...
$ gender    : chr  "Male" "Female" "Female" "Male" ...
$ diagnosis : chr  "Cancer" "Normal" "Cancer" "Normal" ...
$ bmi       : num  22.5 20.3 26.7 23.8 27.1 21.9 25.4 24.2 28.6 19.8 ...
$ smoker    : chr  "Yes" "No" "Yes" "No" ...
> data$gender <- as.factor(data$gender)
> data$diagnosis <- as.factor(data$diagnosis)
> data$smoker <- as.factor(data$smoker)
> str(data)
'data.frame':	20 obs. of  6 variables:
  $ patient_id: chr  "P001" "P002" "P003" "P004" ...
$ age       : int  34 28 45 39 50 30 41 36 55 29 ...
$ gender    : Factor w/ 2 levels "Female","Male": 2 1 1 2 1 2 1 1 2 1 ...
$ diagnosis : Factor w/ 2 levels "Cancer","Normal": 1 2 1 2 1 2 1 2 1 2 ...
$ bmi       : num  22.5 20.3 26.7 23.8 27.1 21.9 25.4 24.2 28.6 19.8 ...
$ smoker    : Factor w/ 2 levels "No","Yes": 2 1 2 1 2 1 2 1 2 1 ...
> data$smoker_num <- ifelse(data$smoker== "Yes", 1, 0)
> view(data)
> # Save cleaned dataset
  > write.csv(data, file = "clean_data/patient_info_clean.csv", row.names = FALSE)