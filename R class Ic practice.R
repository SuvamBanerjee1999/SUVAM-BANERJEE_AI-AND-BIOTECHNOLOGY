cholesterol <- 230
Systolic_bp <- 120

if (cholesterol > 240) {
  print("High Cholesterol")
}

if (Systolic_bp < 120) {
  print("Blood Pressure is normal")
} else { 
  print("Blood Pressure is high")
}
raw_data <- read.csv(file.choose())
str(raw_data)
patient_info_copy <- raw_data

factor_cols <- c("gender", "diagnosis", "smoker")
for (col in factor_cols) {
     patient_info_copy[[col]] <- as.factor(patient_info_copy[[col]])
}
str(patient_info_copy)
raw_mdata <- read.csv(file.choose())
str(raw_mdata)
metadata_copy <- raw_mdata
factor_mcols <- c("height", "gender")

for (col in factor_mcols) {
  metadata_copy[[col]] <- as.factor(metadata_copy[[col]])
}
str(metadata_copy)
patient_info_copy$smoker <- ifelse(patient_info_copy$smoker=="Yes",1L,0L)
str(patient_info_copy)
metadata_copy$gender <- ifelse(metadata_copy$gender=="Male",1L,0L)
str(metadata_copy)
str(raw_data)
str(raw_mdata)







