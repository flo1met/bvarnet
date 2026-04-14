## code to prepare `studentlife` dataset goes here

# Wang, R., Chen, F., Chen, Z., Li, T., Harari, G., Tignor, S., Zhou, X., Ben-Zeev, D., & Campbell, A. T. (2014).
# StudentLife: assessing mental health, academic performance and behavioral trends of college students using smartphones.
# Proceedings of the 2014 ACM International Joint Conference on Pervasive and Ubiquitous Computing, 3–14.
# https://doi.org/10.1145/2632048.2632054

library(openesm)

data <- get_dataset("0004")

studentlife <- data$data

studentlife$anxious <- floor(studentlife$anxious)
studentlife$calm <- floor(studentlife$calm)
studentlife$conventional <- floor(studentlife$conventional)
studentlife$critical <- floor(studentlife$critical)
studentlife$dependable <- floor(studentlife$dependable)
studentlife$stress_level <- floor(studentlife$stress_level)

# binarise: 1 = not happy (original coding 1), 2 = happy (anything > 1)
studentlife$happyornot <- ifelse(studentlife$happyornot > 1, 1L, 0L)

studentlife$difficult_stay_awake <- ifelse(studentlife$difficult_stay_awake > 1, 1L, 0L)

usethis::use_data(studentlife, overwrite = TRUE)
