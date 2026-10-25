# Replicate Table 1 from Sundberg

sundberg_data <- c(4, 5, 8, 11, 20, 29, 35, 40, 66, 70)
alpha = 0.05

# for left column, set t_c <- 50
# for right column, set t_c <- 75
t_c <- 50
r <- sum(sundberg_data < t_c)
S <- sum(pmin(sundberg_data, t_c))
MLE <-  r / S

# results for the parametrization rate = theta
C1_Sundberg(r, alpha = alpha) * MLE
C2_Sundberg(r, alpha = alpha) * MLE
C3_Sundberg(r, alpha = alpha) * MLE
C4_Sundberg(r, alpha = alpha) * MLE
C5_Sundberg(r, alpha = alpha) * MLE
C6_Sundberg(r, alpha = alpha) * MLE
C7_Sundberg(r, alpha = alpha) * MLE

# results for the parametrization rate = 1 / theta (like in Sundberg)
# note that C1 and C2 are swapped
1 / (C1_Sundberg(r, alpha = alpha) * MLE) # C2
1 / (C2_Sundberg(r, alpha = alpha) * MLE) # C1
1 / (C3_Sundberg(r, alpha = alpha) * MLE) 
1 / (C4_Sundberg(r, alpha = alpha) * MLE) # wrong
1 / (C5_Sundberg(r, alpha = alpha) * MLE)
1 / (C6_Sundberg(r, alpha = alpha) * MLE) # wrong
1 / (C7_Sundberg(r, alpha = alpha) * MLE) # wrong


# die CI aus CI_Sundberg.R liefern genau die ergebnisse von Sundberg.
# Was ist unterschiedlich?