# Replicate Table 1 from Sundberg

import numpy as np
import sundberg

sundberg_data = np.array([4, 5, 8, 11, 20, 29, 35, 40, 66, 70])
alpha = 0.05

# for left column, set t_c = 50
# for right column, set t_c = 75
t_c = 75
r = int(np.sum(sundberg_data < t_c))
TTT = np.sum(np.minimum(sundberg_data, t_c)) 
MLE =  TTT / r

# results for the parametrization 
print(f"\nRate = 1 / theta\nCensoring time t_c = {t_c}")
for m in range(1, 8):
    c_lower, c_upper = sundberg.CI(method_nr=m, r=r, alpha=alpha, theta_is_rate=False)
    print(f"C{m}: [{c_lower * MLE:.5f}, {c_upper * MLE:.5f}]")


# t_c = 50
# C1: [18.61, 102.59] 
# C2: [9.67, 53.33]
# C3: [15.75, 62.99]
# C4: [14.33, 58.76] no
# C5: [16.92, 68.93]
# C6: [17.47, 72.96]
# C7: [16.69, 66.63]

# t_c = 75
# Rate = 1 / theta (like Sundberg)
# C1: [17.78, 75.75]
# C2: [10.95, 46.65]
# C3: [15.50, 53.53]
# C4: [14.38, 50.59] no
# C5: [16.42, 57.48]
# C6: [16.86, 60.06]
# C7: [16.24, 56.02]


# TODO check C4, values are off