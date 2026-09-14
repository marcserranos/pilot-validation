# Per-allele ancestry-centroid geometry -- gene HLA-B, cohort `lr (Immuannot long-read, full)`

Ternary components: AFR/EUR/AMR. Tetrahedron components: AFR/EUR/AMR/EAS. Primary centroid renorm: `per_person` (both variants always computed; the alt column shows the other one).

**No minimum-carrier-count floor is applied to centroid computation** -- an allele with 1 carrier is shown with a centroid equal to that person's own admixture proportions. n_carriers plus marker alpha communicate confidence visually instead of a hard cutoff. Bootstrap enrichment (continuous) and Fisher/chi2 enrichment (discrete, on ancestry_pred labels) are two independent significance reads, reported side by side -- neither is used to decide what gets plotted.


## Ternary-space alleles (sorted by bootstrap enrichment p-value)

| allele | kind | n_carriers | centroid | centroid (alt renorm) | bootstrap z | bootstrap p | chi2 p | confidence_tier |
|---|---|---|---|---|---|---|---|---|
| 07:02 | known | 1359 | AFR=0.334, EUR=0.412, AMR=0.254 | AFR=0.349, EUR=0.400, AMR=0.251 | 5.356 | 0.000 | 0.000 |  |
| 07:05 | known | 118 | AFR=0.020, EUR=0.822, AMR=0.157 | AFR=0.017, EUR=0.807, AMR=0.176 | 14.548 | 0.000 | 0.000 |  |
| 07:06 | known | 133 | AFR=0.471, EUR=0.191, AMR=0.337 | AFR=0.504, EUR=0.161, AMR=0.336 | 5.284 | 0.000 | 0.000 |  |
| 08:01 | known | 1029 | AFR=0.256, EUR=0.503, AMR=0.241 | AFR=0.268, EUR=0.489, AMR=0.242 | 13.004 | 0.000 | 0.000 |  |
| 13:01 | known | 115 | AFR=0.072, EUR=0.544, AMR=0.384 | AFR=0.177, EUR=0.429, AMR=0.393 | 7.009 | 0.000 | 0.000 |  |
| 13:02 | known | 380 | AFR=0.175, EUR=0.548, AMR=0.277 | AFR=0.201, EUR=0.490, AMR=0.310 | 10.398 | 0.000 | 0.000 |  |
| 14:02 | known | 693 | AFR=0.204, EUR=0.469, AMR=0.326 | AFR=0.228, EUR=0.415, AMR=0.357 | 9.240 | 0.000 | 0.000 |  |
| 14:03 | known | 16 | AFR=0.916, EUR=0.005, AMR=0.079 | AFR=0.918, EUR=0.004, AMR=0.078 | 5.985 | 0.000 | 0.000 |  |
| 15:01 | known | 566 | AFR=0.168, EUR=0.539, AMR=0.292 | AFR=0.169, EUR=0.535, AMR=0.296 | 13.165 | 0.000 | 0.000 |  |
| 15:02 | known | 200 | AFR=0.096, EUR=0.478, AMR=0.425 | AFR=0.133, EUR=0.551, AMR=0.316 | 8.129 | 0.000 | 0.000 |  |
| 15:03 | known | 477 | AFR=0.763, EUR=0.036, AMR=0.201 | AFR=0.776, EUR=0.025, AMR=0.198 | 31.283 | 0.000 | 0.000 |  |
| 15:10 | known | 266 | AFR=0.755, EUR=0.058, AMR=0.187 | AFR=0.779, EUR=0.038, AMR=0.183 | 21.467 | 0.000 | 0.000 |  |
| 15:11 | known | 33 | AFR=0.047, EUR=0.571, AMR=0.382 | AFR=0.038, EUR=0.479, AMR=0.484 | 3.456 | 0.000 | 0.000 |  |
| 15:15 | known | 35 | AFR=0.032, EUR=0.013, AMR=0.955 | AFR=0.030, EUR=0.009, AMR=0.961 | 10.646 | 0.000 | 0.000 |  |
| 15:16 | known | 167 | AFR=0.714, EUR=0.071, AMR=0.215 | AFR=0.738, EUR=0.049, AMR=0.212 | 15.427 | 0.000 | 0.000 |  |
| 15:17 | known | 158 | AFR=0.200, EUR=0.540, AMR=0.261 | AFR=0.253, EUR=0.437, AMR=0.310 | 5.237 | 0.000 | 0.000 |  |
| 15:25 | known | 52 | AFR=0.167, EUR=0.024, AMR=0.810 | AFR=0.240, EUR=0.010, AMR=0.751 | 10.385 | 0.000 | 0.000 |  |
| 15:27 | known | 17 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 7.168 | 0.000 | 0.000 |  |
| 15:30 | known | 17 | AFR=0.024, EUR=0.007, AMR=0.970 | AFR=0.022, EUR=0.006, AMR=0.972 | 6.775 | 0.000 | 0.000 |  |
| 15:35 | known | 46 | AFR=0.036, EUR=0.656, AMR=0.307 | AFR=0.109, EUR=0.391, AMR=0.500 | 5.598 | 0.000 | 0.000 |  |
| 18:01 | known | 753 | AFR=0.252, EUR=0.472, AMR=0.276 | AFR=0.286, EUR=0.423, AMR=0.291 | 7.571 | 0.000 | 0.000 |  |
| 18:02 | known | 12 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 5.323 | 0.000 | 0.000 |  |
| 27:02 | known | 61 | AFR=0.072, EUR=0.634, AMR=0.293 | AFR=0.081, EUR=0.587, AMR=0.332 | 5.609 | 0.000 | 0.000 |  |
| 27:03 | known | 31 | AFR=0.741, EUR=0.080, AMR=0.179 | AFR=0.801, EUR=0.016, AMR=0.183 | 5.787 | 0.000 | 0.000 |  |
| 27:04 | known | 42 | AFR=0.000, EUR=0.556, AMR=0.444 | AFR=0.000, EUR=0.969, AMR=0.031 | 4.996 | 0.000 | 0.000 |  |
| 27:05 | known | 374 | AFR=0.171, EUR=0.457, AMR=0.372 | AFR=0.172, EUR=0.451, AMR=0.377 | 7.698 | 0.000 | 0.000 |  |
| 27:07 | known | 35 | AFR=0.017, EUR=0.744, AMR=0.239 | AFR=0.014, EUR=0.714, AMR=0.271 | 5.872 | 0.000 | 0.000 |  |
| 35:02 | known | 299 | AFR=0.039, EUR=0.721, AMR=0.240 | AFR=0.040, EUR=0.673, AMR=0.287 | 18.999 | 0.000 | 0.000 |  |
| 35:03 | known | 445 | AFR=0.058, EUR=0.592, AMR=0.350 | AFR=0.063, EUR=0.554, AMR=0.384 | 17.695 | 0.000 | 0.000 |  |
| 35:04 | known | 14 | AFR=0.066, EUR=0.005, AMR=0.928 | AFR=0.065, EUR=0.005, AMR=0.930 | 5.559 | 0.000 | 0.000 |  |
| 35:05 | known | 106 | AFR=0.102, EUR=0.421, AMR=0.477 | AFR=0.228, EUR=0.213, AMR=0.559 | 5.692 | 0.000 | 0.000 |  |
| 35:08 | known | 115 | AFR=0.044, EUR=0.652, AMR=0.304 | AFR=0.050, EUR=0.578, AMR=0.372 | 9.405 | 0.000 | 0.000 |  |
| 35:12 | known | 92 | AFR=0.041, EUR=0.023, AMR=0.936 | AFR=0.038, EUR=0.020, AMR=0.942 | 17.439 | 0.000 | 0.000 |  |
| 35:14 | known | 11 | AFR=0.010, EUR=0.005, AMR=0.986 | AFR=0.008, EUR=0.004, AMR=0.988 | 5.367 | 0.000 | 0.000 |  |
| 35:17 | known | 37 | AFR=0.009, EUR=0.031, AMR=0.960 | AFR=0.008, EUR=0.030, AMR=0.962 | 10.735 | 0.000 | 0.000 |  |
| 35:43 | known | 27 | AFR=0.011, EUR=0.017, AMR=0.972 | AFR=0.010, EUR=0.015, AMR=0.975 | 9.004 | 0.000 | 0.000 |  |
| 37:01 | known | 208 | AFR=0.158, EUR=0.594, AMR=0.248 | AFR=0.174, EUR=0.576, AMR=0.251 | 9.010 | 0.000 | 0.000 |  |
| 38:01 | known | 508 | AFR=0.045, EUR=0.714, AMR=0.241 | AFR=0.051, EUR=0.685, AMR=0.264 | 25.091 | 0.000 | 0.000 |  |
| 38:02 | known | 173 | AFR=0.043, EUR=0.509, AMR=0.448 | AFR=0.073, EUR=0.412, AMR=0.514 | 10.112 | 0.000 | 0.000 |  |
| 39:01 | known | 174 | AFR=0.145, EUR=0.429, AMR=0.426 | AFR=0.154, EUR=0.398, AMR=0.448 | 5.838 | 0.000 | 0.000 |  |
| 39:02 | known | 20 | AFR=0.050, EUR=0.018, AMR=0.933 | AFR=0.048, EUR=0.016, AMR=0.936 | 6.887 | 0.000 | 0.000 |  |
| 39:05 | known | 127 | AFR=0.047, EUR=0.032, AMR=0.922 | AFR=0.045, EUR=0.023, AMR=0.933 | 20.606 | 0.000 | 0.000 |  |
| 39:06 | known | 146 | AFR=0.041, EUR=0.151, AMR=0.808 | AFR=0.039, EUR=0.139, AMR=0.822 | 17.781 | 0.000 | 0.000 |  |
| 39:08 | known | 19 | AFR=0.027, EUR=0.011, AMR=0.962 | AFR=0.026, EUR=0.011, AMR=0.963 | 7.087 | 0.000 | 0.000 |  |
| 39:09 | known | 12 | AFR=0.005, EUR=0.188, AMR=0.807 | AFR=0.006, EUR=0.079, AMR=0.915 | 3.797 | 0.000 | 0.000 |  |
| 39:10 | known | 64 | AFR=0.765, EUR=0.078, AMR=0.157 | AFR=0.842, EUR=0.023, AMR=0.135 | 10.043 | 0.000 | 0.000 |  |
| 39:24 | known | 12 | AFR=0.087, EUR=0.872, AMR=0.041 | AFR=0.155, EUR=0.793, AMR=0.052 | 3.869 | 0.000 | 0.000 |  |
| 40:01 | known | 678 | AFR=0.180, EUR=0.550, AMR=0.270 | AFR=0.195, EUR=0.557, AMR=0.248 | 14.701 | 0.000 | 0.000 |  |
| 40:02 | known | 383 | AFR=0.083, EUR=0.278, AMR=0.639 | AFR=0.082, EUR=0.244, AMR=0.674 | 20.096 | 0.000 | 0.000 |  |
| 40:03 | known | 12 | AFR=0.036, EUR=0.132, AMR=0.832 | AFR=0.034, EUR=0.103, AMR=0.863 | 4.031 | 0.000 | 0.001 |  |
| 40:04 | known | 22 | AFR=0.059, EUR=0.010, AMR=0.931 | AFR=0.058, EUR=0.009, AMR=0.933 | 7.425 | 0.000 | 0.000 |  |
| 40:05 | known | 29 | AFR=0.061, EUR=0.009, AMR=0.930 | AFR=0.056, EUR=0.008, AMR=0.936 | 8.893 | 0.000 | 0.000 |  |
| 40:08 | known | 12 | AFR=0.011, EUR=0.040, AMR=0.949 | AFR=0.011, EUR=0.038, AMR=0.951 | 5.349 | 0.000 | 0.000 |  |
| 41:01 | known | 196 | AFR=0.144, EUR=0.519, AMR=0.337 | AFR=0.188, EUR=0.392, AMR=0.421 | 6.798 | 0.000 | 0.000 |  |
| 41:02 | known | 133 | AFR=0.222, EUR=0.586, AMR=0.192 | AFR=0.277, EUR=0.505, AMR=0.218 | 6.408 | 0.000 | 0.000 |  |
| 42:01 | known | 365 | AFR=0.790, EUR=0.032, AMR=0.178 | AFR=0.806, EUR=0.022, AMR=0.173 | 28.653 | 0.000 | 0.000 |  |
| 42:02 | known | 54 | AFR=0.541, EUR=0.076, AMR=0.384 | AFR=0.562, EUR=0.056, AMR=0.382 | 5.588 | 0.000 | 0.000 |  |
| 44:02 | known | 736 | AFR=0.136, EUR=0.568, AMR=0.296 | AFR=0.142, EUR=0.552, AMR=0.306 | 18.004 | 0.000 | 0.000 |  |
| 44:03 | known | 1139 | AFR=0.320, EUR=0.307, AMR=0.373 | AFR=0.339, EUR=0.279, AMR=0.381 | 5.718 | 0.000 | 0.000 |  |
| 44:10 | known | 13 | AFR=0.950, EUR=0.005, AMR=0.045 | AFR=0.950, EUR=0.005, AMR=0.045 | 5.340 | 0.000 | 0.000 |  |
| 44:27 | known | 20 | AFR=0.053, EUR=0.841, AMR=0.107 | AFR=0.049, EUR=0.847, AMR=0.104 | 4.834 | 0.000 | 0.000 |  |
| 45:01 | known | 421 | AFR=0.631, EUR=0.091, AMR=0.278 | AFR=0.648, EUR=0.075, AMR=0.276 | 20.578 | 0.000 | 0.000 |  |
| 46:01 | known | 216 | AFR=0.032, EUR=0.512, AMR=0.456 | AFR=0.051, EUR=0.482, AMR=0.466 | 12.194 | 0.000 | 0.000 |  |
| 48:01 | known | 167 | AFR=0.063, EUR=0.198, AMR=0.739 | AFR=0.061, EUR=0.124, AMR=0.815 | 16.201 | 0.000 | 0.000 |  |
| 50:01 | known | 328 | AFR=0.210, EUR=0.523, AMR=0.266 | AFR=0.262, EUR=0.439, AMR=0.299 | 7.310 | 0.000 | 0.000 |  |
| 51:01 | known | 923 | AFR=0.199, EUR=0.431, AMR=0.369 | AFR=0.223, EUR=0.376, AMR=0.401 | 10.312 | 0.000 | 0.000 |  |
| 51:02 | known | 66 | AFR=0.125, EUR=0.047, AMR=0.827 | AFR=0.135, EUR=0.019, AMR=0.846 | 12.442 | 0.000 | 0.000 |  |
| 51:06 | known | 31 | AFR=0.093, EUR=0.103, AMR=0.804 | AFR=0.130, EUR=0.174, AMR=0.696 | 6.890 | 0.000 | 0.000 |  |
| 52:01 | known | 595 | AFR=0.234, EUR=0.429, AMR=0.336 | AFR=0.274, EUR=0.359, AMR=0.367 | 4.973 | 0.000 | 0.000 |  |
| 53:01 | known | 825 | AFR=0.733, EUR=0.070, AMR=0.197 | AFR=0.762, EUR=0.043, AMR=0.196 | 37.354 | 0.000 | 0.000 |  |
| 54:01 | known | 84 | AFR=0.023, EUR=0.625, AMR=0.352 | AFR=0.034, EUR=0.756, AMR=0.211 | 7.870 | 0.000 | 0.000 |  |
| 55:01 | known | 233 | AFR=0.129, EUR=0.590, AMR=0.281 | AFR=0.145, EUR=0.557, AMR=0.298 | 10.017 | 0.000 | 0.000 |  |
| 55:02 | known | 49 | AFR=0.000, EUR=0.684, AMR=0.316 | AFR=0.000, EUR=0.872, AMR=0.128 | 6.538 | 0.000 | 0.000 |  |
| 57:01 | known | 497 | AFR=0.123, EUR=0.649, AMR=0.228 | AFR=0.125, EUR=0.644, AMR=0.231 | 19.359 | 0.000 | 0.000 |  |
| 57:02 | known | 30 | AFR=0.637, EUR=0.163, AMR=0.200 | AFR=0.705, EUR=0.102, AMR=0.193 | 3.854 | 0.000 | 0.000 |  |
| 57:03 | known | 270 | AFR=0.723, EUR=0.046, AMR=0.231 | AFR=0.744, EUR=0.028, AMR=0.228 | 21.040 | 0.000 | 0.000 |  |
| 58:01 | known | 704 | AFR=0.502, EUR=0.255, AMR=0.243 | AFR=0.575, EUR=0.190, AMR=0.235 | 14.528 | 0.000 | 0.000 |  |
| 58:02 | known | 257 | AFR=0.816, EUR=0.014, AMR=0.169 | AFR=0.826, EUR=0.010, AMR=0.164 | 25.263 | 0.000 | 0.000 |  |
| 78:01 | known | 66 | AFR=0.763, EUR=0.034, AMR=0.203 | AFR=0.772, EUR=0.030, AMR=0.198 | 10.835 | 0.000 | 0.000 |  |
| 81:01 | known | 142 | AFR=0.796, EUR=0.035, AMR=0.169 | AFR=0.808, EUR=0.027, AMR=0.165 | 17.096 | 0.000 | 0.000 |  |
| 82:01 | known | 17 | AFR=0.837, EUR=0.007, AMR=0.156 | AFR=0.840, EUR=0.007, AMR=0.154 | 5.081 | 0.000 | 0.000 |  |
| HLA-B_nov_2a9ca252 | novel | 52 | AFR=0.661, EUR=0.019, AMR=0.320 | AFR=0.667, EUR=0.014, AMR=0.319 | 7.710 | 0.000 | 0.000 | high |
| HLA-B_nov_d15615ce | novel | 38 | AFR=0.056, EUR=0.018, AMR=0.926 | AFR=0.054, EUR=0.013, AMR=0.932 | 10.448 | 0.000 | 0.000 | high |
| HLA-B_nov_5b19fb36 | novel | 40 | AFR=0.000, EUR=0.662, AMR=0.338 | AFR=0.000, EUR=0.973, AMR=0.027 | 5.357 | 0.000 | 0.000 | high |
| HLA-B_nov_a3e9f87a | novel | 11 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 5.521 | 0.000 | 0.000 | high |
| HLA-B_nov_74d648cd | novel | 39 | AFR=0.012, EUR=0.336, AMR=0.652 | AFR=0.007, EUR=0.309, AMR=0.684 | 5.819 | 0.000 | 0.000 | high |
| HLA-B_nov_7752bb66 | novel | 60 | AFR=0.069, EUR=0.598, AMR=0.333 | AFR=0.079, EUR=0.564, AMR=0.358 | 5.300 | 0.000 | 0.000 | high |
| HLA-B_nov_53f8de08 | novel | 23 | AFR=0.763, EUR=0.121, AMR=0.116 | AFR=0.807, EUR=0.082, AMR=0.111 | 4.940 | 0.000 | 0.000 | high |
| HLA-B_nov_2a7f9a38 | novel | 45 | AFR=0.722, EUR=0.016, AMR=0.261 | AFR=0.732, EUR=0.014, AMR=0.254 | 7.727 | 0.000 | 0.000 | high |
| HLA-B_nov_e056d52f | novel | 39 | AFR=0.033, EUR=0.565, AMR=0.402 | AFR=0.040, EUR=0.555, AMR=0.405 | 4.128 | 0.000 | 0.000 | high |
| HLA-B_nov_2670c4c2 | novel | 57 | AFR=0.059, EUR=0.455, AMR=0.486 | AFR=0.064, EUR=0.398, AMR=0.538 | 4.647 | 0.000 | 0.000 | high |
| HLA-B_nov_ab1351b4 | novel | 77 | AFR=0.752, EUR=0.051, AMR=0.197 | AFR=0.759, EUR=0.048, AMR=0.193 | 10.901 | 0.000 | 0.000 | high |
| HLA-B_nov_7571b4cf | novel | 58 | AFR=0.099, EUR=0.657, AMR=0.243 | AFR=0.069, EUR=0.724, AMR=0.207 | 5.516 | 0.000 | 0.000 | high |
| HLA-B_nov_def7a880 | novel | 96 | AFR=0.676, EUR=0.077, AMR=0.247 | AFR=0.710, EUR=0.046, AMR=0.244 | 10.233 | 0.000 | 0.000 | high |
| HLA-B_nov_abf2f906 | novel | 26 | AFR=0.778, EUR=0.143, AMR=0.078 | AFR=0.819, EUR=0.114, AMR=0.068 | 5.627 | 0.000 | 0.000 | high |
| HLA-B_nov_cd0f9284 | novel | 67 | AFR=0.527, EUR=0.228, AMR=0.245 | AFR=0.588, EUR=0.168, AMR=0.244 | 3.645 | 0.000 | 0.003 | high |
| HLA-B_nov_18232e17 | novel | 35 | AFR=0.869, EUR=0.009, AMR=0.122 | AFR=0.876, EUR=0.008, AMR=0.117 | 8.971 | 0.000 | 0.000 | high |
| HLA-B_nov_af51adc0 | novel | 23 | AFR=0.030, EUR=0.113, AMR=0.857 | AFR=0.037, EUR=0.125, AMR=0.838 | 6.497 | 0.000 | 0.000 | high |
| HLA-B_nov_118d7054 | novel | 44 | AFR=0.116, EUR=0.662, AMR=0.223 | AFR=0.112, EUR=0.674, AMR=0.214 | 4.768 | 0.000 | 0.000 | high |
| HLA-B_nov_bca4b2b9 | novel | 28 | AFR=0.710, EUR=0.013, AMR=0.277 | AFR=0.715, EUR=0.012, AMR=0.273 | 5.534 | 0.000 | 0.000 | high |
| HLA-B_nov_8a2e0836 | novel | 28 | AFR=0.000, EUR=0.957, AMR=0.043 | AFR=0.000, EUR=0.915, AMR=0.085 | 8.296 | 0.000 | 0.000 | high |
| HLA-B_nov_1d3870fc | novel | 44 | AFR=0.846, EUR=0.010, AMR=0.144 | AFR=0.861, EUR=0.008, AMR=0.131 | 10.010 | 0.000 | 0.000 | high |
| HLA-B_nov_d5ae3969 | novel | 60 | AFR=0.069, EUR=0.759, AMR=0.172 | AFR=0.091, EUR=0.793, AMR=0.116 | 8.382 | 0.000 | 0.000 | high |
| HLA-B_nov_a58ba6da | novel | 16 | AFR=0.003, EUR=0.005, AMR=0.992 | AFR=0.003, EUR=0.004, AMR=0.992 | 7.031 | 0.000 | 0.000 | high |
| 15:29 | known | 10 | AFR=0.024, EUR=0.951, AMR=0.024 | AFR=0.024, EUR=0.951, AMR=0.024 | 4.197 | 0.001 | 0.000 |  |
| 35:01 | known | 1294 | AFR=0.344, EUR=0.323, AMR=0.333 | AFR=0.370, EUR=0.279, AMR=0.351 | 3.426 | 0.001 | 0.000 |  |
| 35:16 | known | 7 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 4.285 | 0.001 | 0.000 |  |
| 51:08 | known | 25 | AFR=0.053, EUR=0.714, AMR=0.233 | AFR=0.069, EUR=0.601, AMR=0.329 | 3.824 | 0.001 | 0.000 |  |
| HLA-B_nov_c4e5fc19 | novel | 18 | AFR=0.066, EUR=0.794, AMR=0.140 | AFR=0.092, EUR=0.773, AMR=0.136 | 3.828 | 0.001 | 0.000 | high |
| HLA-B_nov_3891497d | novel | 6 | AFR=0.015, EUR=0.000, AMR=0.985 | AFR=0.015, EUR=0.000, AMR=0.985 | 3.502 | 0.001 | 0.001 | high |
| 73:01 | known | 27 | AFR=0.005, EUR=0.614, AMR=0.381 | AFR=0.006, EUR=0.512, AMR=0.482 | 3.590 | 0.002 | 0.000 |  |
| HLA-B_nov_860e6c9d | novel | 7 | AFR=0.007, EUR=0.054, AMR=0.939 | AFR=0.006, EUR=0.047, AMR=0.947 | 3.717 | 0.002 | 0.000 | high |
| HLA-B_nov_d4f09b7d | novel | 5 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 3.062 | 0.002 | 0.004 | high |
| 40:06 | known | 265 | AFR=0.235, EUR=0.463, AMR=0.302 | AFR=0.310, EUR=0.415, AMR=0.275 | 3.516 | 0.003 | 0.000 |  |
| 48:02 | known | 9 | AFR=0.056, EUR=0.046, AMR=0.897 | AFR=0.054, EUR=0.030, AMR=0.916 | 3.776 | 0.003 | 0.000 |  |
| HLA-B_nov_5957fead | novel | 79 | AFR=0.475, EUR=0.217, AMR=0.308 | AFR=0.498, EUR=0.185, AMR=0.317 | 3.218 | 0.003 | 0.000 | high |
| HLA-B_nov_ba21176a | novel | 21 | AFR=0.112, EUR=0.694, AMR=0.193 | AFR=0.130, EUR=0.669, AMR=0.201 | 2.997 | 0.004 | 0.005 | high |
| 15:39 | known | 4 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 2.525 | 0.004 | 0.018 |  |
| 15:21 | known | 37 | AFR=0.077, EUR=0.532, AMR=0.391 | AFR=0.157, EUR=0.478, AMR=0.366 | 2.963 | 0.005 | 0.000 |  |
| 18:04 | known | 7 | AFR=0.108, EUR=0.032, AMR=0.860 | AFR=0.107, EUR=0.031, AMR=0.862 | 3.030 | 0.005 | 0.000 |  |
| 35:20 | known | 10 | AFR=0.183, EUR=0.015, AMR=0.801 | AFR=0.184, EUR=0.014, AMR=0.802 | 3.347 | 0.005 | 0.000 |  |
| 39:11 | known | 6 | AFR=0.042, EUR=0.009, AMR=0.949 | AFR=0.039, EUR=0.009, AMR=0.952 | 3.218 | 0.005 | 0.001 |  |
| 44:05 | known | 32 | AFR=0.065, EUR=0.600, AMR=0.335 | AFR=0.071, EUR=0.556, AMR=0.373 | 3.459 | 0.005 | 0.000 |  |
| 82:02 | known | 6 | AFR=0.970, EUR=0.005, AMR=0.024 | AFR=0.971, EUR=0.005, AMR=0.024 | 3.332 | 0.005 | 0.004 |  |
| HLA-B_nov_2b2961f7 | novel | 6 | AFR=0.970, EUR=0.005, AMR=0.024 | AFR=0.971, EUR=0.005, AMR=0.024 | 3.332 | 0.005 | 0.004 | high |
| 15:20 | known | 9 | AFR=0.163, EUR=0.018, AMR=0.819 | AFR=0.164, EUR=0.018, AMR=0.818 | 3.120 | 0.006 | 0.000 |  |
| 15:37 | known | 7 | AFR=0.859, EUR=0.069, AMR=0.072 | AFR=0.899, EUR=0.047, AMR=0.054 | 2.925 | 0.006 | 0.012 |  |
| 56:01 | known | 113 | AFR=0.241, EUR=0.509, AMR=0.251 | AFR=0.264, EUR=0.485, AMR=0.251 | 3.041 | 0.006 | 0.029 |  |
| 67:01 | known | 9 | AFR=0.060, EUR=0.119, AMR=0.821 | AFR=0.060, EUR=0.119, AMR=0.821 | 3.025 | 0.006 | 0.000 |  |
| 15:05 | known | 15 | AFR=0.173, EUR=0.123, AMR=0.704 | AFR=0.173, EUR=0.123, AMR=0.704 | 2.960 | 0.007 | 0.000 |  |
| HLA-B_nov_59bcae82 | novel | 30 | AFR=0.192, EUR=0.644, AMR=0.164 | AFR=0.225, EUR=0.608, AMR=0.168 | 3.058 | 0.007 | 0.001 | high |
| HLA-B_nov_f9c021a8 | novel | 8 | AFR=0.831, EUR=0.007, AMR=0.162 | AFR=0.842, EUR=0.006, AMR=0.152 | 3.051 | 0.008 | 0.045 | high |
| 41:03 | known | 6 | AFR=0.922, EUR=0.000, AMR=0.078 | AFR=0.924, EUR=0.000, AMR=0.076 | 2.970 | 0.009 | 0.004 |  |
| 35:11 | known | 6 | AFR=0.098, EUR=0.010, AMR=0.892 | AFR=0.095, EUR=0.010, AMR=0.894 | 2.782 | 0.010 | 0.001 |  |
| HLA-B_nov_79772159 | novel | 6 | AFR=0.098, EUR=0.010, AMR=0.892 | AFR=0.095, EUR=0.010, AMR=0.894 | 2.782 | 0.010 | 0.001 | high |
| 40:27 | known | 7 | AFR=0.126, EUR=0.068, AMR=0.806 | AFR=0.123, EUR=0.067, AMR=0.810 | 2.521 | 0.012 | 0.006 |  |
| 57:04 | known | 18 | AFR=0.651, EUR=0.115, AMR=0.234 | AFR=0.687, EUR=0.068, AMR=0.245 | 2.721 | 0.014 | 0.006 |  |
| HLA-B_nov_207ec506 | novel | 39 | AFR=0.358, EUR=0.158, AMR=0.484 | AFR=0.374, EUR=0.113, AMR=0.514 | 2.650 | 0.014 | 0.000 | high |
| 56:04 | known | 15 | AFR=0.000, EUR=0.389, AMR=0.611 | AFR=0.000, EUR=0.511, AMR=0.489 | 2.503 | 0.017 | 0.000 |  |
| HLA-B_nov_f4ec8517 | novel | 87 | AFR=0.184, EUR=0.423, AMR=0.393 | AFR=0.195, EUR=0.388, AMR=0.417 | 2.273 | 0.017 | 0.101 | high |
| 40:16 | known | 8 | AFR=0.747, EUR=0.001, AMR=0.252 | AFR=0.748, EUR=0.001, AMR=0.251 | 2.440 | 0.018 | 0.036 |  |
| HLA-B_nov_d24423c2 | novel | 19 | AFR=0.412, EUR=0.055, AMR=0.532 | AFR=0.420, EUR=0.051, AMR=0.530 | 2.514 | 0.018 | 0.000 | high |
| 35:10 | known | 3 | AFR=0.003, EUR=0.000, AMR=0.997 | AFR=0.003, EUR=0.000, AMR=0.997 | 1.944 | 0.019 | 0.070 |  |
| 15:08 | known | 12 | AFR=0.346, EUR=0.030, AMR=0.624 | AFR=0.491, EUR=0.044, AMR=0.465 | 2.375 | 0.020 | 0.000 |  |
| 15:04 | known | 5 | AFR=0.068, EUR=0.036, AMR=0.896 | AFR=0.066, EUR=0.033, AMR=0.901 | 2.315 | 0.031 | 0.005 |  |
| 15:07 | known | 11 | AFR=0.171, EUR=0.144, AMR=0.685 | AFR=0.167, EUR=0.146, AMR=0.687 | 2.120 | 0.032 | 0.041 |  |
| 40:12 | known | 6 | AFR=0.783, EUR=0.140, AMR=0.077 | AFR=0.947, EUR=0.013, AMR=0.040 | 1.802 | 0.039 | 0.013 |  |
| 18:145 | known | 8 | AFR=0.702, EUR=0.053, AMR=0.245 | AFR=0.719, EUR=0.041, AMR=0.241 | 1.928 | 0.040 | 0.222 |  |
| HLA-B_nov_5484b26e | novel | 9 | AFR=0.059, EUR=0.745, AMR=0.196 | AFR=0.059, EUR=0.745, AMR=0.196 | 1.853 | 0.043 | 0.000 | high |
| 53:37 | known | 3 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 1.867 | 0.045 | 0.118 |  |
| 35:22 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 1.409 | 0.048 | 0.236 |  |
| 35:26 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 1.409 | 0.048 | 0.236 |  |
| 39:14 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 1.409 | 0.048 | 0.236 |  |
| HLA-B_nov_296d3d00 | novel | 2 | AFR=0.005, EUR=0.000, AMR=0.995 | AFR=0.005, EUR=0.000, AMR=0.995 | 1.380 | 0.048 | 0.236 | recurrent |
| HLA-B_nov_2329b766 | novel | 2 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 1.409 | 0.048 | 0.236 | recurrent |
| HLA-B_nov_224d4cee | novel | 2 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 1.409 | 0.048 | 0.236 | recurrent |
| HLA-B_nov_2d657694 | novel | 21 | AFR=0.121, EUR=0.606, AMR=0.274 | AFR=0.134, EUR=0.589, AMR=0.278 | 1.887 | 0.052 | 0.069 | high |
| HLA-B_nov_8ecf43de | novel | 4 | AFR=0.939, EUR=0.000, AMR=0.061 | AFR=0.939, EUR=0.000, AMR=0.061 | 2.102 | 0.054 | 0.181 | high |
| HLA-B_nov_44b4b59d | novel | 18 | AFR=0.057, EUR=0.551, AMR=0.392 | AFR=0.043, EUR=0.299, AMR=0.658 | 1.721 | 0.057 | 0.000 | high |
| HLA-B_nov_d8cdc4a8 | novel | 7 | AFR=0.247, EUR=0.748, AMR=0.005 | AFR=0.325, EUR=0.669, AMR=0.007 | 1.673 | 0.060 | 0.368 | high |
| 15:54 | known | 3 | AFR=0.973, EUR=0.000, AMR=0.027 | AFR=0.973, EUR=0.000, AMR=0.027 | 1.771 | 0.067 | 0.118 |  |
| 14:05 | known | 3 | AFR=0.973, EUR=0.007, AMR=0.020 | AFR=0.973, EUR=0.007, AMR=0.020 | 1.766 | 0.069 | 0.118 |  |
| HLA-B_nov_a2393038 | novel | 29 | AFR=0.141, EUR=0.529, AMR=0.330 | AFR=0.170, EUR=0.439, AMR=0.392 | 1.525 | 0.078 | 0.007 | high |
| HLA-B_nov_0903559e | novel | 6 | AFR=0.192, EUR=0.771, AMR=0.037 | AFR=0.452, EUR=0.514, AMR=0.033 | 1.393 | 0.086 | 0.000 | high |
| 15:18 | known | 93 | AFR=0.227, EUR=0.364, AMR=0.409 | AFR=0.261, EUR=0.352, AMR=0.387 | 1.489 | 0.087 | 0.000 |  |
| HLA-B_nov_d0894e4e | novel | 25 | AFR=0.135, EUR=0.547, AMR=0.318 | AFR=0.138, EUR=0.515, AMR=0.347 | 1.413 | 0.095 | 0.250 | high |
| 40:10 | known | 7 | AFR=0.056, EUR=0.715, AMR=0.229 | AFR=0.014, EUR=0.876, AMR=0.110 | 1.292 | 0.103 | 0.000 |  |
| 18:25 | known | 2 | AFR=0.995, EUR=0.000, AMR=0.005 | AFR=0.995, EUR=0.000, AMR=0.005 | 1.341 | 0.105 | 0.320 |  |
| 81:02 | known | 3 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 1.645 | 0.106 | 0.001 |  |
| HLA-B_nov_509a423f | novel | 3 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 1.645 | 0.106 | 0.001 | high |
| 35:06 | known | 3 | AFR=0.061, EUR=0.000, AMR=0.939 | AFR=0.060, EUR=0.000, AMR=0.940 | 1.622 | 0.107 | 0.070 |  |
| 15:108 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 15:63 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 35:24 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 35:332 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 39:13 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 39:75 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 44:521 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 51:127 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| 51:237 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 |  |
| HLA-B_nov_03e38808 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_3acc735a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_32a5c027 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_3863c439 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_23ed1117 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_53b6524c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_d40006ff | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_17b10f4b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_fd89aa93 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_4e217d42 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_feb4cad1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_85bdf62f | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_dc7878dd | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_40a4247f | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_e91987c0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_c28cabe3 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_497bf2e5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_ea3116e4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_8f4eaa30 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_45b8df76 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_ea54482e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_6d208db1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_9baf3e13 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_61fbc128 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_79d30be6 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_04733de3 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_c94dd52e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_c75a2c34 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_dace3103 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_0fee0b6b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_45aa1780 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_766af2c2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.208 | flagged_artifact |
| HLA-B_nov_bd37224e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_4bf64393 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_216797a7 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_1a1089bc | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_a22c4f5b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_dd34539a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.208 | flagged_artifact |
| HLA-B_nov_4eb4ae8c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_1f65d474 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_1879609c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_f95cbd0e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_b3d0e4df | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_9abb80b8 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_1aa951a6 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_1aefebe5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_afb99d39 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.208 | flagged_artifact |
| HLA-B_nov_c30ccde0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.208 | flagged_artifact |
| HLA-B_nov_4708e0d4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_f96ccee4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_1da6c7ff | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_921de884 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_19a19987 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_d463760d | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_733cd3d4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_6d42ea59 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_043f30a2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_c0b4cb4c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_b9a1c31c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_366b17a1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_bc903ccf | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_082d2f6e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_5ae13cb0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_4921d5c2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_c646eb03 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_4fdfbb62 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_374a687a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_e3128d14 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | singleton_clean |
| HLA-B_nov_0231069a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_1728bc61 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_92d71921 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_98197dc1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.121 | flagged_artifact |
| HLA-B_nov_b03a9586 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_1652691a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| HLA-B_nov_411bb40d | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000 | AFR=0.000, EUR=0.000, AMR=1.000 | 0.873 | 0.111 | 0.639 | flagged_artifact |
| 58:11 | known | 4 | AFR=0.768, EUR=0.006, AMR=0.226 | AFR=0.788, EUR=0.005, AMR=0.207 | 1.131 | 0.116 | 0.312 |  |
| HLA-B_nov_b732cb1c | novel | 33 | AFR=0.197, EUR=0.529, AMR=0.274 | AFR=0.250, EUR=0.464, AMR=0.286 | 1.264 | 0.121 | 0.000 | high |
| HLA-B_nov_297f6a9c | novel | 4 | AFR=0.034, EUR=0.830, AMR=0.137 | AFR=0.038, EUR=0.808, AMR=0.155 | 1.078 | 0.121 | 0.204 | high |
| 35:421 | known | 2 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 1.312 | 0.123 | 0.320 |  |
| HLA-B_nov_39c31ef1 | novel | 2 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 1.312 | 0.123 | 0.320 | recurrent |
| HLA-B_nov_0d837470 | novel | 4 | AFR=0.740, EUR=0.005, AMR=0.255 | AFR=0.745, EUR=0.005, AMR=0.250 | 1.000 | 0.126 | 0.312 | high |
| 15:38 | known | 3 | AFR=0.022, EUR=0.077, AMR=0.901 | AFR=0.022, EUR=0.077, AMR=0.901 | 1.378 | 0.133 | 0.044 |  |
| HLA-B_nov_6c52313d | novel | 6 | AFR=0.028, EUR=0.727, AMR=0.244 | AFR=0.038, EUR=0.639, AMR=0.323 | 1.147 | 0.135 | 0.484 | high |
| 40:11 | known | 3 | AFR=0.110, EUR=0.000, AMR=0.890 | AFR=0.105, EUR=0.000, AMR=0.895 | 1.355 | 0.138 | 0.070 |  |
| 15:220 | known | 27 | AFR=0.211, EUR=0.542, AMR=0.247 | AFR=0.276, EUR=0.435, AMR=0.289 | 1.042 | 0.147 | 0.000 |  |
| HLA-B_nov_2ec5cc9c | novel | 2 | AFR=0.025, EUR=0.000, AMR=0.975 | AFR=0.025, EUR=0.000, AMR=0.975 | 1.266 | 0.153 | 0.236 | recurrent |
| HLA-B_nov_51a1ac75 | novel | 129 | AFR=0.387, EUR=0.297, AMR=0.317 | AFR=0.413, EUR=0.259, AMR=0.329 | 1.027 | 0.154 | 0.069 | high |
| 15:24 | known | 7 | AFR=0.026, EUR=0.385, AMR=0.589 | AFR=0.026, EUR=0.358, AMR=0.616 | 0.981 | 0.160 | 0.148 |  |
| 49:01 | known | 462 | AFR=0.337, EUR=0.326, AMR=0.337 | AFR=0.384, EUR=0.247, AMR=0.369 | 1.007 | 0.165 | 0.000 |  |
| 51:05 | known | 5 | AFR=0.146, EUR=0.754, AMR=0.100 | AFR=0.250, EUR=0.633, AMR=0.117 | 0.876 | 0.165 | 0.355 |  |
| HLA-B_nov_35db003c | novel | 67 | AFR=0.237, EUR=0.472, AMR=0.291 | AFR=0.251, EUR=0.453, AMR=0.296 | 1.025 | 0.167 | 0.003 | high |
| 47:01 | known | 28 | AFR=0.217, EUR=0.527, AMR=0.256 | AFR=0.231, EUR=0.514, AMR=0.255 | 0.896 | 0.187 | 0.005 |  |
| HLA-B_nov_6f552993 | novel | 7 | AFR=0.560, EUR=0.066, AMR=0.374 | AFR=0.558, EUR=0.052, AMR=0.390 | 0.880 | 0.196 | 0.806 | high |
| HLA-B_nov_3aa8d048 | novel | 2 | AFR=0.957, EUR=0.006, AMR=0.037 | AFR=0.958, EUR=0.005, AMR=0.037 | 1.123 | 0.196 | 0.320 | flagged_artifact |
| 39:20 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 |  |
| 44:28 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 |  |
| 45:07 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 |  |
| HLA-B_nov_15e124e2 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_e048626d | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_f760e3e5 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_5bde10d6 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_93267114 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_6a0d3300 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_91625672 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_9f7d2aaf | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_43bbcf50 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_b4bcc5dc | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_d1683cd8 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_628c88c7 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_df3a30ba | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_84a80f13 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_8dc202a4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_42312a88 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_1a014e8d | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_6f53b6c4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_ffb8b272 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_cb50fa9c | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_cf249fe3 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_60472dca | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_94facc98 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_44eb7f8e | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_42309e46 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_bb3d5087 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_3ecf65c4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_72cdf352 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_920ce8d1 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_2b5013b1 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_603eb428 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_2e8a159a | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_950f781a | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_09b874e9 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_5e41ec2f | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_1e22404d | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_c8698013 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_1202a22c | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_114217e3 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_92c9d9ec | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_565856be | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_db5506d5 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_8343f020 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_a14c370b | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_0fb8eb70 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| HLA-B_nov_b1d82175 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | singleton_clean |
| HLA-B_nov_c2686f98 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000 | AFR=1.000, EUR=0.000, AMR=0.000 | 0.823 | 0.200 | 0.711 | flagged_artifact |
| 51:09 | known | 6 | AFR=0.644, EUR=0.310, AMR=0.047 | AFR=0.643, EUR=0.313, AMR=0.044 | 0.843 | 0.205 | 0.212 |  |
| HLA-B_nov_095f95ad | novel | 17 | AFR=0.115, EUR=0.484, AMR=0.401 | AFR=0.393, EUR=0.308, AMR=0.299 | 0.784 | 0.210 | 0.000 | high |
| HLA-B_nov_e7324a00 | novel | 1 | AFR=0.010, EUR=0.000, AMR=0.990 | AFR=0.010, EUR=0.000, AMR=0.990 | 0.799 | 0.213 | 0.639 | flagged_artifact |
| HLA-B_nov_67071646 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.990 | AFR=0.000, EUR=0.010, AMR=0.990 | 0.796 | 0.221 | 0.639 | flagged_artifact |
| HLA-B_nov_813c65bc | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.990 | AFR=0.000, EUR=0.010, AMR=0.990 | 0.795 | 0.224 | 0.639 | singleton_clean |
| HLA-B_nov_6b34840c | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.990 | AFR=0.000, EUR=0.010, AMR=0.990 | 0.795 | 0.224 | 0.639 | flagged_artifact |
| 07:14 | known | 1 | AFR=0.011, EUR=0.000, AMR=0.989 | AFR=0.011, EUR=0.000, AMR=0.989 | 0.790 | 0.225 | 0.639 |  |
| HLA-B_nov_aaded50a | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.990 | AFR=0.000, EUR=0.010, AMR=0.990 | 0.794 | 0.225 | 0.639 | flagged_artifact |
| 55:04 | known | 2 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 1.088 | 0.244 | 0.014 |  |
| HLA-B_nov_fce99c63 | novel | 2 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 1.088 | 0.244 | 0.014 | recurrent |
| 15:40 | known | 2 | AFR=0.053, EUR=0.005, AMR=0.941 | AFR=0.052, EUR=0.005, AMR=0.943 | 1.073 | 0.247 | 0.236 |  |
| HLA-B_nov_e3fc05f4 | novel | 2 | AFR=0.010, EUR=0.990, AMR=0.000 | AFR=0.010, EUR=0.990, AMR=0.000 | 1.029 | 0.261 | 0.539 | flagged_artifact |
| 14:01 | known | 131 | AFR=0.326, EUR=0.306, AMR=0.368 | AFR=0.341, EUR=0.288, AMR=0.372 | 0.612 | 0.261 | 0.000 |  |
| HLA-B_nov_7dfe68fe | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_1dd91039 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_edcf0038 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_15c1d0aa | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_b277757c | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_5aff1cd6 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_33b86f8b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_c51dd187 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_90c0e76f | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_5c3c76bb | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_9f907b71 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_41152288 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_3e375010 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | singleton_clean |
| HLA-B_nov_995a3ec9 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_eee296c7 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_f2462ce0 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_dba89d4e | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_0fedbd23 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_534167c6 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_bcdf3834 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_4e13f26f | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_93fcb251 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_3dfbd43b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_00d4b295 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_617fa75b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_1f821d4b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_c9864187 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_8b4d2aff | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_72428e2d | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_d31fdb30 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.750 | 0.267 | 0.711 | flagged_artifact |
| HLA-B_nov_f8b10fb2 | novel | 3 | AFR=0.201, EUR=0.084, AMR=0.715 | AFR=0.193, EUR=0.078, AMR=0.730 | 0.349 | 0.268 | 0.070 | high |
| 07:12 | known | 2 | AFR=0.934, EUR=0.011, AMR=0.055 | AFR=0.937, EUR=0.010, AMR=0.052 | 0.990 | 0.273 | 0.320 |  |
| HLA-B_nov_1928d7e3 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_48134f61 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_b7dcf573 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_b49f3bb5 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_0dc0b6f4 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_69a8b14c | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | singleton_clean |
| HLA-B_nov_1d6684b4 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_0c2c9d42 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_541b65b4 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_7f7bcbc8 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| HLA-B_nov_1a2f8bcc | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.749 | 0.273 | 0.711 | flagged_artifact |
| 15:47 | known | 2 | AFR=0.931, EUR=0.000, AMR=0.069 | AFR=0.933, EUR=0.000, AMR=0.067 | 0.980 | 0.274 | 0.320 |  |
| HLA-B_nov_57b7bbb4 | novel | 12 | AFR=0.110, EUR=0.436, AMR=0.454 | AFR=0.130, EUR=0.374, AMR=0.496 | 0.538 | 0.276 | 0.337 | high |
| HLA-B_nov_e3452f6d | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.748 | 0.283 | 0.711 | flagged_artifact |
| HLA-B_nov_4fb34829 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.748 | 0.283 | 0.711 | flagged_artifact |
| HLA-B_nov_9fa331de | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.748 | 0.283 | 0.711 | flagged_artifact |
| HLA-B_nov_0d9622c0 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010 | AFR=0.990, EUR=0.000, AMR=0.010 | 0.747 | 0.286 | 0.711 | flagged_artifact |
| HLA-B_nov_3414ed58 | novel | 1 | AFR=0.990, EUR=0.010, AMR=0.000 | AFR=0.990, EUR=0.010, AMR=0.000 | 0.746 | 0.287 | 0.711 | flagged_artifact |
| HLA-B_nov_f8139703 | novel | 1 | AFR=0.021, EUR=0.000, AMR=0.979 | AFR=0.021, EUR=0.000, AMR=0.979 | 0.722 | 0.294 | 0.639 | flagged_artifact |
| HLA-B_nov_b4af01aa | novel | 1 | AFR=0.010, EUR=0.010, AMR=0.980 | AFR=0.010, EUR=0.010, AMR=0.980 | 0.721 | 0.294 | 0.639 | flagged_artifact |
| HLA-B_nov_83b7b3bf | novel | 1 | AFR=0.020, EUR=0.000, AMR=0.980 | AFR=0.020, EUR=0.000, AMR=0.980 | 0.725 | 0.294 | 0.639 | flagged_artifact |
| HLA-B_nov_d27b7922 | novel | 1 | AFR=0.020, EUR=0.000, AMR=0.980 | AFR=0.020, EUR=0.000, AMR=0.980 | 0.723 | 0.294 | 0.639 | singleton_clean |
| HLA-B_nov_03a90688 | novel | 1 | AFR=0.000, EUR=0.020, AMR=0.980 | AFR=0.000, EUR=0.020, AMR=0.980 | 0.719 | 0.297 | 0.639 | flagged_artifact |
| HLA-B_nov_abcae738 | novel | 1 | AFR=0.000, EUR=0.020, AMR=0.980 | AFR=0.000, EUR=0.020, AMR=0.980 | 0.719 | 0.297 | 0.639 | flagged_artifact |
| HLA-B_nov_562d19f4 | novel | 1 | AFR=0.000, EUR=0.021, AMR=0.979 | AFR=0.000, EUR=0.021, AMR=0.979 | 0.711 | 0.298 | 0.639 | flagged_artifact |
| HLA-B_nov_7716f50b | novel | 100 | AFR=0.370, EUR=0.385, AMR=0.245 | AFR=0.382, EUR=0.373, AMR=0.245 | 0.404 | 0.309 | 0.020 | high |
| HLA-B_nov_4cdba16a | novel | 12 | AFR=0.201, EUR=0.568, AMR=0.231 | AFR=0.426, EUR=0.362, AMR=0.213 | 0.402 | 0.312 | 0.000 | high |
| 15:32 | known | 5 | AFR=0.091, EUR=0.273, AMR=0.636 | AFR=0.091, EUR=0.273, AMR=0.636 | 0.472 | 0.314 | 0.002 |  |
| HLA-B_nov_93bc130b | novel | 5 | AFR=0.000, EUR=0.453, AMR=0.547 | AFR=0.000, EUR=0.389, AMR=0.611 | 0.462 | 0.318 | 0.015 | high |
| 08:12 | known | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 |  |
| HLA-B_nov_a2421080 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_ae28592c | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_95c80fc7 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | singleton_clean |
| HLA-B_nov_45078566 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_a7aafd80 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_af6a8541 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_654c769f | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_5ee74302 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | singleton_clean |
| HLA-B_nov_30e22019 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_4641b6bc | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_82fa1e91 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_14c347fc | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_01d74970 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| HLA-B_nov_a875de23 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | singleton_clean |
| HLA-B_nov_39af4334 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | singleton_clean |
| HLA-B_nov_04e3516e | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.677 | 0.323 | 0.711 | flagged_artifact |
| 51:07 | known | 14 | AFR=0.140, EUR=0.408, AMR=0.452 | AFR=0.163, EUR=0.313, AMR=0.524 | 0.346 | 0.325 | 0.247 |  |
| 07:09 | known | 4 | AFR=0.355, EUR=0.040, AMR=0.605 | AFR=0.342, EUR=0.034, AMR=0.623 | 0.370 | 0.327 | 0.503 |  |
| HLA-B_nov_90aa9c07 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.676 | 0.331 | 0.711 | singleton_clean |
| HLA-B_nov_89c83114 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.676 | 0.331 | 0.711 | flagged_artifact |
| HLA-B_nov_9c8d67db | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.676 | 0.331 | 0.711 | singleton_clean |
| HLA-B_nov_d35b12ee | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.676 | 0.331 | 0.711 | flagged_artifact |
| 18:03 | known | 5 | AFR=0.586, EUR=0.392, AMR=0.023 | AFR=0.587, EUR=0.391, AMR=0.022 | 0.396 | 0.339 | 0.386 |  |
| HLA-B_nov_52b829e0 | novel | 12 | AFR=0.356, EUR=0.517, AMR=0.127 | AFR=0.389, EUR=0.487, AMR=0.124 | 0.351 | 0.339 | 0.081 | high |
| HLA-B_nov_864a10f4 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.674 | 0.340 | 0.711 | flagged_artifact |
| HLA-B_nov_b6e49b3c | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.674 | 0.340 | 0.711 | flagged_artifact |
| HLA-B_nov_94f7382f | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.674 | 0.340 | 0.711 | flagged_artifact |
| HLA-B_nov_a93b0a2b | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.674 | 0.340 | 0.711 | flagged_artifact |
| HLA-B_nov_cda60249 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020 | AFR=0.980, EUR=0.000, AMR=0.020 | 0.674 | 0.340 | 0.711 | flagged_artifact |
| HLA-B_nov_fc728f83 | novel | 2 | AFR=0.111, EUR=0.011, AMR=0.878 | AFR=0.105, EUR=0.011, AMR=0.884 | 0.719 | 0.342 | 0.236 | recurrent |
| HLA-B_nov_4f621ce7 | novel | 1 | AFR=0.979, EUR=0.000, AMR=0.021 | AFR=0.979, EUR=0.000, AMR=0.021 | 0.668 | 0.345 | 0.711 | flagged_artifact |
| HLA-B_nov_8345fbf0 | novel | 1 | AFR=0.979, EUR=0.000, AMR=0.021 | AFR=0.979, EUR=0.000, AMR=0.021 | 0.671 | 0.345 | 0.711 | flagged_artifact |
| HLA-B_nov_1c3ee464 | novel | 1 | AFR=0.980, EUR=0.010, AMR=0.010 | AFR=0.980, EUR=0.010, AMR=0.010 | 0.669 | 0.345 | 0.711 | flagged_artifact |
| HLA-B_nov_a4ba84de | novel | 9 | AFR=0.358, EUR=0.153, AMR=0.489 | AFR=0.404, EUR=0.057, AMR=0.539 | 0.308 | 0.346 | 0.271 | high |
| HLA-B_nov_f468ba32 | novel | 1 | AFR=0.978, EUR=0.000, AMR=0.022 | AFR=0.978, EUR=0.000, AMR=0.022 | 0.664 | 0.346 | 0.711 | flagged_artifact |
| HLA-B_nov_a7c410d9 | novel | 1 | AFR=0.979, EUR=0.010, AMR=0.010 | AFR=0.979, EUR=0.010, AMR=0.010 | 0.668 | 0.346 | 0.711 | flagged_artifact |
| HLA-B_nov_78748980 | novel | 1 | AFR=0.020, EUR=0.010, AMR=0.970 | AFR=0.020, EUR=0.010, AMR=0.970 | 0.647 | 0.349 | 0.639 | singleton_clean |
| HLA-B_nov_158ad9c3 | novel | 1 | AFR=0.000, EUR=0.030, AMR=0.970 | AFR=0.000, EUR=0.030, AMR=0.970 | 0.643 | 0.349 | 0.639 | flagged_artifact |
| HLA-B_nov_f89fa8c2 | novel | 1 | AFR=0.010, EUR=0.020, AMR=0.970 | AFR=0.010, EUR=0.020, AMR=0.970 | 0.643 | 0.349 | 0.639 | flagged_artifact |
| 15:52 | known | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 |  |
| HLA-B_nov_3123b867 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | flagged_artifact |
| HLA-B_nov_41b5eac8 | novel | 1 | AFR=0.010, EUR=0.021, AMR=0.969 | AFR=0.010, EUR=0.021, AMR=0.969 | 0.638 | 0.350 | 0.639 | flagged_artifact |
| HLA-B_nov_656b54ac | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | flagged_artifact |
| HLA-B_nov_9a725575 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | singleton_clean |
| HLA-B_nov_fb966272 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | singleton_clean |
| HLA-B_nov_1ef08883 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | flagged_artifact |
| HLA-B_nov_d4bfb1c6 | novel | 1 | AFR=0.010, EUR=0.020, AMR=0.969 | AFR=0.010, EUR=0.020, AMR=0.969 | 0.641 | 0.350 | 0.639 | flagged_artifact |
| HLA-B_nov_33724e64 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | flagged_artifact |
| HLA-B_nov_e9ce0e06 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | flagged_artifact |
| HLA-B_nov_97890c08 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.605 | 0.350 | 0.711 | singleton_clean |
| HLA-B_nov_5a90260b | novel | 2 | AFR=0.000, EUR=0.929, AMR=0.071 | AFR=0.000, EUR=0.888, AMR=0.112 | 0.675 | 0.357 | 0.031 | flagged_artifact |
| 35:542N | known | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 |  |
| HLA-B_nov_f1d08636 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| HLA-B_nov_3e66ec37 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| HLA-B_nov_a5d0748a | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| HLA-B_nov_85407195 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| HLA-B_nov_fd7ba8ef | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| HLA-B_nov_8604a1ed | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030 | AFR=0.970, EUR=0.000, AMR=0.030 | 0.603 | 0.360 | 0.711 | flagged_artifact |
| 53:15 | known | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.600 | 0.361 | 0.711 |  |
| HLA-B_nov_e74d1f15 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.598 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_1a7cc0c5 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.600 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_8c0e6332 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.600 | 0.361 | 0.711 | singleton_clean |
| HLA-B_nov_415f5cea | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.598 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_af64f753 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.600 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_5e2c4f54 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.598 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_92dad439 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031 | AFR=0.969, EUR=0.000, AMR=0.031 | 0.600 | 0.361 | 0.711 | flagged_artifact |
| HLA-B_nov_442f4a56 | novel | 1 | AFR=0.967, EUR=0.000, AMR=0.033 | AFR=0.967, EUR=0.000, AMR=0.033 | 0.586 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_c177c874 | novel | 1 | AFR=0.968, EUR=0.000, AMR=0.032 | AFR=0.968, EUR=0.000, AMR=0.032 | 0.594 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_249d46ed | novel | 1 | AFR=0.968, EUR=0.000, AMR=0.032 | AFR=0.968, EUR=0.000, AMR=0.032 | 0.594 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_7d255483 | novel | 1 | AFR=0.968, EUR=0.011, AMR=0.022 | AFR=0.968, EUR=0.011, AMR=0.022 | 0.582 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_20b30691 | novel | 1 | AFR=0.968, EUR=0.000, AMR=0.032 | AFR=0.968, EUR=0.000, AMR=0.032 | 0.594 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_17512224 | novel | 1 | AFR=0.968, EUR=0.011, AMR=0.021 | AFR=0.968, EUR=0.011, AMR=0.021 | 0.585 | 0.362 | 0.711 | flagged_artifact |
| HLA-B_nov_3c2934f0 | novel | 41 | AFR=0.305, EUR=0.288, AMR=0.407 | AFR=0.305, EUR=0.289, AMR=0.406 | 0.209 | 0.365 | 0.124 | high |
| HLA-B_nov_95c48066 | novel | 1 | AFR=0.041, EUR=0.000, AMR=0.959 | AFR=0.041, EUR=0.000, AMR=0.959 | 0.576 | 0.366 | 0.639 | flagged_artifact |
| HLA-B_nov_59d5eed5 | novel | 1 | AFR=0.000, EUR=0.040, AMR=0.960 | AFR=0.000, EUR=0.040, AMR=0.960 | 0.567 | 0.367 | 0.639 | flagged_artifact |
| HLA-B_nov_c272771c | novel | 1 | AFR=0.000, EUR=0.040, AMR=0.960 | AFR=0.000, EUR=0.040, AMR=0.960 | 0.567 | 0.367 | 0.639 | flagged_artifact |
| HLA-B_nov_584448c0 | novel | 1 | AFR=0.031, EUR=0.010, AMR=0.959 | AFR=0.031, EUR=0.010, AMR=0.959 | 0.567 | 0.367 | 0.639 | singleton_clean |
| 51:285 | known | 1 | AFR=0.021, EUR=0.021, AMR=0.957 | AFR=0.021, EUR=0.021, AMR=0.957 | 0.552 | 0.368 | 0.639 |  |
| HLA-B_nov_aedf1b6f | novel | 1 | AFR=0.966, EUR=0.023, AMR=0.011 | AFR=0.966, EUR=0.023, AMR=0.011 | 0.561 | 0.368 | 0.711 | flagged_artifact |
| HLA-B_nov_4d179b4c | novel | 1 | AFR=0.010, EUR=0.031, AMR=0.958 | AFR=0.010, EUR=0.031, AMR=0.958 | 0.556 | 0.368 | 0.639 | singleton_clean |
| HLA-B_nov_726dbc56 | novel | 1 | AFR=0.021, EUR=0.021, AMR=0.957 | AFR=0.021, EUR=0.021, AMR=0.957 | 0.552 | 0.368 | 0.639 | singleton_clean |
| HLA-B_nov_b42ddaac | novel | 1 | AFR=0.021, EUR=0.021, AMR=0.959 | AFR=0.021, EUR=0.021, AMR=0.959 | 0.562 | 0.368 | 0.639 | flagged_artifact |
| 47:03 | known | 3 | AFR=0.653, EUR=0.337, AMR=0.010 | AFR=0.946, EUR=0.039, AMR=0.015 | 0.088 | 0.372 | 0.063 |  |
| 08:16 | known | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 |  |
| HLA-B_nov_cecb3b95 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 | flagged_artifact |
| HLA-B_nov_5042bbd5 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 | flagged_artifact |
| HLA-B_nov_f3b4f0b2 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 | singleton_clean |
| HLA-B_nov_3ecd044d | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 | flagged_artifact |
| HLA-B_nov_900e5341 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.533 | 0.378 | 0.711 | flagged_artifact |
| HLA-B_nov_3a379e69 | novel | 2 | AFR=0.014, EUR=0.911, AMR=0.075 | AFR=0.025, EUR=0.838, AMR=0.138 | 0.569 | 0.378 | 0.539 | recurrent |
| HLA-B_nov_dbe98bfc | novel | 2 | AFR=0.006, EUR=0.912, AMR=0.083 | AFR=0.006, EUR=0.917, AMR=0.077 | 0.576 | 0.378 | 0.313 | recurrent |
| 35:28 | known | 2 | AFR=0.077, EUR=0.066, AMR=0.857 | AFR=0.073, EUR=0.063, AMR=0.864 | 0.565 | 0.381 | 0.236 |  |
| HLA-B_nov_6de4a489 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.530 | 0.384 | 0.711 | flagged_artifact |
| HLA-B_nov_89319c90 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.530 | 0.384 | 0.711 | flagged_artifact |
| HLA-B_nov_be69758b | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.530 | 0.384 | 0.711 | flagged_artifact |
| HLA-B_nov_b4e9b4be | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040 | AFR=0.960, EUR=0.000, AMR=0.040 | 0.530 | 0.384 | 0.711 | flagged_artifact |
| HLA-B_nov_709461de | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.041 | AFR=0.959, EUR=0.000, AMR=0.041 | 0.527 | 0.386 | 0.711 | flagged_artifact |
| HLA-B_nov_e6e2087c | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.041 | AFR=0.959, EUR=0.000, AMR=0.041 | 0.524 | 0.387 | 0.711 | flagged_artifact |
| HLA-B_nov_d5349780 | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.041 | AFR=0.959, EUR=0.000, AMR=0.041 | 0.524 | 0.387 | 0.711 | singleton_clean |
| HLA-B_nov_c6839d69 | novel | 1 | AFR=0.958, EUR=0.000, AMR=0.042 | AFR=0.958, EUR=0.000, AMR=0.042 | 0.521 | 0.388 | 0.711 | singleton_clean |
| HLA-B_nov_3d91afb3 | novel | 1 | AFR=0.958, EUR=0.000, AMR=0.042 | AFR=0.958, EUR=0.000, AMR=0.042 | 0.521 | 0.388 | 0.711 | flagged_artifact |
| HLA-B_nov_c3e62ffe | novel | 1 | AFR=0.958, EUR=0.011, AMR=0.032 | AFR=0.958, EUR=0.011, AMR=0.032 | 0.511 | 0.391 | 0.711 | flagged_artifact |
| HLA-B_nov_eb1b4e4b | novel | 1 | AFR=0.957, EUR=0.000, AMR=0.043 | AFR=0.957, EUR=0.000, AMR=0.043 | 0.512 | 0.391 | 0.711 | singleton_clean |
| HLA-B_nov_1caccbe8 | novel | 1 | AFR=0.958, EUR=0.010, AMR=0.031 | AFR=0.958, EUR=0.010, AMR=0.031 | 0.514 | 0.391 | 0.711 | flagged_artifact |
| HLA-B_nov_77a21fe5 | novel | 2 | AFR=0.148, EUR=0.016, AMR=0.837 | AFR=0.148, EUR=0.016, AMR=0.836 | 0.496 | 0.391 | 0.236 | recurrent |
| HLA-B_nov_28799930 | novel | 1 | AFR=0.043, EUR=0.011, AMR=0.946 | AFR=0.043, EUR=0.011, AMR=0.946 | 0.476 | 0.392 | 0.639 | singleton_clean |
| HLA-B_nov_1f3faf7b | novel | 1 | AFR=0.043, EUR=0.011, AMR=0.946 | AFR=0.043, EUR=0.011, AMR=0.946 | 0.476 | 0.392 | 0.639 | flagged_artifact |
| HLA-B_nov_a5149631 | novel | 1 | AFR=0.031, EUR=0.021, AMR=0.948 | AFR=0.031, EUR=0.021, AMR=0.948 | 0.483 | 0.392 | 0.639 | flagged_artifact |
| HLA-B_nov_457dac52 | novel | 1 | AFR=0.000, EUR=0.051, AMR=0.949 | AFR=0.000, EUR=0.051, AMR=0.949 | 0.488 | 0.392 | 0.639 | flagged_artifact |
| HLA-B_nov_a32b50be | novel | 1 | AFR=0.954, EUR=0.023, AMR=0.023 | AFR=0.954, EUR=0.023, AMR=0.023 | 0.477 | 0.392 | 0.711 | flagged_artifact |
| HLA-B_nov_abfaeb71 | novel | 1 | AFR=0.957, EUR=0.011, AMR=0.033 | AFR=0.957, EUR=0.011, AMR=0.033 | 0.501 | 0.392 | 0.711 | singleton_clean |
| HLA-B_nov_4c450148 | novel | 1 | AFR=0.051, EUR=0.000, AMR=0.949 | AFR=0.051, EUR=0.000, AMR=0.949 | 0.503 | 0.392 | 0.639 | flagged_artifact |
| 27:06 | known | 27 | AFR=0.200, EUR=0.442, AMR=0.358 | AFR=0.764, EUR=0.124, AMR=0.112 | 0.123 | 0.394 | 0.000 |  |
| HLA-B_nov_1cad5d70 | novel | 56 | AFR=0.237, EUR=0.407, AMR=0.356 | AFR=0.236, EUR=0.392, AMR=0.373 | 0.133 | 0.400 | 0.037 | high |
| HLA-B_nov_f72ccc05 | novel | 2 | AFR=0.000, EUR=0.889, AMR=0.111 | AFR=0.000, EUR=0.889, AMR=0.111 | 0.455 | 0.402 | 0.436 | flagged_artifact |
| HLA-B_nov_6969be6c | novel | 94 | AFR=0.349, EUR=0.403, AMR=0.247 | AFR=0.400, EUR=0.346, AMR=0.255 | 0.164 | 0.403 | 0.184 | high |
| 50:02 | known | 16 | AFR=0.198, EUR=0.341, AMR=0.461 | AFR=0.228, EUR=0.261, AMR=0.510 | 0.127 | 0.407 | 0.002 |  |
| HLA-B_nov_56b81227 | novel | 20 | AFR=0.178, EUR=0.450, AMR=0.372 | AFR=0.355, EUR=0.395, AMR=0.250 | 0.145 | 0.407 | 0.000 | high |
| HLA-B_nov_34a394d7 | novel | 15 | AFR=0.170, EUR=0.498, AMR=0.332 | AFR=0.182, EUR=0.470, AMR=0.348 | 0.110 | 0.419 | 0.279 | high |
| HLA-B_nov_53f13a25 | novel | 7 | AFR=0.531, EUR=0.178, AMR=0.292 | AFR=0.554, EUR=0.158, AMR=0.288 | 0.084 | 0.425 | 0.459 | high |
| HLA-B_nov_5de4168d | novel | 2 | AFR=0.198, EUR=0.010, AMR=0.792 | AFR=0.198, EUR=0.010, AMR=0.792 | 0.282 | 0.438 | 0.509 | flagged_artifact |
| 39:15 | known | 2 | AFR=0.000, EUR=0.200, AMR=0.800 | AFR=0.000, EUR=0.200, AMR=0.800 | 0.267 | 0.444 | 0.014 |  |
| HLA-B_nov_d0b27b55 | novel | 6 | AFR=0.502, EUR=0.131, AMR=0.367 | AFR=0.509, EUR=0.124, AMR=0.367 | 0.046 | 0.444 | 0.625 | high |
| HLA-B_nov_9a4b1870 | novel | 2 | AFR=0.234, EUR=0.000, AMR=0.766 | AFR=0.234, EUR=0.000, AMR=0.766 | 0.184 | 0.454 | 0.509 | recurrent |
| HLA-B_nov_b5feb9c0 | novel | 2 | AFR=0.311, EUR=0.005, AMR=0.684 | AFR=0.311, EUR=0.005, AMR=0.684 | -0.160 | 0.493 | 0.236 | flagged_artifact |
| 15:09 | known | 7 | AFR=0.162, EUR=0.324, AMR=0.515 | AFR=0.192, EUR=0.204, AMR=0.603 | -0.099 | 0.497 | 0.016 |  |
| HLA-B_nov_7d8c1f06 | novel | 2 | AFR=0.715, EUR=0.183, AMR=0.102 | AFR=0.781, EUR=0.137, AMR=0.081 | -0.323 | 0.510 | 0.320 | recurrent |
| 07:160 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| 15:06 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 |  |
| 27:01 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| 35:41 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| 41:05 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| 41:09 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 |  |
| 44:267N | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| 56:94 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 |  |
| HLA-B_nov_6cbd58c1 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_d8ed128a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_8ae66ad7 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_cdb36446 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_dc21985f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_617ca186 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_03a5909a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_9e57f044 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 | singleton_clean |
| HLA-B_nov_c071b71a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 | flagged_artifact |
| HLA-B_nov_886a10c6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_23ce36a6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_5dd52ff0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_8b074c29 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_f6540be4 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 | flagged_artifact |
| HLA-B_nov_0e304d54 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_ef0da28e | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_f02db3fe | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_adac2bc6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_6edc2a7d | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_6312be7c | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_fd28ebd8 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 | flagged_artifact |
| HLA-B_nov_088fccbd | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_2755f44b | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_beb4db1f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_432883d0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_5e419efa | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_c567a698 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_01bdf7b0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_d2bfdcde | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_7b58694e | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_ee45d146 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_c368473f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_bef98926 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_5e724eec | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_9f3e174f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_b473fb45 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_51f79b75 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_96a3d8c0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_4ed134e4 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_c914b066 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_606353a9 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_7fe285e8 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_65dd51de | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_77dab031 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_c1c31457 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_09c249df | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_6e966e59 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_79869e0a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_92974912 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_8c14ef54 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_a05cc281 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_7373a5dd | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_2df45ea3 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_cc3f5e23 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_7da43fcf | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_09d9054d | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_41e2694d | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_17aef3f8 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_20e7ad3f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_4e8ff699 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.208 | flagged_artifact |
| HLA-B_nov_785793cf | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_d4c5cfff | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_ab3dd779 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_e05f2aaf | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_a8d3b5ae | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_639c4dd0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | singleton_clean |
| HLA-B_nov_67648ed0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | singleton_clean |
| HLA-B_nov_a0d04fdc | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_ed61f8b6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_248c0970 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_27431add | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.000 | flagged_artifact |
| HLA-B_nov_0e39d3d2 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_fb611334 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_f7a1452f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | singleton_clean |
| HLA-B_nov_5678c505 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_0628de16 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000 | AFR=0.000, EUR=1.000, AMR=0.000 | 0.465 | 0.512 | 0.705 | flagged_artifact |
| HLA-B_nov_c30cf4d2 | novel | 2 | AFR=0.173, EUR=0.123, AMR=0.704 | AFR=0.173, EUR=0.123, AMR=0.704 | -0.328 | 0.513 | 0.410 | recurrent |
| HLA-B_nov_13b01f0f | novel | 1 | AFR=0.950, EUR=0.000, AMR=0.050 | AFR=0.950, EUR=0.000, AMR=0.050 | 0.462 | 0.513 | 0.711 | singleton_clean |
| HLA-B_nov_7e151b9e | novel | 1 | AFR=0.949, EUR=0.000, AMR=0.051 | AFR=0.949, EUR=0.000, AMR=0.051 | 0.455 | 0.514 | 0.711 | singleton_clean |
| HLA-B_nov_11794504 | novel | 1 | AFR=0.949, EUR=0.000, AMR=0.051 | AFR=0.949, EUR=0.000, AMR=0.051 | 0.458 | 0.514 | 0.711 | flagged_artifact |
| HLA-B_nov_89424f29 | novel | 1 | AFR=0.948, EUR=0.000, AMR=0.052 | AFR=0.948, EUR=0.000, AMR=0.052 | 0.451 | 0.514 | 0.711 | flagged_artifact |
| HLA-B_nov_18e985ae | novel | 2 | AFR=0.000, EUR=0.340, AMR=0.660 | AFR=0.000, EUR=0.654, AMR=0.346 | -0.378 | 0.518 | 0.014 | recurrent |
| HLA-B_nov_c2365afe | novel | 1 | AFR=0.053, EUR=0.011, AMR=0.937 | AFR=0.053, EUR=0.011, AMR=0.937 | 0.409 | 0.523 | 0.639 | flagged_artifact |
| HLA-B_nov_7263fdef | novel | 1 | AFR=0.063, EUR=0.000, AMR=0.937 | AFR=0.063, EUR=0.000, AMR=0.937 | 0.417 | 0.523 | 0.639 | singleton_clean |
| HLA-B_nov_28d798e5 | novel | 1 | AFR=0.041, EUR=0.021, AMR=0.938 | AFR=0.041, EUR=0.021, AMR=0.938 | 0.412 | 0.523 | 0.639 | flagged_artifact |
| HLA-B_nov_e054ddf9 | novel | 1 | AFR=0.010, EUR=0.052, AMR=0.938 | AFR=0.010, EUR=0.052, AMR=0.938 | 0.403 | 0.524 | 0.639 | flagged_artifact |
| HLA-B_nov_180180d4 | novel | 1 | AFR=0.940, EUR=0.000, AMR=0.060 | AFR=0.940, EUR=0.000, AMR=0.060 | 0.391 | 0.525 | 0.711 | singleton_clean |
| 15:13 | known | 22 | AFR=0.408, EUR=0.386, AMR=0.207 | AFR=0.705, EUR=0.169, AMR=0.126 | -0.221 | 0.535 | 0.000 |  |
| HLA-B_nov_7c4a4f49 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.390 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_2f24cf74 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.390 | 0.535 | 0.705 | singleton_clean |
| HLA-B_nov_b947dcb4 | novel | 1 | AFR=0.010, EUR=0.990, AMR=0.000 | AFR=0.010, EUR=0.990, AMR=0.000 | 0.389 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_5d00e6a1 | novel | 1 | AFR=0.010, EUR=0.990, AMR=0.000 | AFR=0.010, EUR=0.990, AMR=0.000 | 0.389 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_c6d6cdb1 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.390 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_8ddb92d2 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.390 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_40247a71 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.390 | 0.535 | 0.705 | flagged_artifact |
| HLA-B_nov_8e8d3c8f | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.389 | 0.540 | 0.705 | flagged_artifact |
| HLA-B_nov_bf1a69a3 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.388 | 0.544 | 0.705 | flagged_artifact |
| 15:31 | known | 2 | AFR=0.485, EUR=0.000, AMR=0.515 | AFR=0.485, EUR=0.000, AMR=0.515 | -0.488 | 0.545 | 0.826 |  |
| 07:10 | known | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.388 | 0.548 | 0.705 |  |
| HLA-B_nov_318ba7f2 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.388 | 0.548 | 0.705 | flagged_artifact |
| HLA-B_nov_7d2945bd | novel | 1 | AFR=0.939, EUR=0.000, AMR=0.061 | AFR=0.939, EUR=0.000, AMR=0.061 | 0.387 | 0.550 | 0.711 | flagged_artifact |
| HLA-B_nov_d53f646f | novel | 1 | AFR=0.939, EUR=0.000, AMR=0.061 | AFR=0.939, EUR=0.000, AMR=0.061 | 0.383 | 0.551 | 0.711 | flagged_artifact |
| HLA-B_nov_4f805aaa | novel | 1 | AFR=0.000, EUR=0.989, AMR=0.011 | AFR=0.000, EUR=0.989, AMR=0.011 | 0.384 | 0.551 | 0.705 | flagged_artifact |
| HLA-B_nov_ce029c60 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010 | AFR=0.000, EUR=0.990, AMR=0.010 | 0.387 | 0.551 | 0.705 | flagged_artifact |
| HLA-B_nov_672c5e49 | novel | 1 | AFR=0.000, EUR=0.989, AMR=0.011 | AFR=0.000, EUR=0.989, AMR=0.011 | 0.382 | 0.553 | 0.705 | flagged_artifact |
| HLA-B_nov_856c84fc | novel | 1 | AFR=0.044, EUR=0.022, AMR=0.934 | AFR=0.044, EUR=0.022, AMR=0.934 | 0.382 | 0.555 | 0.639 | flagged_artifact |
| HLA-B_nov_ac5a3f46 | novel | 1 | AFR=0.938, EUR=0.000, AMR=0.062 | AFR=0.938, EUR=0.000, AMR=0.062 | 0.378 | 0.556 | 0.711 | flagged_artifact |
| HLA-B_nov_8df7809d | novel | 1 | AFR=0.939, EUR=0.010, AMR=0.051 | AFR=0.939, EUR=0.010, AMR=0.051 | 0.379 | 0.556 | 0.711 | singleton_clean |
| HLA-B_nov_fa36f7c9 | novel | 1 | AFR=0.939, EUR=0.010, AMR=0.051 | AFR=0.939, EUR=0.010, AMR=0.051 | 0.374 | 0.556 | 0.711 | flagged_artifact |
| HLA-B_nov_9455989c | novel | 1 | AFR=0.000, EUR=0.989, AMR=0.011 | AFR=0.000, EUR=0.989, AMR=0.011 | 0.379 | 0.556 | 0.705 | flagged_artifact |
| HLA-B_nov_f12cfa5b | novel | 1 | AFR=0.939, EUR=0.010, AMR=0.051 | AFR=0.939, EUR=0.010, AMR=0.051 | 0.379 | 0.556 | 0.711 | flagged_artifact |
| HLA-B_nov_28a66739 | novel | 9 | AFR=0.478, EUR=0.230, AMR=0.292 | AFR=0.478, EUR=0.230, AMR=0.292 | -0.243 | 0.559 | 0.528 | high |
| 39:31 | known | 1 | AFR=0.000, EUR=0.983, AMR=0.017 | AFR=0.000, EUR=0.983, AMR=0.017 | 0.340 | 0.560 | 0.705 |  |
| HLA-B_nov_8fec16fc | novel | 1 | AFR=0.062, EUR=0.010, AMR=0.928 | AFR=0.062, EUR=0.010, AMR=0.928 | 0.345 | 0.560 | 0.639 | flagged_artifact |
| HLA-B_nov_ad8338b7 | novel | 1 | AFR=0.934, EUR=0.044, AMR=0.022 | AFR=0.934, EUR=0.044, AMR=0.022 | 0.323 | 0.561 | 0.711 | flagged_artifact |
| HLA-B_nov_92b7a05e | novel | 1 | AFR=0.031, EUR=0.042, AMR=0.927 | AFR=0.031, EUR=0.042, AMR=0.927 | 0.322 | 0.561 | 0.639 | flagged_artifact |
| HLA-B_nov_3893b442 | novel | 1 | AFR=0.064, EUR=0.011, AMR=0.926 | AFR=0.064, EUR=0.011, AMR=0.926 | 0.328 | 0.561 | 0.639 | flagged_artifact |
| HLA-B_nov_c60316d0 | novel | 1 | AFR=0.064, EUR=0.011, AMR=0.926 | AFR=0.064, EUR=0.011, AMR=0.926 | 0.328 | 0.561 | 0.639 | flagged_artifact |
| HLA-B_nov_f46b406f | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.071 | AFR=0.929, EUR=0.000, AMR=0.071 | 0.316 | 0.562 | 0.711 | flagged_artifact |
| HLA-B_nov_b44619c9 | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.020 | AFR=0.000, EUR=0.980, AMR=0.020 | 0.316 | 0.562 | 0.705 | flagged_artifact |
| HLA-B_nov_2868079e | novel | 1 | AFR=0.054, EUR=0.022, AMR=0.925 | AFR=0.054, EUR=0.022, AMR=0.925 | 0.314 | 0.562 | 0.639 | singleton_clean |
| HLA-B_nov_36457ef8 | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.071 | AFR=0.929, EUR=0.000, AMR=0.071 | 0.316 | 0.562 | 0.711 | flagged_artifact |
| HLA-B_nov_c75cde5b | novel | 4 | AFR=0.031, EUR=0.508, AMR=0.461 | AFR=0.031, EUR=0.517, AMR=0.453 | -0.094 | 0.563 | 0.498 | high |
| 44:302 | known | 1 | AFR=0.010, EUR=0.980, AMR=0.010 | AFR=0.010, EUR=0.980, AMR=0.010 | 0.314 | 0.564 | 0.705 |  |
| HLA-B_nov_183b774a | novel | 1 | AFR=0.010, EUR=0.980, AMR=0.010 | AFR=0.010, EUR=0.980, AMR=0.010 | 0.314 | 0.564 | 0.705 | flagged_artifact |
| HLA-B_nov_ea5e695b | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.071 | AFR=0.929, EUR=0.000, AMR=0.071 | 0.311 | 0.566 | 0.711 | flagged_artifact |
| HLA-B_nov_09eb855b | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.020 | AFR=0.000, EUR=0.980, AMR=0.020 | 0.313 | 0.566 | 0.705 | flagged_artifact |
| HLA-B_nov_1a6efd18 | novel | 1 | AFR=0.000, EUR=0.979, AMR=0.021 | AFR=0.000, EUR=0.979, AMR=0.021 | 0.309 | 0.569 | 0.705 | flagged_artifact |
| HLA-B_nov_f4735665 | novel | 1 | AFR=0.000, EUR=0.979, AMR=0.021 | AFR=0.000, EUR=0.979, AMR=0.021 | 0.308 | 0.569 | 0.705 | flagged_artifact |
| HLA-B_nov_cbae2f89 | novel | 1 | AFR=0.000, EUR=0.979, AMR=0.021 | AFR=0.000, EUR=0.979, AMR=0.021 | 0.308 | 0.569 | 0.705 | flagged_artifact |
| HLA-B_nov_75cec853 | novel | 1 | AFR=0.000, EUR=0.979, AMR=0.021 | AFR=0.000, EUR=0.979, AMR=0.021 | 0.309 | 0.569 | 0.705 | flagged_artifact |
| HLA-B_nov_3995636d | novel | 1 | AFR=0.000, EUR=0.978, AMR=0.022 | AFR=0.000, EUR=0.978, AMR=0.022 | 0.299 | 0.571 | 0.705 | flagged_artifact |
| HLA-B_nov_e8b0f3b7 | novel | 1 | AFR=0.011, EUR=0.978, AMR=0.011 | AFR=0.011, EUR=0.978, AMR=0.011 | 0.297 | 0.571 | 0.705 | flagged_artifact |
| HLA-B_nov_3be52441 | novel | 1 | AFR=0.081, EUR=0.000, AMR=0.919 | AFR=0.081, EUR=0.000, AMR=0.919 | 0.294 | 0.574 | 0.639 | flagged_artifact |
| HLA-B_nov_1f850e58 | novel | 2 | AFR=0.505, EUR=0.000, AMR=0.495 | AFR=0.513, EUR=0.000, AMR=0.487 | -0.492 | 0.576 | 0.826 | recurrent |
| HLA-B_nov_0453be48 | novel | 1 | AFR=0.000, EUR=0.973, AMR=0.027 | AFR=0.000, EUR=0.973, AMR=0.027 | 0.266 | 0.580 | 0.705 | flagged_artifact |
| HLA-B_nov_950fe3aa | novel | 1 | AFR=0.052, EUR=0.031, AMR=0.917 | AFR=0.052, EUR=0.031, AMR=0.917 | 0.250 | 0.582 | 0.639 | flagged_artifact |
| HLA-B_nov_f3cd8945 | novel | 1 | AFR=0.920, EUR=0.000, AMR=0.080 | AFR=0.920, EUR=0.000, AMR=0.080 | 0.252 | 0.582 | 0.711 | flagged_artifact |
| HLA-B_nov_fde1221b | novel | 1 | AFR=0.918, EUR=0.000, AMR=0.082 | AFR=0.918, EUR=0.000, AMR=0.082 | 0.241 | 0.586 | 0.711 | flagged_artifact |
| HLA-B_nov_dadd5742 | novel | 1 | AFR=0.000, EUR=0.971, AMR=0.029 | AFR=0.000, EUR=0.971, AMR=0.029 | 0.246 | 0.586 | 0.000 | singleton_clean |
| HLA-B_nov_a915fad3 | novel | 1 | AFR=0.010, EUR=0.970, AMR=0.020 | AFR=0.010, EUR=0.970, AMR=0.020 | 0.237 | 0.586 | 0.705 | singleton_clean |
| HLA-B_nov_cb340e43 | novel | 1 | AFR=0.000, EUR=0.969, AMR=0.031 | AFR=0.000, EUR=0.969, AMR=0.031 | 0.235 | 0.590 | 0.705 | flagged_artifact |
| HLA-B_nov_c0671265 | novel | 1 | AFR=0.000, EUR=0.969, AMR=0.031 | AFR=0.000, EUR=0.969, AMR=0.031 | 0.235 | 0.590 | 0.705 | flagged_artifact |
| HLA-B_nov_f07b2d2c | novel | 1 | AFR=0.000, EUR=0.968, AMR=0.032 | AFR=0.000, EUR=0.968, AMR=0.032 | 0.230 | 0.591 | 0.705 | flagged_artifact |
| HLA-B_nov_6ca27245 | novel | 1 | AFR=0.010, EUR=0.969, AMR=0.021 | AFR=0.010, EUR=0.969, AMR=0.021 | 0.232 | 0.591 | 0.705 | flagged_artifact |
| HLA-B_nov_f50f5dda | novel | 1 | AFR=0.022, EUR=0.967, AMR=0.011 | AFR=0.022, EUR=0.967, AMR=0.011 | 0.219 | 0.598 | 0.705 | flagged_artifact |
| HLA-B_nov_5f30e7ec | novel | 1 | AFR=0.000, EUR=0.965, AMR=0.035 | AFR=0.000, EUR=0.965, AMR=0.035 | 0.206 | 0.600 | 0.705 | singleton_clean |
| HLA-B_nov_cbda84f5 | novel | 1 | AFR=0.000, EUR=0.964, AMR=0.036 | AFR=0.000, EUR=0.964, AMR=0.036 | 0.199 | 0.600 | 0.705 | flagged_artifact |
| HLA-B_nov_b0ad9c35 | novel | 1 | AFR=0.011, EUR=0.967, AMR=0.022 | AFR=0.011, EUR=0.967, AMR=0.022 | 0.214 | 0.600 | 0.705 | flagged_artifact |
| HLA-B_nov_d3ddc74c | novel | 1 | AFR=0.917, EUR=0.036, AMR=0.048 | AFR=0.917, EUR=0.036, AMR=0.048 | 0.199 | 0.602 | 0.711 | flagged_artifact |
| HLA-B_nov_633f823f | novel | 1 | AFR=0.000, EUR=0.963, AMR=0.037 | AFR=0.000, EUR=0.963, AMR=0.037 | 0.190 | 0.602 | 0.705 | flagged_artifact |
| HLA-B_nov_ff49c3b2 | novel | 1 | AFR=0.000, EUR=0.963, AMR=0.037 | AFR=0.000, EUR=0.963, AMR=0.037 | 0.190 | 0.602 | 0.705 | flagged_artifact |
| HLA-B_nov_4c946edc | novel | 1 | AFR=0.000, EUR=0.963, AMR=0.037 | AFR=0.000, EUR=0.963, AMR=0.037 | 0.190 | 0.602 | 0.000 | flagged_artifact |
| HLA-B_nov_bf76c177 | novel | 1 | AFR=0.013, EUR=0.961, AMR=0.026 | AFR=0.013, EUR=0.961, AMR=0.026 | 0.168 | 0.604 | 0.705 | singleton_clean |
| HLA-B_nov_a46e00a7 | novel | 1 | AFR=0.913, EUR=0.054, AMR=0.033 | AFR=0.913, EUR=0.054, AMR=0.033 | 0.164 | 0.606 | 0.711 | flagged_artifact |
| HLA-B_nov_90596be6 | novel | 1 | AFR=0.907, EUR=0.000, AMR=0.093 | AFR=0.907, EUR=0.000, AMR=0.093 | 0.164 | 0.607 | 0.711 | singleton_clean |
| HLA-B_nov_2a7661a6 | novel | 1 | AFR=0.010, EUR=0.960, AMR=0.030 | AFR=0.010, EUR=0.960, AMR=0.030 | 0.162 | 0.608 | 0.705 | flagged_artifact |
| HLA-B_nov_e524bf6b | novel | 1 | AFR=0.010, EUR=0.960, AMR=0.030 | AFR=0.010, EUR=0.960, AMR=0.030 | 0.162 | 0.608 | 0.705 | flagged_artifact |
| HLA-B_nov_274e126e | novel | 1 | AFR=0.076, EUR=0.022, AMR=0.902 | AFR=0.076, EUR=0.022, AMR=0.902 | 0.153 | 0.609 | 0.639 | flagged_artifact |
| HLA-B_nov_f1efcbf8 | novel | 1 | AFR=0.000, EUR=0.959, AMR=0.041 | AFR=0.000, EUR=0.959, AMR=0.041 | 0.159 | 0.609 | 0.705 | flagged_artifact |
| 08:35 | known | 1 | AFR=0.010, EUR=0.958, AMR=0.031 | AFR=0.010, EUR=0.958, AMR=0.031 | 0.152 | 0.611 | 0.705 |  |
| 44:06 | known | 1 | AFR=0.000, EUR=0.957, AMR=0.043 | AFR=0.000, EUR=0.957, AMR=0.043 | 0.143 | 0.616 | 0.705 |  |
| HLA-B_nov_df3dd115 | novel | 1 | AFR=0.000, EUR=0.954, AMR=0.046 | AFR=0.000, EUR=0.954, AMR=0.046 | 0.124 | 0.616 | 0.705 | flagged_artifact |
| 44:04 | known | 2 | AFR=0.516, EUR=0.021, AMR=0.464 | AFR=0.518, EUR=0.021, AMR=0.462 | -0.608 | 0.618 | 0.826 |  |
| 44:21 | known | 1 | AFR=0.000, EUR=0.952, AMR=0.048 | AFR=0.000, EUR=0.952, AMR=0.048 | 0.112 | 0.619 | 0.000 |  |
| HLA-B_nov_9fbe8078 | novel | 1 | AFR=0.012, EUR=0.952, AMR=0.036 | AFR=0.012, EUR=0.952, AMR=0.036 | 0.104 | 0.619 | 0.705 | flagged_artifact |
| HLA-B_nov_fd593752 | novel | 1 | AFR=0.000, EUR=0.949, AMR=0.051 | AFR=0.000, EUR=0.949, AMR=0.051 | 0.086 | 0.621 | 0.705 | flagged_artifact |
| HLA-B_nov_c3a4bc1a | novel | 1 | AFR=0.012, EUR=0.950, AMR=0.037 | AFR=0.012, EUR=0.950, AMR=0.037 | 0.090 | 0.621 | 0.705 | flagged_artifact |
| HLA-B_nov_3cdacf93 | novel | 1 | AFR=0.012, EUR=0.950, AMR=0.037 | AFR=0.012, EUR=0.950, AMR=0.037 | 0.090 | 0.621 | 0.705 | singleton_clean |
| HLA-B_nov_70d6526e | novel | 2 | AFR=0.038, EUR=0.672, AMR=0.290 | AFR=0.034, EUR=0.706, AMR=0.260 | -0.651 | 0.623 | 0.823 | flagged_artifact |
| 18:20 | known | 1 | AFR=0.000, EUR=0.947, AMR=0.053 | AFR=0.000, EUR=0.947, AMR=0.053 | 0.076 | 0.627 | 0.705 |  |
| HLA-B_nov_60277aac | novel | 1 | AFR=0.032, EUR=0.074, AMR=0.895 | AFR=0.032, EUR=0.074, AMR=0.895 | 0.076 | 0.627 | 0.639 | flagged_artifact |
| HLA-B_nov_73714420 | novel | 1 | AFR=0.000, EUR=0.946, AMR=0.054 | AFR=0.000, EUR=0.946, AMR=0.054 | 0.069 | 0.628 | 0.705 | flagged_artifact |
| HLA-B_nov_b99b4797 | novel | 57 | AFR=0.265, EUR=0.377, AMR=0.358 | AFR=0.304, EUR=0.297, AMR=0.399 | -0.462 | 0.630 | 0.001 | high |
| 35:516 | known | 1 | AFR=0.000, EUR=0.944, AMR=0.056 | AFR=0.000, EUR=0.944, AMR=0.056 | 0.050 | 0.634 | 0.705 |  |
| HLA-B_nov_198ccf04 | novel | 1 | AFR=0.000, EUR=0.944, AMR=0.056 | AFR=0.000, EUR=0.944, AMR=0.056 | 0.050 | 0.634 | 0.705 | singleton_clean |
| HLA-B_nov_59632b99 | novel | 1 | AFR=0.000, EUR=0.944, AMR=0.056 | AFR=0.000, EUR=0.944, AMR=0.056 | 0.050 | 0.634 | 0.705 | flagged_artifact |
| HLA-B_nov_99ead9b9 | novel | 1 | AFR=0.000, EUR=0.943, AMR=0.057 | AFR=0.000, EUR=0.943, AMR=0.057 | 0.046 | 0.635 | 0.705 | singleton_clean |
| HLA-B_nov_2d37ddd9 | novel | 1 | AFR=0.000, EUR=0.944, AMR=0.056 | AFR=0.000, EUR=0.944, AMR=0.056 | 0.049 | 0.635 | 0.705 | flagged_artifact |
| HLA-B_nov_bad7b9fc | novel | 1 | AFR=0.000, EUR=0.941, AMR=0.059 | AFR=0.000, EUR=0.941, AMR=0.059 | 0.031 | 0.639 | 0.705 | flagged_artifact |
| HLA-B_nov_7d4604b5 | novel | 1 | AFR=0.000, EUR=0.939, AMR=0.061 | AFR=0.000, EUR=0.939, AMR=0.061 | 0.018 | 0.642 | 0.705 | flagged_artifact |
| HLA-B_nov_a8e4136a | novel | 1 | AFR=0.000, EUR=0.939, AMR=0.061 | AFR=0.000, EUR=0.939, AMR=0.061 | 0.015 | 0.643 | 0.705 | flagged_artifact |
| 18:34 | known | 1 | AFR=0.000, EUR=0.938, AMR=0.062 | AFR=0.000, EUR=0.938, AMR=0.062 | 0.009 | 0.644 | 0.705 |  |
| HLA-B_nov_e2b93788 | novel | 1 | AFR=0.000, EUR=0.938, AMR=0.062 | AFR=0.000, EUR=0.938, AMR=0.062 | 0.004 | 0.644 | 0.705 | singleton_clean |
| HLA-B_nov_69ee4ede | novel | 1 | AFR=0.031, EUR=0.939, AMR=0.031 | AFR=0.031, EUR=0.939, AMR=0.031 | 0.004 | 0.648 | 0.705 | flagged_artifact |
| HLA-B_nov_3444d893 | novel | 1 | AFR=0.112, EUR=0.010, AMR=0.878 | AFR=0.112, EUR=0.010, AMR=0.878 | -0.005 | 0.650 | 0.639 | flagged_artifact |
| HLA-B_nov_997b815e | novel | 1 | AFR=0.087, EUR=0.033, AMR=0.880 | AFR=0.087, EUR=0.033, AMR=0.880 | -0.010 | 0.651 | 0.639 | flagged_artifact |
| HLA-B_nov_149a7a9c | novel | 1 | AFR=0.890, EUR=0.066, AMR=0.044 | AFR=0.890, EUR=0.066, AMR=0.044 | -0.008 | 0.651 | 0.711 | flagged_artifact |
| HLA-B_nov_2424728f | novel | 1 | AFR=0.043, EUR=0.074, AMR=0.883 | AFR=0.043, EUR=0.074, AMR=0.883 | -0.012 | 0.653 | 0.639 | flagged_artifact |
| 14:85N | known | 1 | AFR=0.013, EUR=0.935, AMR=0.052 | AFR=0.013, EUR=0.935, AMR=0.052 | -0.020 | 0.654 | 0.705 |  |
| 27:08 | known | 1 | AFR=0.000, EUR=0.120, AMR=0.880 | AFR=0.000, EUR=0.120, AMR=0.880 | -0.024 | 0.654 | 0.639 |  |
| HLA-B_nov_c4f7dbd5 | novel | 1 | AFR=0.126, EUR=0.000, AMR=0.874 | AFR=0.126, EUR=0.000, AMR=0.874 | -0.016 | 0.654 | 0.639 | flagged_artifact |
| HLA-B_nov_fdc9e59d | novel | 1 | AFR=0.879, EUR=0.000, AMR=0.121 | AFR=0.879, EUR=0.000, AMR=0.121 | -0.028 | 0.654 | 0.711 | flagged_artifact |
| HLA-B_nov_9a3a8c00 | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.071 | AFR=0.000, EUR=0.929, AMR=0.071 | -0.054 | 0.659 | 0.705 | flagged_artifact |
| HLA-B_nov_d661c981 | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.071 | AFR=0.000, EUR=0.929, AMR=0.071 | -0.054 | 0.659 | 0.705 | flagged_artifact |
| HLA-B_nov_67886062 | novel | 1 | AFR=0.000, EUR=0.931, AMR=0.069 | AFR=0.000, EUR=0.931, AMR=0.069 | -0.042 | 0.659 | 0.705 | singleton_clean |
| HLA-B_nov_7302f76d | novel | 1 | AFR=0.874, EUR=0.000, AMR=0.126 | AFR=0.874, EUR=0.000, AMR=0.126 | -0.062 | 0.660 | 0.711 | flagged_artifact |
| HLA-B_nov_7ae61930 | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.071 | AFR=0.000, EUR=0.929, AMR=0.071 | -0.060 | 0.660 | 0.705 | flagged_artifact |
| HLA-B_nov_d9f13efa | novel | 1 | AFR=0.041, EUR=0.082, AMR=0.876 | AFR=0.041, EUR=0.082, AMR=0.876 | -0.063 | 0.660 | 0.639 | flagged_artifact |
| HLA-B_nov_85fca092 | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.071 | AFR=0.000, EUR=0.929, AMR=0.071 | -0.060 | 0.660 | 0.000 | flagged_artifact |
| 35:55 | known | 1 | AFR=0.000, EUR=0.927, AMR=0.073 | AFR=0.000, EUR=0.927, AMR=0.073 | -0.072 | 0.661 | 0.705 |  |
| HLA-B_nov_c8196e09 | novel | 1 | AFR=0.000, EUR=0.927, AMR=0.073 | AFR=0.000, EUR=0.927, AMR=0.073 | -0.072 | 0.661 | 0.705 | flagged_artifact |
| HLA-B_nov_c2d85fff | novel | 1 | AFR=0.021, EUR=0.105, AMR=0.874 | AFR=0.021, EUR=0.105, AMR=0.874 | -0.080 | 0.665 | 0.639 | flagged_artifact |
| HLA-B_nov_8101def2 | novel | 1 | AFR=0.064, EUR=0.064, AMR=0.872 | AFR=0.064, EUR=0.064, AMR=0.872 | -0.088 | 0.667 | 0.639 | flagged_artifact |
| HLA-B_nov_d687a2c5 | novel | 1 | AFR=0.870, EUR=0.000, AMR=0.130 | AFR=0.870, EUR=0.000, AMR=0.130 | -0.086 | 0.667 | 0.711 | flagged_artifact |
| HLA-B_nov_d9ef93f6 | novel | 1 | AFR=0.138, EUR=0.000, AMR=0.862 | AFR=0.138, EUR=0.000, AMR=0.862 | -0.096 | 0.667 | 0.639 | singleton_clean |
| HLA-B_nov_30cc397d | novel | 1 | AFR=0.000, EUR=0.923, AMR=0.077 | AFR=0.000, EUR=0.923, AMR=0.077 | -0.099 | 0.669 | 0.705 | flagged_artifact |
| HLA-B_nov_0ae312f6 | novel | 1 | AFR=0.113, EUR=0.021, AMR=0.866 | AFR=0.113, EUR=0.021, AMR=0.866 | -0.097 | 0.669 | 0.639 | flagged_artifact |
| 18:05 | known | 1 | AFR=0.000, EUR=0.134, AMR=0.866 | AFR=0.000, EUR=0.134, AMR=0.866 | -0.124 | 0.672 | 0.639 |  |
| HLA-B_nov_81b45220 | novel | 1 | AFR=0.000, EUR=0.919, AMR=0.081 | AFR=0.000, EUR=0.919, AMR=0.081 | -0.131 | 0.676 | 0.705 | flagged_artifact |
| HLA-B_nov_abbf46cb | novel | 1 | AFR=0.000, EUR=0.919, AMR=0.081 | AFR=0.000, EUR=0.919, AMR=0.081 | -0.131 | 0.676 | 0.705 | flagged_artifact |
| HLA-B_nov_98c711e5 | novel | 1 | AFR=0.870, EUR=0.033, AMR=0.098 | AFR=0.870, EUR=0.033, AMR=0.098 | -0.133 | 0.676 | 0.711 | singleton_clean |
| HLA-B_nov_f718667b | novel | 1 | AFR=0.000, EUR=0.917, AMR=0.083 | AFR=0.000, EUR=0.917, AMR=0.083 | -0.145 | 0.677 | 0.000 | flagged_artifact |
| 27:17 | known | 1 | AFR=0.041, EUR=0.918, AMR=0.041 | AFR=0.041, EUR=0.918, AMR=0.041 | -0.157 | 0.684 | 0.705 |  |
| HLA-B_nov_708e0c17 | novel | 1 | AFR=0.868, EUR=0.044, AMR=0.088 | AFR=0.868, EUR=0.044, AMR=0.088 | -0.154 | 0.684 | 0.711 | flagged_artifact |
| HLA-B_nov_0786781c | novel | 1 | AFR=0.041, EUR=0.918, AMR=0.041 | AFR=0.041, EUR=0.918, AMR=0.041 | -0.157 | 0.684 | 0.705 | flagged_artifact |
| HLA-B_nov_69b8d710 | novel | 1 | AFR=0.021, EUR=0.116, AMR=0.863 | AFR=0.021, EUR=0.116, AMR=0.863 | -0.157 | 0.684 | 0.639 | flagged_artifact |
| HLA-B_nov_7dc0e55c | novel | 1 | AFR=0.000, EUR=0.914, AMR=0.086 | AFR=0.000, EUR=0.914, AMR=0.086 | -0.164 | 0.685 | 0.705 | singleton_clean |
| HLA-B_nov_8a010d86 | novel | 1 | AFR=0.014, EUR=0.915, AMR=0.070 | AFR=0.014, EUR=0.915, AMR=0.070 | -0.163 | 0.685 | 0.705 | singleton_clean |
| HLA-B_nov_2ec3ff0e | novel | 1 | AFR=0.000, EUR=0.914, AMR=0.086 | AFR=0.000, EUR=0.914, AMR=0.086 | -0.164 | 0.685 | 0.705 | flagged_artifact |
| HLA-B_nov_a97e900d | novel | 1 | AFR=0.000, EUR=0.915, AMR=0.085 | AFR=0.000, EUR=0.915, AMR=0.085 | -0.159 | 0.685 | 0.705 | flagged_artifact |
| HLA-B_nov_8738690d | novel | 1 | AFR=0.000, EUR=0.915, AMR=0.085 | AFR=0.000, EUR=0.915, AMR=0.085 | -0.159 | 0.685 | 0.705 | singleton_clean |
| HLA-B_nov_d12e5487 | novel | 1 | AFR=0.031, EUR=0.917, AMR=0.052 | AFR=0.031, EUR=0.917, AMR=0.052 | -0.161 | 0.685 | 0.705 | flagged_artifact |
| HLA-B_nov_bf61ed92 | novel | 1 | AFR=0.000, EUR=0.912, AMR=0.088 | AFR=0.000, EUR=0.912, AMR=0.088 | -0.174 | 0.686 | 0.705 | flagged_artifact |
| HLA-B_nov_e14db612 | novel | 1 | AFR=0.112, EUR=0.031, AMR=0.857 | AFR=0.112, EUR=0.031, AMR=0.857 | -0.170 | 0.686 | 0.639 | flagged_artifact |
| HLA-B_nov_fa264b30 | novel | 1 | AFR=0.000, EUR=0.912, AMR=0.088 | AFR=0.000, EUR=0.912, AMR=0.088 | -0.174 | 0.686 | 0.705 | flagged_artifact |
| HLA-B_nov_67b2eeab | novel | 1 | AFR=0.065, EUR=0.075, AMR=0.860 | AFR=0.065, EUR=0.075, AMR=0.860 | -0.181 | 0.687 | 0.639 | flagged_artifact |
| HLA-B_nov_dfac340c | novel | 1 | AFR=0.013, EUR=0.911, AMR=0.076 | AFR=0.013, EUR=0.911, AMR=0.076 | -0.192 | 0.690 | 0.705 | flagged_artifact |
| HLA-B_nov_0213941a | novel | 1 | AFR=0.000, EUR=0.911, AMR=0.089 | AFR=0.000, EUR=0.911, AMR=0.089 | -0.184 | 0.690 | 0.705 | flagged_artifact |
| HLA-B_nov_160154e6 | novel | 1 | AFR=0.012, EUR=0.912, AMR=0.075 | AFR=0.012, EUR=0.912, AMR=0.075 | -0.183 | 0.690 | 0.705 | flagged_artifact |
| HLA-B_nov_d63ce2a9 | novel | 1 | AFR=0.018, EUR=0.912, AMR=0.070 | AFR=0.018, EUR=0.912, AMR=0.070 | -0.188 | 0.690 | 0.705 | singleton_clean |
| HLA-B_nov_abeace88 | novel | 1 | AFR=0.859, EUR=0.022, AMR=0.120 | AFR=0.859, EUR=0.022, AMR=0.120 | -0.194 | 0.690 | 0.711 | flagged_artifact |
| HLA-B_nov_ff663f66 | novel | 1 | AFR=0.860, EUR=0.023, AMR=0.116 | AFR=0.860, EUR=0.023, AMR=0.116 | -0.184 | 0.690 | 0.711 | flagged_artifact |
| HLA-B_nov_5971e459 | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.092 | AFR=0.000, EUR=0.908, AMR=0.092 | -0.207 | 0.695 | 0.705 | singleton_clean |
| HLA-B_nov_4951f11c | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.092 | AFR=0.000, EUR=0.908, AMR=0.092 | -0.207 | 0.695 | 0.705 | flagged_artifact |
| HLA-B_nov_64107ce0 | novel | 1 | AFR=0.000, EUR=0.907, AMR=0.093 | AFR=0.000, EUR=0.907, AMR=0.093 | -0.215 | 0.698 | 0.705 | flagged_artifact |
| HLA-B_nov_54da8a2c | novel | 1 | AFR=0.013, EUR=0.908, AMR=0.079 | AFR=0.013, EUR=0.908, AMR=0.079 | -0.217 | 0.700 | 0.705 | flagged_artifact |
| HLA-B_nov_54cd9fa3 | novel | 1 | AFR=0.000, EUR=0.906, AMR=0.094 | AFR=0.000, EUR=0.906, AMR=0.094 | -0.218 | 0.700 | 0.705 | flagged_artifact |
| HLA-B_nov_89cf2b08 | novel | 1 | AFR=0.000, EUR=0.906, AMR=0.094 | AFR=0.000, EUR=0.906, AMR=0.094 | -0.221 | 0.700 | 0.705 | flagged_artifact |
| HLA-B_nov_2fdaced0 | novel | 1 | AFR=0.855, EUR=0.048, AMR=0.096 | AFR=0.855, EUR=0.048, AMR=0.096 | -0.248 | 0.703 | 0.711 | flagged_artifact |
| HLA-B_nov_bfb6333b | novel | 1 | AFR=0.000, EUR=0.901, AMR=0.099 | AFR=0.000, EUR=0.901, AMR=0.099 | -0.252 | 0.703 | 0.705 | flagged_artifact |
| HLA-B_nov_893848f1 | novel | 1 | AFR=0.014, EUR=0.901, AMR=0.085 | AFR=0.014, EUR=0.901, AMR=0.085 | -0.264 | 0.707 | 0.705 | flagged_artifact |
| HLA-B_nov_beff239d | novel | 1 | AFR=0.057, EUR=0.900, AMR=0.043 | AFR=0.057, EUR=0.900, AMR=0.043 | -0.289 | 0.711 | 0.705 | singleton_clean |
| HLA-B_nov_f46291ae | novel | 1 | AFR=0.839, EUR=0.000, AMR=0.161 | AFR=0.839, EUR=0.000, AMR=0.161 | -0.288 | 0.711 | 0.711 | flagged_artifact |
| HLA-B_nov_a80fc10b | novel | 1 | AFR=0.000, EUR=0.896, AMR=0.104 | AFR=0.000, EUR=0.896, AMR=0.104 | -0.289 | 0.712 | 0.705 | flagged_artifact |
| 35:30 | known | 4 | AFR=0.385, EUR=0.500, AMR=0.115 | AFR=0.762, EUR=0.010, AMR=0.228 | -0.663 | 0.718 | 0.009 |  |
| HLA-B_nov_88cce230 | novel | 1 | AFR=0.000, EUR=0.894, AMR=0.106 | AFR=0.000, EUR=0.894, AMR=0.106 | -0.306 | 0.721 | 0.705 | flagged_artifact |
| HLA-B_nov_c31c522f | novel | 1 | AFR=0.012, EUR=0.895, AMR=0.093 | AFR=0.012, EUR=0.895, AMR=0.093 | -0.305 | 0.721 | 0.705 | flagged_artifact |
| HLA-B_nov_3cf50b0e | novel | 1 | AFR=0.835, EUR=0.000, AMR=0.165 | AFR=0.835, EUR=0.000, AMR=0.165 | -0.311 | 0.722 | 0.711 | flagged_artifact |
| HLA-B_nov_b757e4ca | novel | 1 | AFR=0.030, EUR=0.896, AMR=0.075 | AFR=0.030, EUR=0.896, AMR=0.075 | -0.316 | 0.724 | 0.705 | flagged_artifact |
| HLA-B_nov_94a057b6 | novel | 1 | AFR=0.015, EUR=0.894, AMR=0.091 | AFR=0.015, EUR=0.894, AMR=0.091 | -0.318 | 0.726 | 0.705 | singleton_clean |
| HLA-B_nov_ebe632a1 | novel | 1 | AFR=0.013, EUR=0.893, AMR=0.093 | AFR=0.013, EUR=0.893, AMR=0.093 | -0.321 | 0.727 | 0.705 | flagged_artifact |
| HLA-B_nov_4f658252 | novel | 1 | AFR=0.147, EUR=0.021, AMR=0.832 | AFR=0.147, EUR=0.021, AMR=0.832 | -0.328 | 0.727 | 0.639 | singleton_clean |
| HLA-B_nov_f4e609fe | novel | 1 | AFR=0.000, EUR=0.890, AMR=0.110 | AFR=0.000, EUR=0.890, AMR=0.110 | -0.330 | 0.727 | 0.705 | flagged_artifact |
| HLA-B_nov_970e6c4a | novel | 1 | AFR=0.000, EUR=0.890, AMR=0.110 | AFR=0.000, EUR=0.890, AMR=0.110 | -0.330 | 0.727 | 0.705 | flagged_artifact |
| HLA-B_nov_2423c642 | novel | 1 | AFR=0.000, EUR=0.889, AMR=0.111 | AFR=0.000, EUR=0.889, AMR=0.111 | -0.339 | 0.728 | 0.705 | singleton_clean |
| HLA-B_nov_d73ebc0c | novel | 1 | AFR=0.165, EUR=0.010, AMR=0.825 | AFR=0.165, EUR=0.010, AMR=0.825 | -0.353 | 0.731 | 0.639 | flagged_artifact |
| HLA-B_nov_903af24a | novel | 1 | AFR=0.000, EUR=0.886, AMR=0.114 | AFR=0.000, EUR=0.886, AMR=0.114 | -0.356 | 0.731 | 0.000 | flagged_artifact |
| HLA-B_nov_470c7b6a | novel | 1 | AFR=0.010, EUR=0.156, AMR=0.833 | AFR=0.010, EUR=0.156, AMR=0.833 | -0.363 | 0.732 | 0.639 | flagged_artifact |
| HLA-B_nov_3ac46116 | novel | 1 | AFR=0.000, EUR=0.885, AMR=0.115 | AFR=0.000, EUR=0.885, AMR=0.115 | -0.369 | 0.733 | 0.705 | flagged_artifact |
| HLA-B_nov_560d3c37 | novel | 1 | AFR=0.000, EUR=0.882, AMR=0.118 | AFR=0.000, EUR=0.882, AMR=0.118 | -0.388 | 0.734 | 0.705 | flagged_artifact |
| HLA-B_nov_d48371e1 | novel | 1 | AFR=0.827, EUR=0.010, AMR=0.163 | AFR=0.827, EUR=0.010, AMR=0.163 | -0.385 | 0.734 | 0.711 | flagged_artifact |
| HLA-B_nov_f8a83d5c | novel | 1 | AFR=0.185, EUR=0.000, AMR=0.815 | AFR=0.185, EUR=0.000, AMR=0.815 | -0.392 | 0.734 | 0.639 | flagged_artifact |
| HLA-B_nov_e258e8f9 | novel | 1 | AFR=0.000, EUR=0.880, AMR=0.120 | AFR=0.000, EUR=0.880, AMR=0.120 | -0.403 | 0.735 | 0.705 | flagged_artifact |
| HLA-B_nov_9914427e | novel | 1 | AFR=0.010, EUR=0.163, AMR=0.827 | AFR=0.010, EUR=0.163, AMR=0.827 | -0.411 | 0.735 | 0.639 | flagged_artifact |
| HLA-B_nov_6209ee87 | novel | 1 | AFR=0.010, EUR=0.163, AMR=0.827 | AFR=0.010, EUR=0.163, AMR=0.827 | -0.411 | 0.735 | 0.639 | flagged_artifact |
| HLA-B_nov_9479734a | novel | 1 | AFR=0.147, EUR=0.032, AMR=0.821 | AFR=0.147, EUR=0.032, AMR=0.821 | -0.415 | 0.737 | 0.639 | flagged_artifact |
| HLA-B_nov_75c0c1a4 | novel | 1 | AFR=0.033, EUR=0.880, AMR=0.087 | AFR=0.033, EUR=0.880, AMR=0.087 | -0.427 | 0.737 | 0.705 | flagged_artifact |
| HLA-B_nov_bd6c67aa | novel | 1 | AFR=0.835, EUR=0.129, AMR=0.035 | AFR=0.835, EUR=0.129, AMR=0.035 | -0.420 | 0.737 | 0.711 | flagged_artifact |
| HLA-B_nov_009d823b | novel | 1 | AFR=0.012, EUR=0.880, AMR=0.108 | AFR=0.012, EUR=0.880, AMR=0.108 | -0.417 | 0.737 | 0.705 | flagged_artifact |
| 38:91 | known | 1 | AFR=0.014, EUR=0.877, AMR=0.110 | AFR=0.014, EUR=0.877, AMR=0.110 | -0.438 | 0.739 | 0.705 |  |
| HLA-B_nov_a0e01a96 | novel | 1 | AFR=0.830, EUR=0.057, AMR=0.114 | AFR=0.830, EUR=0.057, AMR=0.114 | -0.437 | 0.739 | 0.711 | singleton_clean |
| HLA-B_nov_2920fc70 | novel | 1 | AFR=0.000, EUR=0.875, AMR=0.125 | AFR=0.000, EUR=0.875, AMR=0.125 | -0.434 | 0.739 | 0.705 | flagged_artifact |
| HLA-B_nov_f83c8686 | novel | 1 | AFR=0.178, EUR=0.011, AMR=0.811 | AFR=0.178, EUR=0.011, AMR=0.811 | -0.441 | 0.739 | 0.639 | flagged_artifact |
| HLA-B_nov_0a3ca459 | novel | 1 | AFR=0.088, EUR=0.088, AMR=0.824 | AFR=0.088, EUR=0.088, AMR=0.824 | -0.451 | 0.744 | 0.639 | flagged_artifact |
| HLA-B_nov_ce78bc9f | novel | 1 | AFR=0.000, EUR=0.870, AMR=0.130 | AFR=0.000, EUR=0.870, AMR=0.130 | -0.467 | 0.748 | 0.705 | flagged_artifact |
| HLA-B_nov_35f7e815 | novel | 1 | AFR=0.000, EUR=0.869, AMR=0.131 | AFR=0.000, EUR=0.869, AMR=0.131 | -0.475 | 0.749 | 0.705 | flagged_artifact |
| HLA-B_nov_3f16cac8 | novel | 1 | AFR=0.000, EUR=0.868, AMR=0.132 | AFR=0.000, EUR=0.868, AMR=0.132 | -0.479 | 0.750 | 0.705 | flagged_artifact |
| HLA-B_nov_eef5bdfe | novel | 1 | AFR=0.000, EUR=0.868, AMR=0.132 | AFR=0.000, EUR=0.868, AMR=0.132 | -0.479 | 0.750 | 0.705 | flagged_artifact |
| HLA-B_nov_0087a08e | novel | 1 | AFR=0.060, EUR=0.119, AMR=0.821 | AFR=0.060, EUR=0.119, AMR=0.821 | -0.480 | 0.752 | 0.639 | singleton_clean |
| HLA-B_nov_257a6409 | novel | 1 | AFR=0.020, EUR=0.163, AMR=0.816 | AFR=0.020, EUR=0.163, AMR=0.816 | -0.492 | 0.755 | 0.639 | flagged_artifact |
| HLA-B_nov_d4b14b9a | novel | 1 | AFR=0.160, EUR=0.032, AMR=0.809 | AFR=0.160, EUR=0.032, AMR=0.809 | -0.498 | 0.756 | 0.639 | flagged_artifact |
| HLA-B_nov_ca434506 | novel | 1 | AFR=0.000, EUR=0.865, AMR=0.135 | AFR=0.000, EUR=0.865, AMR=0.135 | -0.503 | 0.758 | 0.705 | flagged_artifact |
| HLA-B_nov_80c40461 | novel | 1 | AFR=0.821, EUR=0.090, AMR=0.090 | AFR=0.821, EUR=0.090, AMR=0.090 | -0.528 | 0.764 | 0.711 | flagged_artifact |
| HLA-B_nov_10d36ffc | novel | 1 | AFR=0.167, EUR=0.031, AMR=0.802 | AFR=0.167, EUR=0.031, AMR=0.802 | -0.538 | 0.768 | 0.639 | flagged_artifact |
| HLA-B_nov_6f459aa7 | novel | 1 | AFR=0.000, EUR=0.859, AMR=0.141 | AFR=0.000, EUR=0.859, AMR=0.141 | -0.543 | 0.768 | 0.705 | singleton_clean |
| HLA-B_nov_bbbd239d | novel | 1 | AFR=0.000, EUR=0.859, AMR=0.141 | AFR=0.000, EUR=0.859, AMR=0.141 | -0.543 | 0.768 | 0.705 | flagged_artifact |
| HLA-B_nov_a023356a | novel | 1 | AFR=0.013, EUR=0.863, AMR=0.125 | AFR=0.013, EUR=0.863, AMR=0.125 | -0.535 | 0.768 | 0.705 | flagged_artifact |
| HLA-B_nov_d8da2570 | novel | 1 | AFR=0.194, EUR=0.010, AMR=0.796 | AFR=0.194, EUR=0.010, AMR=0.796 | -0.533 | 0.768 | 0.639 | singleton_clean |
| HLA-B_nov_f05a082c | novel | 1 | AFR=0.000, EUR=0.859, AMR=0.141 | AFR=0.000, EUR=0.859, AMR=0.141 | -0.542 | 0.768 | 0.705 | singleton_clean |
| HLA-B_nov_b6619d8d | novel | 1 | AFR=0.013, EUR=0.861, AMR=0.127 | AFR=0.013, EUR=0.861, AMR=0.127 | -0.547 | 0.769 | 0.705 | flagged_artifact |
| 27:09 | known | 1 | AFR=0.000, EUR=0.857, AMR=0.143 | AFR=0.000, EUR=0.857, AMR=0.143 | -0.555 | 0.770 | 0.000 |  |
| HLA-B_nov_f29be7f4 | novel | 1 | AFR=0.000, EUR=0.857, AMR=0.143 | AFR=0.000, EUR=0.857, AMR=0.143 | -0.555 | 0.770 | 0.705 | flagged_artifact |
| HLA-B_nov_7ab154bd | novel | 1 | AFR=0.000, EUR=0.857, AMR=0.143 | AFR=0.000, EUR=0.857, AMR=0.143 | -0.555 | 0.770 | 0.705 | flagged_artifact |
| HLA-B_nov_cfb772eb | novel | 1 | AFR=0.012, EUR=0.859, AMR=0.129 | AFR=0.012, EUR=0.859, AMR=0.129 | -0.560 | 0.771 | 0.705 | flagged_artifact |
| HLA-B_nov_5ca36281 | novel | 1 | AFR=0.108, EUR=0.084, AMR=0.807 | AFR=0.108, EUR=0.084, AMR=0.807 | -0.571 | 0.774 | 0.639 | singleton_clean |
| HLA-B_nov_a612b678 | novel | 1 | AFR=0.811, EUR=0.067, AMR=0.122 | AFR=0.811, EUR=0.067, AMR=0.122 | -0.575 | 0.778 | 0.711 | flagged_artifact |
| HLA-B_nov_d25ed3df | novel | 1 | AFR=0.811, EUR=0.067, AMR=0.122 | AFR=0.811, EUR=0.067, AMR=0.122 | -0.575 | 0.778 | 0.711 | flagged_artifact |
| HLA-B_nov_03cd7f7f | novel | 1 | AFR=0.000, EUR=0.853, AMR=0.147 | AFR=0.000, EUR=0.853, AMR=0.147 | -0.580 | 0.780 | 0.705 | flagged_artifact |
| HLA-B_nov_2d7e9819 | novel | 1 | AFR=0.786, EUR=0.000, AMR=0.214 | AFR=0.786, EUR=0.000, AMR=0.214 | -0.609 | 0.784 | 0.711 | flagged_artifact |
| HLA-B_nov_2783493a | novel | 1 | AFR=0.000, EUR=0.845, AMR=0.155 | AFR=0.000, EUR=0.845, AMR=0.155 | -0.633 | 0.790 | 0.705 | flagged_artifact |
| HLA-B_nov_2035a20d | novel | 1 | AFR=0.067, EUR=0.853, AMR=0.080 | AFR=0.067, EUR=0.853, AMR=0.080 | -0.639 | 0.790 | 0.705 | flagged_artifact |
| HLA-B_nov_e81ee869 | novel | 1 | AFR=0.014, EUR=0.847, AMR=0.139 | AFR=0.014, EUR=0.847, AMR=0.139 | -0.641 | 0.791 | 0.705 | singleton_clean |
| HLA-B_nov_3c32a88a | novel | 1 | AFR=0.214, EUR=0.010, AMR=0.776 | AFR=0.214, EUR=0.010, AMR=0.776 | -0.656 | 0.791 | 0.639 | flagged_artifact |
| HLA-B_nov_62fbec54 | novel | 1 | AFR=0.215, EUR=0.011, AMR=0.774 | AFR=0.215, EUR=0.011, AMR=0.774 | -0.665 | 0.791 | 0.639 | flagged_artifact |
| HLA-B_nov_05c6957c | novel | 1 | AFR=0.216, EUR=0.010, AMR=0.773 | AFR=0.216, EUR=0.010, AMR=0.773 | -0.670 | 0.794 | 0.639 | flagged_artifact |
| HLA-B_nov_dd140996 | novel | 1 | AFR=0.012, EUR=0.841, AMR=0.146 | AFR=0.012, EUR=0.841, AMR=0.146 | -0.678 | 0.796 | 0.705 | flagged_artifact |
| HLA-B_nov_f0ecc128 | novel | 1 | AFR=0.770, EUR=0.000, AMR=0.230 | AFR=0.770, EUR=0.000, AMR=0.230 | -0.698 | 0.801 | 0.711 | singleton_clean |
| HLA-B_nov_41499f6d | novel | 1 | AFR=0.000, EUR=0.835, AMR=0.165 | AFR=0.000, EUR=0.835, AMR=0.165 | -0.697 | 0.801 | 0.705 | flagged_artifact |
| HLA-B_nov_a5806dd7 | novel | 1 | AFR=0.205, EUR=0.023, AMR=0.773 | AFR=0.205, EUR=0.023, AMR=0.773 | -0.704 | 0.802 | 0.639 | flagged_artifact |
| HLA-B_nov_e4455130 | novel | 19 | AFR=0.353, EUR=0.293, AMR=0.354 | AFR=0.539, EUR=0.096, AMR=0.365 | -0.895 | 0.809 | 0.000 | high |
| HLA-B_nov_75d764e9 | novel | 1 | AFR=0.011, EUR=0.832, AMR=0.158 | AFR=0.011, EUR=0.832, AMR=0.158 | -0.741 | 0.812 | 0.705 | flagged_artifact |
| HLA-B_nov_b1992d3b | novel | 1 | AFR=0.041, EUR=0.837, AMR=0.122 | AFR=0.041, EUR=0.837, AMR=0.122 | -0.744 | 0.813 | 0.208 | flagged_artifact |
| HLA-B_nov_e7abe413 | novel | 1 | AFR=0.760, EUR=0.010, AMR=0.229 | AFR=0.760, EUR=0.010, AMR=0.229 | -0.783 | 0.815 | 0.711 | flagged_artifact |
| HLA-B_nov_c734a2f0 | novel | 1 | AFR=0.259, EUR=0.000, AMR=0.741 | AFR=0.259, EUR=0.000, AMR=0.741 | -0.818 | 0.818 | 0.639 | singleton_clean |
| HLA-B_nov_5a5869cb | novel | 1 | AFR=0.753, EUR=0.010, AMR=0.237 | AFR=0.753, EUR=0.010, AMR=0.237 | -0.826 | 0.821 | 0.711 | flagged_artifact |
| HLA-B_nov_38e03e62 | novel | 1 | AFR=0.013, EUR=0.818, AMR=0.169 | AFR=0.013, EUR=0.818, AMR=0.169 | -0.832 | 0.821 | 0.705 | flagged_artifact |
| HLA-B_nov_8b9103a4 | novel | 1 | AFR=0.000, EUR=0.816, AMR=0.184 | AFR=0.000, EUR=0.816, AMR=0.184 | -0.819 | 0.821 | 0.705 | singleton_clean |
| HLA-B_nov_451add33 | novel | 1 | AFR=0.000, EUR=0.812, AMR=0.188 | AFR=0.000, EUR=0.812, AMR=0.188 | -0.843 | 0.821 | 0.000 | flagged_artifact |
| HLA-B_nov_c4efb5c4 | novel | 1 | AFR=0.231, EUR=0.022, AMR=0.747 | AFR=0.231, EUR=0.022, AMR=0.747 | -0.852 | 0.821 | 0.639 | flagged_artifact |
| HLA-B_nov_7409cfcb | novel | 1 | AFR=0.014, EUR=0.819, AMR=0.167 | AFR=0.014, EUR=0.819, AMR=0.167 | -0.826 | 0.821 | 0.705 | flagged_artifact |
| HLA-B_nov_9080fd27 | novel | 1 | AFR=0.018, EUR=0.818, AMR=0.164 | AFR=0.018, EUR=0.818, AMR=0.164 | -0.841 | 0.821 | 0.705 | flagged_artifact |
| HLA-B_nov_e609b2a7 | novel | 1 | AFR=0.018, EUR=0.818, AMR=0.164 | AFR=0.018, EUR=0.818, AMR=0.164 | -0.841 | 0.821 | 0.705 | flagged_artifact |
| HLA-B_nov_f924859c | novel | 1 | AFR=0.103, EUR=0.126, AMR=0.770 | AFR=0.103, EUR=0.126, AMR=0.770 | -0.862 | 0.822 | 0.639 | flagged_artifact |
| HLA-B_nov_fb8dafaa | novel | 1 | AFR=0.270, EUR=0.000, AMR=0.730 | AFR=0.270, EUR=0.000, AMR=0.730 | -0.876 | 0.823 | 0.639 | flagged_artifact |
| HLA-B_nov_f62b3ed6 | novel | 1 | AFR=0.021, EUR=0.223, AMR=0.755 | AFR=0.021, EUR=0.223, AMR=0.755 | -0.907 | 0.826 | 0.639 | flagged_artifact |
| HLA-B_nov_587a4cc2 | novel | 1 | AFR=0.260, EUR=0.010, AMR=0.729 | AFR=0.260, EUR=0.010, AMR=0.729 | -0.916 | 0.827 | 0.639 | flagged_artifact |
| HLA-B_nov_9ee0b8f8 | novel | 2 | AFR=0.491, EUR=0.479, AMR=0.031 | AFR=0.503, EUR=0.467, AMR=0.031 | -0.979 | 0.828 | 0.856 | flagged_artifact |
| 18:13 | known | 1 | AFR=0.176, EUR=0.071, AMR=0.753 | AFR=0.176, EUR=0.071, AMR=0.753 | -0.927 | 0.829 | 0.639 |  |
| HLA-B_nov_3afdde33 | novel | 1 | AFR=0.283, EUR=0.000, AMR=0.717 | AFR=0.283, EUR=0.000, AMR=0.717 | -0.941 | 0.832 | 0.639 | flagged_artifact |
| HLA-B_nov_6e4db735 | novel | 1 | AFR=0.156, EUR=0.091, AMR=0.753 | AFR=0.156, EUR=0.091, AMR=0.753 | -0.955 | 0.833 | 0.639 | flagged_artifact |
| HLA-B_nov_10aabe32 | novel | 1 | AFR=0.186, EUR=0.070, AMR=0.744 | AFR=0.186, EUR=0.070, AMR=0.744 | -0.982 | 0.834 | 0.639 | flagged_artifact |
| HLA-B_nov_5600bf31 | novel | 1 | AFR=0.734, EUR=0.043, AMR=0.223 | AFR=0.734, EUR=0.043, AMR=0.223 | -1.024 | 0.838 | 0.711 | flagged_artifact |
| HLA-B_nov_d30b9a84 | novel | 1 | AFR=0.145, EUR=0.108, AMR=0.747 | AFR=0.145, EUR=0.108, AMR=0.747 | -1.018 | 0.838 | 0.639 | flagged_artifact |
| HLA-B_nov_8cc7819b | novel | 1 | AFR=0.735, EUR=0.048, AMR=0.217 | AFR=0.735, EUR=0.048, AMR=0.217 | -1.034 | 0.841 | 0.711 | flagged_artifact |
| HLA-B_nov_048d0d89 | novel | 1 | AFR=0.714, EUR=0.022, AMR=0.264 | AFR=0.714, EUR=0.022, AMR=0.264 | -1.069 | 0.843 | 0.711 | flagged_artifact |
| HLA-B_nov_736a8631 | novel | 1 | AFR=0.043, EUR=0.217, AMR=0.739 | AFR=0.043, EUR=0.217, AMR=0.739 | -1.049 | 0.843 | 0.639 | singleton_clean |
| HLA-B_nov_5f4e9a68 | novel | 4 | AFR=0.262, EUR=0.506, AMR=0.232 | AFR=0.267, EUR=0.516, AMR=0.217 | -1.072 | 0.844 | 0.842 | high |
| HLA-B_nov_21f78371 | novel | 1 | AFR=0.293, EUR=0.011, AMR=0.696 | AFR=0.293, EUR=0.011, AMR=0.696 | -1.088 | 0.845 | 0.639 | flagged_artifact |
| HLA-B_nov_6de6f05c | novel | 1 | AFR=0.010, EUR=0.776, AMR=0.214 | AFR=0.010, EUR=0.776, AMR=0.214 | -1.093 | 0.845 | 0.705 | flagged_artifact |
| HLA-B_nov_1745f037 | novel | 1 | AFR=0.031, EUR=0.781, AMR=0.188 | AFR=0.031, EUR=0.781, AMR=0.188 | -1.104 | 0.848 | 0.121 | flagged_artifact |
| HLA-B_nov_88f0f010 | novel | 1 | AFR=0.110, EUR=0.154, AMR=0.736 | AFR=0.110, EUR=0.154, AMR=0.736 | -1.118 | 0.849 | 0.639 | singleton_clean |
| HLA-B_nov_aa4fcd4b | novel | 1 | AFR=0.224, EUR=0.059, AMR=0.718 | AFR=0.224, EUR=0.059, AMR=0.718 | -1.123 | 0.849 | 0.639 | singleton_clean |
| HLA-B_nov_65b3cfad | novel | 1 | AFR=0.000, EUR=0.286, AMR=0.714 | AFR=0.000, EUR=0.286, AMR=0.714 | -1.110 | 0.849 | 0.208 | flagged_artifact |
| HLA-B_nov_fb8ee4e5 | novel | 1 | AFR=0.684, EUR=0.000, AMR=0.316 | AFR=0.684, EUR=0.000, AMR=0.316 | -1.128 | 0.850 | 0.711 | singleton_clean |
| HLA-B_nov_a4143085 | novel | 1 | AFR=0.684, EUR=0.000, AMR=0.316 | AFR=0.684, EUR=0.000, AMR=0.316 | -1.128 | 0.850 | 0.711 | flagged_artifact |
| HLA-B_nov_f8b8518b | novel | 1 | AFR=0.740, EUR=0.143, AMR=0.117 | AFR=0.740, EUR=0.143, AMR=0.117 | -1.137 | 0.850 | 0.711 | flagged_artifact |
| 08:09 | known | 1 | AFR=0.337, EUR=0.000, AMR=0.663 | AFR=0.337, EUR=0.000, AMR=0.663 | -1.186 | 0.854 | 0.639 |  |
| HLA-B_nov_88b5e881 | novel | 1 | AFR=0.316, EUR=0.010, AMR=0.673 | AFR=0.316, EUR=0.010, AMR=0.673 | -1.187 | 0.854 | 0.639 | flagged_artifact |
| HLA-B_nov_d5cd981a | novel | 2 | AFR=0.150, EUR=0.633, AMR=0.217 | AFR=0.250, EUR=0.500, AMR=0.250 | -1.067 | 0.859 | 0.000 | recurrent |
| 55:97N | known | 1 | AFR=0.094, EUR=0.188, AMR=0.719 | AFR=0.094, EUR=0.188, AMR=0.719 | -1.243 | 0.864 | 0.639 |  |
| HLA-B_nov_29d1f9f7 | novel | 1 | AFR=0.650, EUR=0.000, AMR=0.350 | AFR=0.650, EUR=0.000, AMR=0.350 | -1.262 | 0.865 | 0.711 | flagged_artifact |
| HLA-B_nov_978b3f67 | novel | 1 | AFR=0.640, EUR=0.000, AMR=0.360 | AFR=0.640, EUR=0.000, AMR=0.360 | -1.298 | 0.869 | 0.711 | flagged_artifact |
| HLA-B_nov_95cbb061 | novel | 1 | AFR=0.067, EUR=0.756, AMR=0.178 | AFR=0.067, EUR=0.756, AMR=0.178 | -1.332 | 0.873 | 0.208 | flagged_artifact |
| HLA-B_nov_2516c759 | novel | 1 | AFR=0.107, EUR=0.190, AMR=0.702 | AFR=0.107, EUR=0.190, AMR=0.702 | -1.369 | 0.876 | 0.639 | flagged_artifact |
| HLA-B_nov_52cc5c74 | novel | 1 | AFR=0.000, EUR=0.719, AMR=0.281 | AFR=0.000, EUR=0.719, AMR=0.281 | -1.373 | 0.876 | 0.208 | flagged_artifact |
| HLA-B_nov_97d4c44c | novel | 1 | AFR=0.000, EUR=0.714, AMR=0.286 | AFR=0.000, EUR=0.714, AMR=0.286 | -1.395 | 0.878 | 0.705 | flagged_artifact |
| HLA-B_nov_128bc54d | novel | 1 | AFR=0.707, EUR=0.160, AMR=0.133 | AFR=0.707, EUR=0.160, AMR=0.133 | -1.390 | 0.878 | 0.711 | flagged_artifact |
| HLA-B_nov_5e92b57a | novel | 1 | AFR=0.072, EUR=0.747, AMR=0.181 | AFR=0.072, EUR=0.747, AMR=0.181 | -1.398 | 0.878 | 0.705 | flagged_artifact |
| HLA-B_nov_4fa3d83a | novel | 1 | AFR=0.083, EUR=0.750, AMR=0.167 | AFR=0.083, EUR=0.750, AMR=0.167 | -1.392 | 0.878 | 0.000 | singleton_clean |
| HLA-B_nov_3e25c4fd | novel | 1 | AFR=0.053, EUR=0.737, AMR=0.211 | AFR=0.053, EUR=0.737, AMR=0.211 | -1.427 | 0.880 | 0.705 | singleton_clean |
| HLA-B_nov_3b05e11b | novel | 1 | AFR=0.693, EUR=0.114, AMR=0.193 | AFR=0.693, EUR=0.114, AMR=0.193 | -1.435 | 0.880 | 0.711 | singleton_clean |
| HLA-B_nov_fc9e0671 | novel | 1 | AFR=0.596, EUR=0.000, AMR=0.404 | AFR=0.596, EUR=0.000, AMR=0.404 | -1.427 | 0.880 | 0.711 | flagged_artifact |
| HLA-B_nov_feaa531e | novel | 1 | AFR=0.025, EUR=0.725, AMR=0.250 | AFR=0.025, EUR=0.725, AMR=0.250 | -1.425 | 0.880 | 0.121 | flagged_artifact |
| HLA-B_nov_1db89393 | novel | 1 | AFR=0.384, EUR=0.010, AMR=0.606 | AFR=0.384, EUR=0.010, AMR=0.606 | -1.440 | 0.881 | 0.639 | flagged_artifact |
| HLA-B_nov_723315d8 | novel | 1 | AFR=0.420, EUR=0.000, AMR=0.580 | AFR=0.420, EUR=0.000, AMR=0.580 | -1.449 | 0.883 | 0.639 | flagged_artifact |
| HLA-B_nov_80d1d7b0 | novel | 1 | AFR=0.385, EUR=0.010, AMR=0.604 | AFR=0.385, EUR=0.010, AMR=0.604 | -1.447 | 0.883 | 0.639 | flagged_artifact |
| HLA-B_nov_608ea697 | novel | 1 | AFR=0.426, EUR=0.000, AMR=0.574 | AFR=0.426, EUR=0.000, AMR=0.574 | -1.460 | 0.884 | 0.639 | flagged_artifact |
| HLA-B_nov_e8dacf34 | novel | 1 | AFR=0.697, EUR=0.184, AMR=0.118 | AFR=0.697, EUR=0.184, AMR=0.118 | -1.462 | 0.885 | 0.711 | flagged_artifact |
| HLA-B_nov_32ecffa1 | novel | 1 | AFR=0.646, EUR=0.042, AMR=0.312 | AFR=0.646, EUR=0.042, AMR=0.312 | -1.474 | 0.886 | 0.208 | flagged_artifact |
| HLA-B_nov_4d67a1e9 | novel | 1 | AFR=0.029, EUR=0.714, AMR=0.257 | AFR=0.029, EUR=0.714, AMR=0.257 | -1.495 | 0.891 | 0.208 | flagged_artifact |
| HLA-B_nov_bf18667a | novel | 1 | AFR=0.449, EUR=0.000, AMR=0.551 | AFR=0.449, EUR=0.000, AMR=0.551 | -1.502 | 0.892 | 0.639 | flagged_artifact |
| 35:49 | known | 1 | AFR=0.455, EUR=0.000, AMR=0.545 | AFR=0.455, EUR=0.000, AMR=0.545 | -1.510 | 0.894 | 0.639 |  |
| HLA-B_nov_4a3fccf1 | novel | 1 | AFR=0.690, EUR=0.190, AMR=0.119 | AFR=0.690, EUR=0.190, AMR=0.119 | -1.514 | 0.894 | 0.711 | flagged_artifact |
| HLA-B_nov_d082bee6 | novel | 1 | AFR=0.085, EUR=0.732, AMR=0.183 | AFR=0.085, EUR=0.732, AMR=0.183 | -1.519 | 0.895 | 0.705 | flagged_artifact |
| HLA-B_nov_3974baa8 | novel | 1 | AFR=0.680, EUR=0.120, AMR=0.200 | AFR=0.680, EUR=0.120, AMR=0.200 | -1.532 | 0.896 | 0.711 | singleton_clean |
| HLA-B_nov_4379b659 | novel | 1 | AFR=0.100, EUR=0.222, AMR=0.678 | AFR=0.100, EUR=0.222, AMR=0.678 | -1.539 | 0.897 | 0.639 | flagged_artifact |
| HLA-B_nov_8e6f67ac | novel | 1 | AFR=0.510, EUR=0.000, AMR=0.490 | AFR=0.510, EUR=0.000, AMR=0.490 | -1.542 | 0.900 | 0.711 | flagged_artifact |
| HLA-B_nov_6fe23a07 | novel | 1 | AFR=0.027, EUR=0.703, AMR=0.270 | AFR=0.027, EUR=0.703, AMR=0.270 | -1.551 | 0.901 | 0.208 | flagged_artifact |
| HLA-B_nov_fd264f91 | novel | 1 | AFR=0.391, EUR=0.022, AMR=0.587 | AFR=0.391, EUR=0.022, AMR=0.587 | -1.561 | 0.902 | 0.639 | singleton_clean |
| HLA-B_nov_124cb7d8 | novel | 1 | AFR=0.667, EUR=0.100, AMR=0.233 | AFR=0.667, EUR=0.100, AMR=0.233 | -1.570 | 0.902 | 0.121 | flagged_artifact |
| 38:09 | known | 2 | AFR=0.470, EUR=0.464, AMR=0.066 | AFR=0.511, EUR=0.423, AMR=0.066 | -1.193 | 0.904 | 0.856 |  |
| HLA-B_nov_121201a0 | novel | 1 | AFR=0.128, EUR=0.198, AMR=0.674 | AFR=0.128, EUR=0.198, AMR=0.674 | -1.582 | 0.906 | 0.639 | flagged_artifact |
| HLA-B_nov_2ef23553 | novel | 1 | AFR=0.000, EUR=0.667, AMR=0.333 | AFR=0.000, EUR=0.667, AMR=0.333 | -1.605 | 0.910 | 0.000 | singleton_clean |
| HLA-B_nov_3654bb38 | novel | 1 | AFR=0.016, EUR=0.677, AMR=0.306 | AFR=0.016, EUR=0.677, AMR=0.306 | -1.632 | 0.912 | 0.705 | flagged_artifact |
| HLA-B_nov_530c73c3 | novel | 1 | AFR=0.047, EUR=0.698, AMR=0.256 | AFR=0.047, EUR=0.698, AMR=0.256 | -1.642 | 0.916 | 0.208 | flagged_artifact |
| 40:139 | known | 1 | AFR=0.667, EUR=0.210, AMR=0.123 | AFR=0.667, EUR=0.210, AMR=0.123 | -1.690 | 0.917 | 0.711 |  |
| HLA-B_nov_d9fdfca4 | novel | 1 | AFR=0.667, EUR=0.210, AMR=0.123 | AFR=0.667, EUR=0.210, AMR=0.123 | -1.690 | 0.917 | 0.711 | singleton_clean |
| 51:14 | known | 1 | AFR=0.091, EUR=0.261, AMR=0.648 | AFR=0.091, EUR=0.261, AMR=0.648 | -1.729 | 0.921 | 0.639 |  |
| HLA-B_nov_983d840a | novel | 1 | AFR=0.033, EUR=0.667, AMR=0.300 | AFR=0.033, EUR=0.667, AMR=0.300 | -1.753 | 0.922 | 0.208 | flagged_artifact |
| HLA-B_nov_68a9e512 | novel | 1 | AFR=0.021, EUR=0.379, AMR=0.600 | AFR=0.021, EUR=0.379, AMR=0.600 | -1.756 | 0.923 | 0.639 | flagged_artifact |
| HLA-B_nov_42efe47f | novel | 1 | AFR=0.000, EUR=0.439, AMR=0.561 | AFR=0.000, EUR=0.439, AMR=0.561 | -1.776 | 0.926 | 0.639 | singleton_clean |
| HLA-B_nov_42b8e4ab | novel | 1 | AFR=0.091, EUR=0.273, AMR=0.636 | AFR=0.091, EUR=0.273, AMR=0.636 | -1.801 | 0.927 | 0.208 | flagged_artifact |
| HLA-B_nov_f02e6ac0 | novel | 1 | AFR=0.000, EUR=0.607, AMR=0.393 | AFR=0.000, EUR=0.607, AMR=0.393 | -1.796 | 0.927 | 0.208 | flagged_artifact |
| HLA-B_nov_3a91837f | novel | 1 | AFR=0.161, EUR=0.195, AMR=0.644 | AFR=0.161, EUR=0.195, AMR=0.644 | -1.815 | 0.930 | 0.639 | flagged_artifact |
| HLA-B_nov_e368bf67 | novel | 1 | AFR=0.646, EUR=0.183, AMR=0.171 | AFR=0.646, EUR=0.183, AMR=0.171 | -1.841 | 0.931 | 0.711 | singleton_clean |
| HLA-B_nov_7cdbfe7c | novel | 1 | AFR=0.130, EUR=0.688, AMR=0.182 | AFR=0.130, EUR=0.688, AMR=0.182 | -1.870 | 0.934 | 0.705 | flagged_artifact |
| 53:05 | known | 1 | AFR=0.000, EUR=0.556, AMR=0.444 | AFR=0.000, EUR=0.556, AMR=0.444 | -1.882 | 0.936 | 0.000 |  |
| HLA-B_nov_a98e7bc5 | novel | 1 | AFR=0.000, EUR=0.500, AMR=0.500 | AFR=0.000, EUR=0.500, AMR=0.500 | -1.884 | 0.937 | 0.208 | flagged_artifact |
| 15:34 | known | 2 | AFR=0.104, EUR=0.517, AMR=0.379 | AFR=0.102, EUR=0.527, AMR=0.371 | -1.355 | 0.940 | 0.823 |  |
| HLA-B_nov_6d041ba3 | novel | 1 | AFR=0.429, EUR=0.055, AMR=0.516 | AFR=0.429, EUR=0.055, AMR=0.516 | -1.920 | 0.943 | 0.639 | flagged_artifact |
| HLA-B_nov_72abb4d5 | novel | 1 | AFR=0.548, EUR=0.065, AMR=0.387 | AFR=0.548, EUR=0.065, AMR=0.387 | -1.931 | 0.946 | 0.711 | singleton_clean |
| HLA-B_nov_b6c41ddd | novel | 1 | AFR=0.471, EUR=0.059, AMR=0.471 | AFR=0.471, EUR=0.059, AMR=0.471 | -1.985 | 0.946 | 0.121 | flagged_artifact |
| HLA-B_nov_cb221017 | novel | 1 | AFR=0.010, EUR=0.495, AMR=0.495 | AFR=0.010, EUR=0.495, AMR=0.495 | -1.961 | 0.946 | 0.639 | flagged_artifact |
| HLA-B_nov_5947d4af | novel | 1 | AFR=0.617, EUR=0.148, AMR=0.235 | AFR=0.617, EUR=0.148, AMR=0.235 | -1.985 | 0.946 | 0.711 | singleton_clean |
| HLA-B_nov_8c596552 | novel | 1 | AFR=0.163, EUR=0.225, AMR=0.612 | AFR=0.163, EUR=0.225, AMR=0.612 | -2.050 | 0.948 | 0.639 | flagged_artifact |
| HLA-B_nov_4e2d3568 | novel | 1 | AFR=0.617, EUR=0.198, AMR=0.185 | AFR=0.617, EUR=0.198, AMR=0.185 | -2.059 | 0.949 | 0.711 | flagged_artifact |
| HLA-B_nov_34f5e896 | novel | 1 | AFR=0.617, EUR=0.198, AMR=0.185 | AFR=0.617, EUR=0.198, AMR=0.185 | -2.059 | 0.949 | 0.711 | flagged_artifact |
| 27:10 | known | 1 | AFR=0.068, EUR=0.364, AMR=0.568 | AFR=0.068, EUR=0.364, AMR=0.568 | -2.101 | 0.951 | 0.639 |  |
| HLA-B_nov_9ad05e51 | novel | 1 | AFR=0.031, EUR=0.552, AMR=0.417 | AFR=0.031, EUR=0.552, AMR=0.417 | -2.101 | 0.951 | 0.705 | flagged_artifact |
| HLA-B_nov_915a65b7 | novel | 1 | AFR=0.091, EUR=0.636, AMR=0.273 | AFR=0.091, EUR=0.636, AMR=0.273 | -2.115 | 0.954 | 0.000 | flagged_artifact |
| HLA-B_nov_f1cd89ab | novel | 1 | AFR=0.071, EUR=0.607, AMR=0.321 | AFR=0.071, EUR=0.607, AMR=0.321 | -2.178 | 0.960 | 0.121 | flagged_artifact |
| HLA-B_nov_1c72f545 | novel | 1 | AFR=0.105, EUR=0.632, AMR=0.263 | AFR=0.105, EUR=0.632, AMR=0.263 | -2.186 | 0.962 | 0.000 | flagged_artifact |
| HLA-B_nov_13b4e364 | novel | 1 | AFR=0.043, EUR=0.478, AMR=0.478 | AFR=0.043, EUR=0.478, AMR=0.478 | -2.209 | 0.963 | 0.639 | flagged_artifact |
| HLA-B_nov_3d00e477 | novel | 1 | AFR=0.043, EUR=0.478, AMR=0.478 | AFR=0.043, EUR=0.478, AMR=0.478 | -2.209 | 0.963 | 0.639 | flagged_artifact |
| 27:14 | known | 1 | AFR=0.120, EUR=0.304, AMR=0.576 | AFR=0.120, EUR=0.304, AMR=0.576 | -2.236 | 0.964 | 0.639 |  |
| HLA-B_nov_ae080ac9 | novel | 1 | AFR=0.067, EUR=0.400, AMR=0.533 | AFR=0.067, EUR=0.400, AMR=0.533 | -2.231 | 0.964 | 0.639 | flagged_artifact |
| HLA-B_nov_5daeb7bc | novel | 1 | AFR=0.062, EUR=0.438, AMR=0.500 | AFR=0.062, EUR=0.438, AMR=0.500 | -2.301 | 0.965 | 0.121 | flagged_artifact |
| HLA-B_nov_ceb2f866 | novel | 1 | AFR=0.585, EUR=0.185, AMR=0.231 | AFR=0.585, EUR=0.185, AMR=0.231 | -2.266 | 0.965 | 0.711 | flagged_artifact |
| HLA-B_nov_e6962c16 | novel | 1 | AFR=0.062, EUR=0.500, AMR=0.438 | AFR=0.062, EUR=0.500, AMR=0.438 | -2.365 | 0.968 | 0.121 | flagged_artifact |
| HLA-B_nov_ae76889a | novel | 1 | AFR=0.076, EUR=0.435, AMR=0.489 | AFR=0.076, EUR=0.435, AMR=0.489 | -2.409 | 0.970 | 0.639 | flagged_artifact |
| HLA-B_nov_bf505b36 | novel | 1 | AFR=0.156, EUR=0.578, AMR=0.266 | AFR=0.156, EUR=0.578, AMR=0.266 | -2.635 | 0.974 | 0.705 | flagged_artifact |
| HLA-B_nov_087e3604 | novel | 1 | AFR=0.535, EUR=0.310, AMR=0.155 | AFR=0.535, EUR=0.310, AMR=0.155 | -2.626 | 0.974 | 0.711 | flagged_artifact |
| HLA-B_nov_4fac5a8f | novel | 1 | AFR=0.556, EUR=0.222, AMR=0.222 | AFR=0.556, EUR=0.222, AMR=0.222 | -2.516 | 0.974 | 0.711 | flagged_artifact |
| 35:205 | known | 1 | AFR=0.447, EUR=0.158, AMR=0.395 | AFR=0.447, EUR=0.158, AMR=0.395 | -2.719 | 0.977 | 0.121 |  |
| HLA-B_nov_1db2a33a | novel | 1 | AFR=0.179, EUR=0.308, AMR=0.513 | AFR=0.179, EUR=0.308, AMR=0.513 | -2.761 | 0.977 | 0.639 | flagged_artifact |
| HLA-B_nov_81116cfa | novel | 1 | AFR=0.382, EUR=0.176, AMR=0.441 | AFR=0.382, EUR=0.176, AMR=0.441 | -2.844 | 0.977 | 0.639 | flagged_artifact |
| HLA-B_nov_0f47e2cd | novel | 1 | AFR=0.382, EUR=0.176, AMR=0.441 | AFR=0.382, EUR=0.176, AMR=0.441 | -2.844 | 0.977 | 0.639 | flagged_artifact |
| HLA-B_nov_d39c0c1c | novel | 1 | AFR=0.115, EUR=0.423, AMR=0.462 | AFR=0.115, EUR=0.423, AMR=0.462 | -2.713 | 0.977 | 0.208 | flagged_artifact |
| HLA-B_nov_6b37665d | novel | 1 | AFR=0.507, EUR=0.338, AMR=0.155 | AFR=0.507, EUR=0.338, AMR=0.155 | -2.788 | 0.977 | 0.711 | flagged_artifact |
| 56:02 | known | 1 | AFR=0.250, EUR=0.250, AMR=0.500 | AFR=0.250, EUR=0.250, AMR=0.500 | -2.884 | 0.978 | 0.208 |  |
| HLA-B_nov_6fe77099 | novel | 1 | AFR=0.250, EUR=0.250, AMR=0.500 | AFR=0.250, EUR=0.250, AMR=0.500 | -2.884 | 0.978 | 0.208 | singleton_clean |
| HLA-B_nov_8877fa30 | novel | 1 | AFR=0.486, EUR=0.333, AMR=0.181 | AFR=0.486, EUR=0.333, AMR=0.181 | -2.982 | 0.981 | 0.711 | flagged_artifact |
| HLA-B_nov_164a08f8 | novel | 1 | AFR=0.159, EUR=0.492, AMR=0.349 | AFR=0.159, EUR=0.492, AMR=0.349 | -3.030 | 0.981 | 0.705 | flagged_artifact |
| HLA-B_nov_b7d1f5da | novel | 1 | AFR=0.448, EUR=0.388, AMR=0.164 | AFR=0.448, EUR=0.388, AMR=0.164 | -3.091 | 0.985 | 0.711 | flagged_artifact |
| HLA-B_nov_7d093379 | novel | 1 | AFR=0.300, EUR=0.500, AMR=0.200 | AFR=0.300, EUR=0.500, AMR=0.200 | -3.224 | 0.986 | 0.705 | flagged_artifact |
| HLA-B_nov_e3a8b2e1 | novel | 1 | AFR=0.205, EUR=0.474, AMR=0.321 | AFR=0.205, EUR=0.474, AMR=0.321 | -3.337 | 0.986 | 0.705 | flagged_artifact |
| HLA-B_nov_877845eb | novel | 1 | AFR=0.430, EUR=0.405, AMR=0.165 | AFR=0.430, EUR=0.405, AMR=0.165 | -3.143 | 0.986 | 0.711 | flagged_artifact |
| HLA-B_nov_cc11ebef | novel | 1 | AFR=0.317, EUR=0.467, AMR=0.217 | AFR=0.317, EUR=0.467, AMR=0.217 | -3.451 | 0.989 | 0.705 | flagged_artifact |
| HLA-B_nov_c012ebe3 | novel | 1 | AFR=0.432, EUR=0.333, AMR=0.235 | AFR=0.432, EUR=0.333, AMR=0.235 | -3.439 | 0.989 | 0.711 | flagged_artifact |
| HLA-B_nov_6c14ecd9 | novel | 1 | AFR=0.225, EUR=0.451, AMR=0.324 | AFR=0.225, EUR=0.451, AMR=0.324 | -3.525 | 0.990 | 0.705 | flagged_artifact |
| HLA-B_nov_c0d643d8 | novel | 1 | AFR=0.355, EUR=0.421, AMR=0.224 | AFR=0.355, EUR=0.421, AMR=0.224 | -3.626 | 0.994 | 0.705 | flagged_artifact |
| HLA-B_nov_a58e54e7 | novel | 1 | AFR=0.391, EUR=0.312, AMR=0.297 | AFR=0.391, EUR=0.312, AMR=0.297 | -3.749 | 0.998 | 0.711 | flagged_artifact |
| HLA-B_nov_bfa755fa | novel | 1 | AFR=0.288, EUR=0.441, AMR=0.271 | AFR=0.288, EUR=0.441, AMR=0.271 | -3.748 | 0.998 | 0.705 | flagged_artifact |
| HLA-B_nov_eaa1f00a | novel | 1 | AFR=0.333, EUR=0.420, AMR=0.246 | AFR=0.333, EUR=0.420, AMR=0.246 | -3.773 | 0.999 | 0.705 | flagged_artifact |
| HLA-B_nov_35865647 | novel | 1 | AFR=0.317, EUR=0.333, AMR=0.349 | AFR=0.317, EUR=0.333, AMR=0.349 | -4.003 | 0.999 | 0.208 | flagged_artifact |
| 15:113 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 |  |
| 15:12 | known | 14 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 |  |
| 15:75 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 |  |
| 27:25 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 |  |
| 39:04 | known | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.014 |  |
| 40:229 | known | 3 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.001 |  |
| 40:23 | known | 4 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 |  |
| 48:03 | known | 3 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.001 |  |
| 48:04 | known | 3 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 |  |
| 50:04 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 |  |
| 51:361N | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 |  |
| 52:04 | known | 5 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 |  |
| 52:07 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 |  |
| 55:07 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 |  |
| 56:03 | known | 3 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.001 |  |
| 56:05 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 |  |
| 57:61 | known | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 |  |
| 59:01 | known | 5 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 |  |
| HLA-B_nov_82017cd4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_a6fe07fd | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_99e3a911 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_4cd85998 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_8c43039b | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_7a05cd54 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_e30cee8e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_dddc4b92 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_b8962c6b | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_fc3bd25b | novel | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.225 | recurrent |
| HLA-B_nov_a625005e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_9d08ae6f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ede46040 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_af8b859b | novel | 13 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 | high |
| HLA-B_nov_bb5af146 | novel | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.014 | recurrent |
| HLA-B_nov_8ec2f62f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_296e2bcf | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_5e16dd02 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_63aa1ae4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_8976677d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_f98f02f4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_3f8f404d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_6a5cc862 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_bfd7e2ea | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_db3f6940 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_2e4f87ea | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_a9522e91 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_4c3ab01e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_eaffd835 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_15670a99 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_9a1bd411 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_aec7bc6f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_fef19996 | novel | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.014 | recurrent |
| HLA-B_nov_24d4ad14 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_942ad49d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_d8a464bb | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_d6153940 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ef674346 | novel | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.014 | recurrent |
| HLA-B_nov_429fc1b7 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_d9bb7679 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_d857b794 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_4b101839 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_9cc4537a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_515802b9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_1e84ffb4 | novel | 2 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.014 | recurrent |
| HLA-B_nov_80e14459 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_bcbd6621 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_c0369ebd | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_313a73ca | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_263d6ace | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_7227aa04 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_d62c1684 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_61b1817d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ab7a1970 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_cca9288a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_6a4e0ed5 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_8c6fe73c | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ab839558 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_d2dafa0e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_9a779eb0 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_337f87e9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_cbd8313f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_57717e15 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_b0dbc94d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ab01cadc | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_beb4d124 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_026afb40 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_296c34c7 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_0e2a36dc | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_8f012d91 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_6a8fd8e2 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_088dcf0d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_3a5db847 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_f6e168be | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_2e864900 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_67256714 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_a0acfc2b | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_843fa887 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_17621e32 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_6afdcf32 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_6f6320ab | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_6e536944 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_9f0c5e12 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_163a5aba | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_1904755d | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_a9c9570f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_79c3a6b5 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_6afb2f6a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_d284e465 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_6beaadc6 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_6f36e254 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.000 | flagged_artifact |
| HLA-B_nov_02ffb226 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_009a400c | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_960090af | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_f06cf905 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_ab73e528 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_1c5e6303 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_3e88c405 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_40b6e9d2 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_faaa0ea8 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_ead802fb | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_bc1468c9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_95a6cace | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_e0f96dec | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_44d61bb0 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_be498733 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_0c948668 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_38683399 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_57aa1303 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_f010c0ce | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_97003b12 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_f65a63ce | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_fcebb155 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_c17448dc | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_b57de30e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_027777b6 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_348960dc | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_f18f7081 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_5d23f27a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_bdf27ef9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_694dd122 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_be336e5a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_021f48ad | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_964805a7 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_c7440247 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_e5aa62e4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_b8f68360 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_d02f84ab | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_387c65d5 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_b85da312 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_a6e4e751 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_e4a3c143 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | singleton_clean |
| HLA-B_nov_a3f62ec9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_49b0379f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_b0c6c4a8 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_663f338e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_8e0e0943 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_4145800f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_6df32fb8 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_af30735e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_7d6911e6 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_b5602c21 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_41eb4588 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_746946c1 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_9e1df2cd | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_24a8297a | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_63c119eb | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_5ea37b02 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_dc9657cb | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_f374e762 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | singleton_clean |
| HLA-B_nov_843c5aa2 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_4a14d990 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_091f0220 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_b196ca06 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_16472b9f | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_2284d26e | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_d1996f44 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_fd209ef1 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.208 | flagged_artifact |
| HLA-B_nov_67f61a63 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_85ee38b4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |
| HLA-B_nov_001d03b1 | novel | 1 | AFR=nan, EUR=nan, AMR=nan | AFR=0.000, EUR=0.000, AMR=0.000 | NA | NA | 0.121 | flagged_artifact |

## Tetrahedron-space alleles

| allele | kind | n_carriers | centroid | bootstrap p | chi2 p |
|---|---|---|---|---|---|
| 07:02 | known | 1359 | AFR=0.317, EUR=0.376, AMR=0.231, EAS=0.076 | 0.000 | 0.000 |
| 07:05 | known | 118 | AFR=0.011, EUR=0.485, AMR=0.084, EAS=0.419 | 0.000 | 0.000 |
| 07:06 | known | 133 | AFR=0.396, EUR=0.139, AMR=0.267, EAS=0.198 | 0.000 | 0.000 |
| 08:01 | known | 1029 | AFR=0.246, EUR=0.474, AMR=0.226, EAS=0.053 | 0.000 | 0.000 |
| 13:01 | known | 115 | AFR=0.008, EUR=0.022, AMR=0.027, EAS=0.943 | 0.000 | 0.000 |
| 13:02 | known | 380 | AFR=0.144, EUR=0.449, AMR=0.225, EAS=0.182 | 0.000 | 0.000 |
| 14:02 | known | 693 | AFR=0.199, EUR=0.445, AMR=0.316, EAS=0.041 | 0.000 | 0.000 |
| 14:03 | known | 16 | AFR=0.909, EUR=0.005, AMR=0.078, EAS=0.008 | 0.000 | 0.000 |
| 15:01 | known | 566 | AFR=0.135, EUR=0.428, AMR=0.233, EAS=0.204 | 0.000 | 0.000 |
| 15:02 | known | 200 | AFR=0.011, EUR=0.029, AMR=0.027, EAS=0.933 | 0.000 | 0.000 |
| 15:03 | known | 477 | AFR=0.752, EUR=0.033, AMR=0.195, EAS=0.019 | 0.000 | 0.000 |
| 15:10 | known | 266 | AFR=0.741, EUR=0.052, AMR=0.179, EAS=0.027 | 0.000 | 0.000 |
| 15:11 | known | 33 | AFR=0.003, EUR=0.034, AMR=0.033, EAS=0.931 | 0.000 | 0.000 |
| 15:12 | known | 14 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.000 | 0.000 |
| 15:13 | known | 22 | AFR=0.218, EUR=0.056, AMR=0.043, EAS=0.683 | 0.000 | 0.000 |
| 15:15 | known | 35 | AFR=0.030, EUR=0.009, AMR=0.945, EAS=0.016 | 0.000 | 0.000 |
| 15:16 | known | 167 | AFR=0.707, EUR=0.066, AMR=0.208, EAS=0.018 | 0.000 | 0.000 |
| 15:17 | known | 158 | AFR=0.197, EUR=0.509, AMR=0.253, EAS=0.042 | 0.000 | 0.000 |
| 15:18 | known | 93 | AFR=0.154, EUR=0.218, AMR=0.253, EAS=0.376 | 0.000 | 0.000 |
| 15:21 | known | 37 | AFR=0.014, EUR=0.041, AMR=0.039, EAS=0.906 | 0.000 | 0.000 |
| 15:25 | known | 52 | AFR=0.025, EUR=0.001, AMR=0.088, EAS=0.886 | 0.000 | 0.000 |
| 15:27 | known | 17 | AFR=0.000, EUR=0.000, AMR=0.001, EAS=0.999 | 0.000 | 0.000 |
| 15:30 | known | 17 | AFR=0.022, EUR=0.006, AMR=0.960, EAS=0.011 | 0.000 | 0.000 |
| 15:32 | known | 5 | AFR=0.008, EUR=0.025, AMR=0.058, EAS=0.908 | 0.000 | 0.002 |
| 15:35 | known | 46 | AFR=0.012, EUR=0.043, AMR=0.054, EAS=0.891 | 0.000 | 0.000 |
| 18:01 | known | 753 | AFR=0.240, EUR=0.429, AMR=0.253, EAS=0.079 | 0.000 | 0.000 |
| 18:02 | known | 12 | AFR=0.000, EUR=0.002, AMR=0.000, EAS=0.998 | 0.000 | 0.000 |
| 27:02 | known | 61 | AFR=0.070, EUR=0.611, AMR=0.284, EAS=0.034 | 0.000 | 0.000 |
| 27:03 | known | 31 | AFR=0.731, EUR=0.078, AMR=0.175, EAS=0.017 | 0.000 | 0.000 |
| 27:04 | known | 42 | AFR=0.000, EUR=0.027, AMR=0.001, EAS=0.972 | 0.000 | 0.000 |
| 27:05 | known | 374 | AFR=0.162, EUR=0.432, AMR=0.355, EAS=0.051 | 0.000 | 0.000 |
| 27:06 | known | 27 | AFR=0.073, EUR=0.015, AMR=0.019, EAS=0.892 | 0.000 | 0.000 |
| 27:07 | known | 35 | AFR=0.012, EUR=0.670, AMR=0.199, EAS=0.118 | 0.000 | 0.000 |
| 35:01 | known | 1294 | AFR=0.320, EUR=0.282, AMR=0.306, EAS=0.092 | 0.000 | 0.000 |
| 35:02 | known | 299 | AFR=0.036, EUR=0.684, AMR=0.231, EAS=0.049 | 0.000 | 0.000 |
| 35:03 | known | 445 | AFR=0.051, EUR=0.521, AMR=0.300, EAS=0.128 | 0.000 | 0.000 |
| 35:04 | known | 14 | AFR=0.066, EUR=0.005, AMR=0.922, EAS=0.007 | 0.000 | 0.000 |
| 35:05 | known | 106 | AFR=0.036, EUR=0.039, AMR=0.088, EAS=0.837 | 0.000 | 0.000 |
| 35:08 | known | 115 | AFR=0.041, EUR=0.620, AMR=0.292, EAS=0.047 | 0.000 | 0.000 |
| 35:12 | known | 92 | AFR=0.039, EUR=0.020, AMR=0.926, EAS=0.014 | 0.000 | 0.000 |
| 35:14 | known | 11 | AFR=0.008, EUR=0.004, AMR=0.968, EAS=0.020 | 0.000 | 0.000 |
| 35:16 | known | 7 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.000 | 0.000 |
| 35:17 | known | 37 | AFR=0.008, EUR=0.030, AMR=0.956, EAS=0.006 | 0.000 | 0.000 |
| 35:20 | known | 10 | AFR=0.183, EUR=0.014, AMR=0.795, EAS=0.008 | 0.000 | 0.000 |
| 35:43 | known | 27 | AFR=0.010, EUR=0.015, AMR=0.966, EAS=0.009 | 0.000 | 0.000 |
| 37:01 | known | 208 | AFR=0.139, EUR=0.516, AMR=0.212, EAS=0.133 | 0.000 | 0.000 |
| 38:01 | known | 508 | AFR=0.043, EUR=0.667, AMR=0.225, EAS=0.066 | 0.000 | 0.000 |
| 38:02 | known | 173 | AFR=0.006, EUR=0.029, AMR=0.036, EAS=0.929 | 0.000 | 0.000 |
| 39:01 | known | 174 | AFR=0.110, EUR=0.306, AMR=0.315, EAS=0.268 | 0.000 | 0.000 |
| 39:02 | known | 20 | AFR=0.047, EUR=0.015, AMR=0.880, EAS=0.059 | 0.000 | 0.000 |
| 39:05 | known | 127 | AFR=0.044, EUR=0.024, AMR=0.903, EAS=0.029 | 0.000 | 0.000 |
| 39:06 | known | 146 | AFR=0.039, EUR=0.141, AMR=0.798, EAS=0.023 | 0.000 | 0.000 |
| 39:08 | known | 19 | AFR=0.026, EUR=0.011, AMR=0.957, EAS=0.006 | 0.000 | 0.000 |
| 39:09 | known | 12 | AFR=0.002, EUR=0.022, AMR=0.241, EAS=0.735 | 0.000 | 0.000 |
| 39:10 | known | 64 | AFR=0.754, EUR=0.074, AMR=0.147, EAS=0.024 | 0.000 | 0.000 |
| 39:24 | known | 12 | AFR=0.087, EUR=0.866, AMR=0.040, EAS=0.007 | 0.000 | 0.000 |
| 40:01 | known | 678 | AFR=0.116, EUR=0.335, AMR=0.147, EAS=0.402 | 0.000 | 0.000 |
| 40:02 | known | 383 | AFR=0.064, EUR=0.197, AMR=0.508, EAS=0.231 | 0.000 | 0.000 |
| 40:04 | known | 22 | AFR=0.057, EUR=0.009, AMR=0.920, EAS=0.013 | 0.000 | 0.000 |
| 40:05 | known | 29 | AFR=0.058, EUR=0.008, AMR=0.924, EAS=0.009 | 0.000 | 0.000 |
| 40:06 | known | 265 | AFR=0.091, EUR=0.132, AMR=0.098, EAS=0.679 | 0.000 | 0.000 |
| 40:08 | known | 12 | AFR=0.011, EUR=0.038, AMR=0.947, EAS=0.004 | 0.000 | 0.000 |
| 41:01 | known | 196 | AFR=0.140, EUR=0.496, AMR=0.325, EAS=0.038 | 0.000 | 0.000 |
| 41:02 | known | 133 | AFR=0.218, EUR=0.562, AMR=0.183, EAS=0.037 | 0.000 | 0.000 |
| 42:01 | known | 365 | AFR=0.780, EUR=0.029, AMR=0.172, EAS=0.019 | 0.000 | 0.000 |
| 42:02 | known | 54 | AFR=0.532, EUR=0.069, AMR=0.368, EAS=0.031 | 0.000 | 0.000 |
| 44:02 | known | 736 | AFR=0.128, EUR=0.532, AMR=0.276, EAS=0.064 | 0.000 | 0.000 |
| 44:03 | known | 1139 | AFR=0.294, EUR=0.270, AMR=0.334, EAS=0.102 | 0.000 | 0.000 |
| 44:05 | known | 32 | AFR=0.062, EUR=0.573, AMR=0.320, EAS=0.045 | 0.000 | 0.000 |
| 44:10 | known | 13 | AFR=0.942, EUR=0.005, AMR=0.045, EAS=0.009 | 0.000 | 0.000 |
| 44:27 | known | 20 | AFR=0.045, EUR=0.793, AMR=0.098, EAS=0.065 | 0.000 | 0.000 |
| 45:01 | known | 421 | AFR=0.619, EUR=0.082, AMR=0.267, EAS=0.032 | 0.000 | 0.000 |
| 46:01 | known | 216 | AFR=0.003, EUR=0.025, AMR=0.022, EAS=0.951 | 0.000 | 0.000 |
| 48:01 | known | 167 | AFR=0.038, EUR=0.084, AMR=0.485, EAS=0.393 | 0.000 | 0.000 |
| 48:02 | known | 9 | AFR=0.051, EUR=0.030, AMR=0.867, EAS=0.051 | 0.000 | 0.000 |
| 49:01 | known | 462 | AFR=0.330, EUR=0.306, AMR=0.325, EAS=0.038 | 0.000 | 0.000 |
| 50:01 | known | 328 | AFR=0.200, EUR=0.485, AMR=0.249, EAS=0.066 | 0.000 | 0.000 |
| 51:01 | known | 923 | AFR=0.164, EUR=0.340, AMR=0.298, EAS=0.198 | 0.000 | 0.000 |
| 51:02 | known | 66 | AFR=0.078, EUR=0.012, AMR=0.489, EAS=0.421 | 0.000 | 0.000 |
| 51:06 | known | 31 | AFR=0.010, EUR=0.005, AMR=0.063, EAS=0.922 | 0.000 | 0.000 |
| 51:08 | known | 25 | AFR=0.042, EUR=0.658, AMR=0.208, EAS=0.092 | 0.000 | 0.000 |
| 52:01 | known | 595 | AFR=0.202, EUR=0.355, AMR=0.278, EAS=0.166 | 0.000 | 0.000 |
| 53:01 | known | 825 | AFR=0.722, EUR=0.064, AMR=0.191, EAS=0.023 | 0.000 | 0.000 |
| 54:01 | known | 84 | AFR=0.002, EUR=0.035, AMR=0.012, EAS=0.951 | 0.000 | 0.000 |
| 55:01 | known | 233 | AFR=0.122, EUR=0.552, AMR=0.254, EAS=0.072 | 0.000 | 0.000 |
| 55:02 | known | 49 | AFR=0.000, EUR=0.015, AMR=0.002, EAS=0.983 | 0.000 | 0.000 |
| 56:04 | known | 15 | AFR=0.000, EUR=0.019, AMR=0.030, EAS=0.951 | 0.000 | 0.000 |
| 57:01 | known | 497 | AFR=0.112, EUR=0.587, AMR=0.204, EAS=0.096 | 0.000 | 0.000 |
| 57:02 | known | 30 | AFR=0.632, EUR=0.152, AMR=0.193, EAS=0.023 | 0.000 | 0.000 |
| 57:03 | known | 270 | AFR=0.715, EUR=0.040, AMR=0.224, EAS=0.020 | 0.000 | 0.000 |
| 58:01 | known | 704 | AFR=0.365, EUR=0.156, AMR=0.156, EAS=0.323 | 0.000 | 0.000 |
| 58:02 | known | 257 | AFR=0.807, EUR=0.012, AMR=0.164, EAS=0.017 | 0.000 | 0.000 |
| 59:01 | known | 5 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.000 | 0.000 |
| 67:01 | known | 9 | AFR=0.005, EUR=0.009, AMR=0.064, EAS=0.922 | 0.000 | 0.000 |
| 78:01 | known | 66 | AFR=0.756, EUR=0.031, AMR=0.197, EAS=0.016 | 0.000 | 0.000 |
| 81:01 | known | 142 | AFR=0.785, EUR=0.029, AMR=0.163, EAS=0.023 | 0.000 | 0.000 |
| 82:01 | known | 17 | AFR=0.832, EUR=0.007, AMR=0.155, EAS=0.006 | 0.000 | 0.000 |
| HLA-B_nov_2a9ca252 | novel | 52 | AFR=0.656, EUR=0.017, AMR=0.317, EAS=0.010 | 0.000 | 0.000 |
| HLA-B_nov_d15615ce | novel | 38 | AFR=0.055, EUR=0.015, AMR=0.922, EAS=0.008 | 0.000 | 0.000 |
| HLA-B_nov_5b19fb36 | novel | 40 | AFR=0.000, EUR=0.111, AMR=0.113, EAS=0.777 | 0.000 | 0.000 |
| HLA-B_nov_44b4b59d | novel | 18 | AFR=0.005, EUR=0.041, AMR=0.073, EAS=0.880 | 0.000 | 0.000 |
| HLA-B_nov_a3e9f87a | novel | 11 | AFR=0.000, EUR=0.000, AMR=0.001, EAS=0.999 | 0.000 | 0.000 |
| HLA-B_nov_74d648cd | novel | 39 | AFR=0.005, EUR=0.237, AMR=0.512, EAS=0.246 | 0.000 | 0.000 |
| HLA-B_nov_7752bb66 | novel | 60 | AFR=0.066, EUR=0.561, AMR=0.307, EAS=0.066 | 0.000 | 0.000 |
| HLA-B_nov_53f8de08 | novel | 23 | AFR=0.757, EUR=0.115, AMR=0.112, EAS=0.015 | 0.000 | 0.000 |
| HLA-B_nov_2a7f9a38 | novel | 45 | AFR=0.715, EUR=0.015, AMR=0.254, EAS=0.016 | 0.000 | 0.000 |
| HLA-B_nov_2670c4c2 | novel | 57 | AFR=0.056, EUR=0.419, AMR=0.466, EAS=0.058 | 0.000 | 0.000 |
| HLA-B_nov_ab1351b4 | novel | 77 | AFR=0.748, EUR=0.048, AMR=0.193, EAS=0.011 | 0.000 | 0.000 |
| HLA-B_nov_7571b4cf | novel | 58 | AFR=0.030, EUR=0.275, AMR=0.079, EAS=0.615 | 0.000 | 0.000 |
| HLA-B_nov_def7a880 | novel | 96 | AFR=0.667, EUR=0.070, AMR=0.239, EAS=0.025 | 0.000 | 0.000 |
| HLA-B_nov_56b81227 | novel | 20 | AFR=0.045, EUR=0.034, AMR=0.025, EAS=0.897 | 0.000 | 0.000 |
| HLA-B_nov_abf2f906 | novel | 26 | AFR=0.771, EUR=0.136, AMR=0.075, EAS=0.018 | 0.000 | 0.000 |
| HLA-B_nov_cd0f9284 | novel | 67 | AFR=0.414, EUR=0.144, AMR=0.178, EAS=0.265 | 0.000 | 0.003 |
| HLA-B_nov_18232e17 | novel | 35 | AFR=0.865, EUR=0.008, AMR=0.119, EAS=0.008 | 0.000 | 0.000 |
| HLA-B_nov_860e6c9d | novel | 7 | AFR=0.006, EUR=0.050, AMR=0.932, EAS=0.012 | 0.000 | 0.000 |
| HLA-B_nov_af51adc0 | novel | 23 | AFR=0.016, EUR=0.056, AMR=0.344, EAS=0.584 | 0.000 | 0.000 |
| HLA-B_nov_118d7054 | novel | 44 | AFR=0.103, EUR=0.613, AMR=0.198, EAS=0.085 | 0.000 | 0.000 |
| HLA-B_nov_bca4b2b9 | novel | 28 | AFR=0.707, EUR=0.012, AMR=0.272, EAS=0.008 | 0.000 | 0.000 |
| HLA-B_nov_8a2e0836 | novel | 28 | AFR=0.000, EUR=0.069, AMR=0.006, EAS=0.925 | 0.000 | 0.000 |
| HLA-B_nov_1d3870fc | novel | 44 | AFR=0.836, EUR=0.009, AMR=0.134, EAS=0.021 | 0.000 | 0.000 |
| HLA-B_nov_c4e5fc19 | novel | 18 | AFR=0.066, EUR=0.746, AMR=0.128, EAS=0.061 | 0.000 | 0.000 |
| HLA-B_nov_d5ae3969 | novel | 60 | AFR=0.052, EUR=0.556, AMR=0.086, EAS=0.306 | 0.000 | 0.000 |
| HLA-B_nov_3891497d | novel | 6 | AFR=0.015, EUR=0.000, AMR=0.980, EAS=0.005 | 0.000 | 0.001 |
| HLA-B_nov_ba21176a | novel | 21 | AFR=0.109, EUR=0.662, AMR=0.181, EAS=0.048 | 0.000 | 0.005 |
| HLA-B_nov_4cdba16a | novel | 12 | AFR=0.044, EUR=0.033, AMR=0.020, EAS=0.903 | 0.000 | 0.000 |
| HLA-B_nov_095f95ad | novel | 17 | AFR=0.028, EUR=0.023, AMR=0.021, EAS=0.929 | 0.000 | 0.000 |
| HLA-B_nov_d4f09b7d | novel | 5 | AFR=0.000, EUR=0.000, AMR=0.007, EAS=0.993 | 0.000 | 0.004 |
| HLA-B_nov_a58ba6da | novel | 16 | AFR=0.003, EUR=0.004, AMR=0.987, EAS=0.005 | 0.000 | 0.000 |
| HLA-B_nov_af8b859b | novel | 13 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.000 | 0.000 |
| HLA-B_nov_5484b26e | novel | 9 | AFR=0.004, EUR=0.045, AMR=0.012, EAS=0.940 | 0.000 | 0.000 |
| 15:20 | known | 9 | AFR=0.161, EUR=0.018, AMR=0.809, EAS=0.011 | 0.001 | 0.000 |
| 40:03 | known | 12 | AFR=0.017, EUR=0.054, AMR=0.413, EAS=0.516 | 0.001 | 0.001 |
| 40:10 | known | 7 | AFR=0.003, EUR=0.192, AMR=0.028, EAS=0.777 | 0.001 | 0.000 |
| 73:01 | known | 27 | AFR=0.005, EUR=0.517, AMR=0.348, EAS=0.130 | 0.001 | 0.000 |
| HLA-B_nov_5957fead | novel | 79 | AFR=0.460, EUR=0.203, AMR=0.296, EAS=0.040 | 0.001 | 0.000 |
| 18:04 | known | 7 | AFR=0.105, EUR=0.031, AMR=0.844, EAS=0.019 | 0.002 | 0.000 |
| 39:11 | known | 6 | AFR=0.039, EUR=0.008, AMR=0.936, EAS=0.017 | 0.002 | 0.001 |
| 57:04 | known | 18 | AFR=0.642, EUR=0.115, AMR=0.231, EAS=0.012 | 0.002 | 0.006 |
| 82:02 | known | 6 | AFR=0.957, EUR=0.005, AMR=0.024, EAS=0.014 | 0.002 | 0.004 |
| HLA-B_nov_2b2961f7 | novel | 6 | AFR=0.957, EUR=0.005, AMR=0.024, EAS=0.014 | 0.002 | 0.004 |
| 15:37 | known | 7 | AFR=0.842, EUR=0.062, AMR=0.066, EAS=0.030 | 0.003 | 0.012 |
| 15:39 | known | 4 | AFR=0.000, EUR=0.000, AMR=0.998, EAS=0.003 | 0.003 | 0.018 |
| HLA-B_nov_e056d52f | novel | 39 | AFR=0.028, EUR=0.457, AMR=0.284, EAS=0.232 | 0.003 | 0.000 |
| HLA-B_nov_d24423c2 | novel | 19 | AFR=0.406, EUR=0.050, AMR=0.515, EAS=0.029 | 0.003 | 0.000 |
| HLA-B_nov_f9c021a8 | novel | 8 | AFR=0.827, EUR=0.006, AMR=0.156, EAS=0.011 | 0.003 | 0.045 |
| HLA-B_nov_207ec506 | novel | 39 | AFR=0.350, EUR=0.155, AMR=0.475, EAS=0.019 | 0.004 | 0.000 |
| 40:27 | known | 7 | AFR=0.124, EUR=0.067, AMR=0.797, EAS=0.013 | 0.005 | 0.006 |
| 40:229 | known | 3 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.005 | 0.001 |
| 48:03 | known | 3 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.005 | 0.001 |
| 56:03 | known | 3 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.005 | 0.001 |
| 15:29 | known | 10 | AFR=0.012, EUR=0.453, AMR=0.012, EAS=0.523 | 0.006 | 0.000 |
| 41:03 | known | 6 | AFR=0.921, EUR=0.000, AMR=0.077, EAS=0.002 | 0.006 | 0.004 |
| 81:02 | known | 3 | AFR=0.000, EUR=0.003, AMR=0.000, EAS=0.997 | 0.006 | 0.001 |
| HLA-B_nov_509a423f | novel | 3 | AFR=0.000, EUR=0.003, AMR=0.000, EAS=0.997 | 0.006 | 0.001 |
| 35:11 | known | 6 | AFR=0.095, EUR=0.010, AMR=0.876, EAS=0.019 | 0.007 | 0.001 |
| HLA-B_nov_79772159 | novel | 6 | AFR=0.095, EUR=0.010, AMR=0.876, EAS=0.019 | 0.007 | 0.001 |
| 15:04 | known | 5 | AFR=0.065, EUR=0.032, AMR=0.881, EAS=0.022 | 0.009 | 0.005 |
| 40:16 | known | 8 | AFR=0.744, EUR=0.001, AMR=0.249, EAS=0.006 | 0.009 | 0.036 |
| 56:01 | known | 113 | AFR=0.201, EUR=0.398, AMR=0.197, EAS=0.204 | 0.009 | 0.029 |
| HLA-B_nov_35db003c | novel | 67 | AFR=0.229, EUR=0.454, AMR=0.275, EAS=0.042 | 0.009 | 0.003 |
| HLA-B_nov_59bcae82 | novel | 30 | AFR=0.165, EUR=0.552, AMR=0.137, EAS=0.146 | 0.011 | 0.001 |
| 15:05 | known | 15 | AFR=0.151, EUR=0.108, AMR=0.613, EAS=0.129 | 0.012 | 0.000 |
| 14:01 | known | 131 | AFR=0.314, EUR=0.282, AMR=0.344, EAS=0.060 | 0.013 | 0.000 |
| HLA-B_nov_f4ec8517 | novel | 87 | AFR=0.164, EUR=0.366, AMR=0.349, EAS=0.120 | 0.013 | 0.101 |
| 35:10 | known | 3 | AFR=0.003, EUR=0.000, AMR=0.993, EAS=0.003 | 0.014 | 0.070 |
| 18:145 | known | 8 | AFR=0.692, EUR=0.041, AMR=0.235, EAS=0.031 | 0.024 | 0.222 |
| 53:37 | known | 3 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.025 | 0.118 |
| 39:04 | known | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.014 |
| HLA-B_nov_fc3bd25b | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.225 |
| HLA-B_nov_bb5af146 | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.014 |
| HLA-B_nov_fef19996 | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.014 |
| HLA-B_nov_ef674346 | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.014 |
| HLA-B_nov_1e84ffb4 | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.026 | 0.014 |
| 55:04 | known | 2 | AFR=0.000, EUR=0.005, AMR=0.000, EAS=0.995 | 0.027 | 0.014 |
| HLA-B_nov_fce99c63 | novel | 2 | AFR=0.000, EUR=0.005, AMR=0.000, EAS=0.995 | 0.027 | 0.014 |
| 39:15 | known | 2 | AFR=0.000, EUR=0.012, AMR=0.049, EAS=0.938 | 0.029 | 0.014 |
| 15:08 | known | 12 | AFR=0.244, EUR=0.022, AMR=0.272, EAS=0.463 | 0.030 | 0.000 |
| 15:220 | known | 27 | AFR=0.204, EUR=0.526, AMR=0.242, EAS=0.028 | 0.035 | 0.000 |
| HLA-B_nov_b732cb1c | novel | 33 | AFR=0.194, EUR=0.495, AMR=0.264, EAS=0.047 | 0.035 | 0.000 |
| HLA-B_nov_51a1ac75 | novel | 129 | AFR=0.359, EUR=0.269, AMR=0.291, EAS=0.080 | 0.037 | 0.069 |
| 15:54 | known | 3 | AFR=0.970, EUR=0.000, AMR=0.027, EAS=0.003 | 0.041 | 0.118 |
| HLA-B_nov_0903559e | novel | 6 | AFR=0.192, EUR=0.767, AMR=0.037, EAS=0.004 | 0.044 | 0.000 |
| 35:30 | known | 4 | AFR=0.193, EUR=0.003, AMR=0.058, EAS=0.748 | 0.046 | 0.009 |
| HLA-B_nov_d0894e4e | novel | 25 | AFR=0.119, EUR=0.504, AMR=0.291, EAS=0.086 | 0.048 | 0.250 |
| 14:05 | known | 3 | AFR=0.956, EUR=0.007, AMR=0.020, EAS=0.017 | 0.050 | 0.118 |
| 35:06 | known | 3 | AFR=0.061, EUR=0.000, AMR=0.939, EAS=0.000 | 0.051 | 0.070 |
| 40:12 | known | 6 | AFR=0.739, EUR=0.138, AMR=0.066, EAS=0.057 | 0.054 | 0.013 |
| HLA-B_nov_8ecf43de | novel | 4 | AFR=0.310, EUR=0.000, AMR=0.020, EAS=0.670 | 0.055 | 0.181 |
| 15:07 | known | 11 | AFR=0.108, EUR=0.090, AMR=0.429, EAS=0.373 | 0.060 | 0.041 |
| 35:22 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.060 | 0.236 |
| 35:26 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.060 | 0.236 |
| 39:14 | known | 2 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.060 | 0.236 |
| HLA-B_nov_2329b766 | novel | 2 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.060 | 0.236 |
| HLA-B_nov_224d4cee | novel | 2 | AFR=0.000, EUR=0.000, AMR=0.995, EAS=0.005 | 0.060 | 0.236 |
| HLA-B_nov_7716f50b | novel | 100 | AFR=0.350, EUR=0.350, AMR=0.227, EAS=0.072 | 0.061 | 0.020 |
| HLA-B_nov_296d3d00 | novel | 2 | AFR=0.005, EUR=0.000, AMR=0.990, EAS=0.005 | 0.061 | 0.236 |
| 47:01 | known | 28 | AFR=0.215, EUR=0.502, AMR=0.240, EAS=0.043 | 0.063 | 0.005 |
| 40:11 | known | 3 | AFR=0.109, EUR=0.000, AMR=0.888, EAS=0.004 | 0.071 | 0.070 |
| HLA-B_nov_b99b4797 | novel | 57 | AFR=0.262, EUR=0.359, AMR=0.347, EAS=0.032 | 0.072 | 0.001 |
| HLA-B_nov_2d657694 | novel | 21 | AFR=0.104, EUR=0.498, AMR=0.226, EAS=0.172 | 0.080 | 0.069 |
| 58:11 | known | 4 | AFR=0.766, EUR=0.005, AMR=0.206, EAS=0.023 | 0.084 | 0.312 |
| 15:38 | known | 3 | AFR=0.021, EUR=0.072, AMR=0.845, EAS=0.062 | 0.093 | 0.044 |
| HLA-B_nov_0d837470 | novel | 4 | AFR=0.740, EUR=0.005, AMR=0.252, EAS=0.003 | 0.095 | 0.312 |
| HLA-B_nov_297f6a9c | novel | 4 | AFR=0.028, EUR=0.784, AMR=0.116, EAS=0.072 | 0.095 | 0.204 |
| 18:25 | known | 2 | AFR=0.995, EUR=0.000, AMR=0.005, EAS=0.000 | 0.096 | 0.320 |
| HLA-B_nov_2ec5cc9c | novel | 2 | AFR=0.025, EUR=0.000, AMR=0.975, EAS=0.000 | 0.097 | 0.236 |
| 15:24 | known | 7 | AFR=0.025, EUR=0.372, AMR=0.579, EAS=0.025 | 0.097 | 0.148 |
| 27:25 | known | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| 51:361N | known | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| 55:07 | known | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| 56:05 | known | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_8c43039b | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_7a05cd54 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_dddc4b92 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b8962c6b | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_a625005e | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_9d08ae6f | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_8ec2f62f | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_296e2bcf | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_5e16dd02 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_2e4f87ea | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_a9522e91 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_4c3ab01e | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_eaffd835 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_aec7bc6f | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_24d4ad14 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_942ad49d | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_d8a464bb | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_d6153940 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_429fc1b7 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_9cc4537a | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_515802b9 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_313a73ca | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_263d6ace | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_7227aa04 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_d62c1684 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_61b1817d | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_ab7a1970 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_cca9288a | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_6a4e0ed5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.121 |
| HLA-B_nov_8c6fe73c | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_ab839558 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_337f87e9 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_57717e15 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b0dbc94d | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_ab01cadc | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_0e2a36dc | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_8f012d91 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_088dcf0d | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_3a5db847 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_f6e168be | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_2e864900 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_67256714 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_843fa887 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_6afdcf32 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_6f6320ab | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_6e536944 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_163a5aba | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_1904755d | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_a9c9570f | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_79c3a6b5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_02ffb226 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_960090af | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_f06cf905 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_1c5e6303 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_3e88c405 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_40b6e9d2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_ead802fb | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_bc1468c9 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_95a6cace | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_e0f96dec | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_be498733 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_0c948668 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_38683399 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_57aa1303 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_97003b12 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_f65a63ce | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_fcebb155 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_c17448dc | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b57de30e | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_694dd122 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_c7440247 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_e5aa62e4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b8f68360 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b85da312 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_a6e4e751 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_e4a3c143 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_4145800f | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_6df32fb8 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_7d6911e6 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b5602c21 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_24a8297a | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_63c119eb | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_dc9657cb | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_843c5aa2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_b196ca06 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_fd209ef1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.000, EAS=1.000 | 0.101 | 0.208 |
| HLA-B_nov_766af2c2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.010, EAS=0.990 | 0.104 | 0.208 |
| HLA-B_nov_afb99d39 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.010, EAS=0.990 | 0.104 | 0.208 |
| HLA-B_nov_c30ccde0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.010, EAS=0.990 | 0.104 | 0.208 |
| 35:421 | known | 2 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.105 | 0.320 |
| HLA-B_nov_39c31ef1 | novel | 2 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.105 | 0.320 |
| 15:06 | known | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| HLA-B_nov_9e57f044 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| HLA-B_nov_c071b71a | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| HLA-B_nov_f6540be4 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| HLA-B_nov_fd28ebd8 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| HLA-B_nov_dd34539a | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.020, EAS=0.980 | 0.110 | 0.208 |
| HLA-B_nov_4e8ff699 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.000, EAS=0.990 | 0.110 | 0.208 |
| 56:02 | known | 1 | AFR=0.010, EUR=0.010, AMR=0.020, EAS=0.960 | 0.111 | 0.208 |
| HLA-B_nov_6fe77099 | novel | 1 | AFR=0.010, EUR=0.010, AMR=0.020, EAS=0.960 | 0.111 | 0.208 |
| HLA-B_nov_a98e7bc5 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.010, EAS=0.980 | 0.111 | 0.208 |
| HLA-B_nov_a2393038 | novel | 29 | AFR=0.114, EUR=0.418, AMR=0.266, EAS=0.203 | 0.122 | 0.007 |
| 51:07 | known | 14 | AFR=0.139, EUR=0.406, AMR=0.448, EAS=0.007 | 0.145 | 0.247 |
| 51:09 | known | 6 | AFR=0.623, EUR=0.301, AMR=0.043, EAS=0.033 | 0.160 | 0.212 |
| HLA-B_nov_3c2934f0 | novel | 41 | AFR=0.290, EUR=0.271, AMR=0.383, EAS=0.056 | 0.167 | 0.124 |
| HLA-B_nov_d8cdc4a8 | novel | 7 | AFR=0.165, EUR=0.497, AMR=0.003, EAS=0.335 | 0.169 | 0.368 |
| 15:40 | known | 2 | AFR=0.052, EUR=0.005, AMR=0.932, EAS=0.010 | 0.186 | 0.236 |
| HLA-B_nov_93bc130b | novel | 5 | AFR=0.000, EUR=0.436, AMR=0.545, EAS=0.019 | 0.199 | 0.015 |
| HLA-B_nov_3aa8d048 | novel | 2 | AFR=0.932, EUR=0.005, AMR=0.036, EAS=0.026 | 0.203 | 0.320 |
| HLA-B_nov_18e985ae | novel | 2 | AFR=0.000, EUR=0.123, AMR=0.063, EAS=0.814 | 0.205 | 0.014 |
| 15:108 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 15:63 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 35:24 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 35:332 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 39:75 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 44:521 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| 51:127 | known | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_03e38808 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_32a5c027 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_3863c439 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_23ed1117 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_53b6524c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_d40006ff | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_17b10f4b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_fd89aa93 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_4e217d42 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_feb4cad1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_85bdf62f | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_dc7878dd | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_40a4247f | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_e91987c0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_c28cabe3 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_497bf2e5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_ea3116e4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_8f4eaa30 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_45b8df76 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_ea54482e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_6d208db1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_9baf3e13 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_61fbc128 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_79d30be6 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_04733de3 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_c94dd52e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_c75a2c34 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_dace3103 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_0fee0b6b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_45aa1780 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_bd37224e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_4bf64393 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_216797a7 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1a1089bc | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_a22c4f5b | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_4eb4ae8c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1f65d474 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1879609c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_f95cbd0e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_b3d0e4df | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_9abb80b8 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1aa951a6 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1aefebe5 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_4708e0d4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_f96ccee4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1da6c7ff | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_921de884 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_19a19987 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_d463760d | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_733cd3d4 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_6d42ea59 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_043f30a2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_c0b4cb4c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_b9a1c31c | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_366b17a1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_bc903ccf | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_082d2f6e | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_5ae13cb0 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_4921d5c2 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_c646eb03 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_374a687a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_0231069a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1728bc61 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_92d71921 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_b03a9586 | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_1652691a | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_411bb40d | novel | 1 | AFR=0.000, EUR=0.000, AMR=1.000, EAS=0.000 | 0.207 | 0.639 |
| HLA-B_nov_52b829e0 | novel | 12 | AFR=0.355, EUR=0.498, AMR=0.123, EAS=0.024 | 0.207 | 0.081 |
| HLA-B_nov_4fdfbb62 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.990, EAS=0.010 | 0.211 | 0.639 |
| HLA-B_nov_e3128d14 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.990, EAS=0.010 | 0.211 | 0.639 |
| 39:13 | known | 1 | AFR=0.000, EUR=0.000, AMR=0.990, EAS=0.010 | 0.213 | 0.639 |
| HLA-B_nov_34a394d7 | novel | 15 | AFR=0.166, EUR=0.484, AMR=0.321, EAS=0.029 | 0.216 | 0.279 |
| 15:47 | known | 2 | AFR=0.913, EUR=0.000, AMR=0.067, EAS=0.021 | 0.222 | 0.320 |
| HLA-B_nov_e7324a00 | novel | 1 | AFR=0.010, EUR=0.000, AMR=0.990, EAS=0.000 | 0.224 | 0.639 |
| 07:12 | known | 2 | AFR=0.907, EUR=0.010, AMR=0.052, EAS=0.031 | 0.225 | 0.320 |
| HLA-B_nov_57b7bbb4 | novel | 12 | AFR=0.109, EUR=0.400, AMR=0.433, EAS=0.057 | 0.225 | 0.337 |
| HLA-B_nov_67071646 | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.990, EAS=0.000 | 0.230 | 0.639 |
| HLA-B_nov_6f552993 | novel | 7 | AFR=0.510, EUR=0.054, AMR=0.356, EAS=0.081 | 0.233 | 0.806 |
| HLA-B_nov_5a90260b | novel | 2 | AFR=0.000, EUR=0.914, AMR=0.069, EAS=0.017 | 0.240 | 0.031 |
| 18:03 | known | 5 | AFR=0.576, EUR=0.384, AMR=0.022, EAS=0.018 | 0.248 | 0.386 |
| HLA-B_nov_6c52313d | novel | 6 | AFR=0.022, EUR=0.571, AMR=0.192, EAS=0.214 | 0.251 | 0.484 |
| HLA-B_nov_6969be6c | novel | 94 | AFR=0.329, EUR=0.348, AMR=0.218, EAS=0.105 | 0.256 | 0.184 |
| HLA-B_nov_fc728f83 | novel | 2 | AFR=0.102, EUR=0.010, AMR=0.847, EAS=0.041 | 0.262 | 0.236 |
| HLA-B_nov_dbe98bfc | novel | 2 | AFR=0.006, EUR=0.882, AMR=0.080, EAS=0.032 | 0.269 | 0.313 |
| 50:02 | known | 16 | AFR=0.188, EUR=0.330, AMR=0.439, EAS=0.043 | 0.274 | 0.002 |
| HLA-B_nov_a4ba84de | novel | 9 | AFR=0.312, EUR=0.045, AMR=0.417, EAS=0.226 | 0.275 | 0.271 |
| 35:28 | known | 2 | AFR=0.071, EUR=0.061, AMR=0.828, EAS=0.040 | 0.278 | 0.236 |
| HLA-B_nov_77a21fe5 | novel | 2 | AFR=0.143, EUR=0.015, AMR=0.811, EAS=0.031 | 0.290 | 0.236 |
| HLA-B_nov_f72ccc05 | novel | 2 | AFR=0.000, EUR=0.847, AMR=0.106, EAS=0.047 | 0.291 | 0.436 |
| 07:09 | known | 4 | AFR=0.340, EUR=0.036, AMR=0.597, EAS=0.027 | 0.292 | 0.503 |
| HLA-B_nov_1cad5d70 | novel | 56 | AFR=0.184, EUR=0.312, AMR=0.277, EAS=0.227 | 0.299 | 0.037 |
| 39:20 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| 44:28 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| 45:07 | known | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| 51:237 | known | 1 | AFR=0.000, EUR=0.000, AMR=0.980, EAS=0.020 | 0.309 | 0.639 |
| HLA-B_nov_15e124e2 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_3acc735a | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.980, EAS=0.020 | 0.309 | 0.639 |
| HLA-B_nov_e048626d | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_f760e3e5 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_5bde10d6 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_93267114 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_6a0d3300 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_91625672 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_9f7d2aaf | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_43bbcf50 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_b4bcc5dc | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_d1683cd8 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_628c88c7 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_df3a30ba | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_84a80f13 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_8dc202a4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_42312a88 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_1a014e8d | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_6f53b6c4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_ffb8b272 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_cb50fa9c | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_cf249fe3 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_60472dca | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_94facc98 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_44eb7f8e | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_42309e46 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_bb3d5087 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_3ecf65c4 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_72cdf352 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_920ce8d1 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_2b5013b1 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_603eb428 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_2e8a159a | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_950f781a | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_09b874e9 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_5e41ec2f | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_c8698013 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_1202a22c | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_114217e3 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_92c9d9ec | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_565856be | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_db5506d5 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_8343f020 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_a14c370b | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_0fb8eb70 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_b1d82175 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_c2686f98 | novel | 1 | AFR=1.000, EUR=0.000, AMR=0.000, EAS=0.000 | 0.309 | 0.711 |
| HLA-B_nov_813c65bc | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.980, EAS=0.010 | 0.312 | 0.639 |
| HLA-B_nov_6b34840c | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.980, EAS=0.010 | 0.312 | 0.639 |
| HLA-B_nov_aaded50a | novel | 1 | AFR=0.000, EUR=0.010, AMR=0.980, EAS=0.010 | 0.313 | 0.639 |
| HLA-B_nov_83b7b3bf | novel | 1 | AFR=0.020, EUR=0.000, AMR=0.980, EAS=0.000 | 0.317 | 0.639 |
| HLA-B_nov_f8b10fb2 | novel | 3 | AFR=0.184, EUR=0.076, AMR=0.682, EAS=0.057 | 0.317 | 0.070 |
| HLA-B_nov_03a90688 | novel | 1 | AFR=0.000, EUR=0.020, AMR=0.980, EAS=0.000 | 0.320 | 0.639 |
| HLA-B_nov_1e22404d | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.000, EAS=0.010 | 0.320 | 0.711 |
| HLA-B_nov_abcae738 | novel | 1 | AFR=0.000, EUR=0.020, AMR=0.980, EAS=0.000 | 0.320 | 0.639 |
| HLA-B_nov_7dfe68fe | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_1dd91039 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_edcf0038 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_15c1d0aa | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_b277757c | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_5aff1cd6 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_33b86f8b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_c51dd187 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_90c0e76f | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_5c3c76bb | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_9f907b71 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_41152288 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_3e375010 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_995a3ec9 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_eee296c7 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_f2462ce0 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_dba89d4e | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_0fedbd23 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_534167c6 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_bcdf3834 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_4e13f26f | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_93fcb251 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_3dfbd43b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_00d4b295 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_617fa75b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_1f821d4b | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_c9864187 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_8b4d2aff | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_72428e2d | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_d31fdb30 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.359 | 0.711 |
| HLA-B_nov_1928d7e3 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_b7dcf573 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_0dc0b6f4 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_0c2c9d42 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_541b65b4 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_1a2f8bcc | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.363 | 0.711 |
| HLA-B_nov_0d9622c0 | novel | 1 | AFR=0.990, EUR=0.000, AMR=0.010, EAS=0.000 | 0.366 | 0.711 |
| HLA-B_nov_3414ed58 | novel | 1 | AFR=0.990, EUR=0.010, AMR=0.000, EAS=0.000 | 0.367 | 0.711 |
| 15:09 | known | 7 | AFR=0.162, EUR=0.307, AMR=0.512, EAS=0.019 | 0.369 | 0.016 |
| HLA-B_nov_b4af01aa | novel | 1 | AFR=0.010, EUR=0.010, AMR=0.970, EAS=0.010 | 0.370 | 0.639 |
| HLA-B_nov_158ad9c3 | novel | 1 | AFR=0.000, EUR=0.030, AMR=0.970, EAS=0.000 | 0.370 | 0.639 |
| HLA-B_nov_d27b7922 | novel | 1 | AFR=0.020, EUR=0.000, AMR=0.970, EAS=0.010 | 0.370 | 0.639 |
| HLA-B_nov_48134f61 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.372 | 0.711 |
| HLA-B_nov_b49f3bb5 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.372 | 0.711 |
| HLA-B_nov_69a8b14c | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.372 | 0.711 |
| HLA-B_nov_1d6684b4 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.372 | 0.711 |
| HLA-B_nov_7f7bcbc8 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.010, EAS=0.010 | 0.372 | 0.711 |
| HLA-B_nov_53f13a25 | novel | 7 | AFR=0.525, EUR=0.167, AMR=0.275, EAS=0.033 | 0.382 | 0.459 |
| HLA-B_nov_d0b27b55 | novel | 6 | AFR=0.502, EUR=0.122, AMR=0.363, EAS=0.013 | 0.386 | 0.625 |
| HLA-B_nov_c75cde5b | novel | 4 | AFR=0.030, EUR=0.508, AMR=0.449, EAS=0.013 | 0.390 | 0.498 |
| 08:12 | known | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_a2421080 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_ae28592c | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_95c80fc7 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_45078566 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_a7aafd80 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_af6a8541 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_654c769f | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_5ee74302 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_30e22019 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_4641b6bc | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_82fa1e91 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_14c347fc | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_01d74970 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_a875de23 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_39af4334 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_04e3516e | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.397 | 0.711 |
| HLA-B_nov_89c83114 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.399 | 0.711 |
| HLA-B_nov_d35b12ee | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.399 | 0.711 |
| HLA-B_nov_864a10f4 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.403 | 0.711 |
| HLA-B_nov_1c3ee464 | novel | 1 | AFR=0.980, EUR=0.010, AMR=0.010, EAS=0.000 | 0.403 | 0.711 |
| HLA-B_nov_a93b0a2b | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.403 | 0.711 |
| HLA-B_nov_cda60249 | novel | 1 | AFR=0.980, EUR=0.000, AMR=0.020, EAS=0.000 | 0.403 | 0.711 |
| HLA-B_nov_78748980 | novel | 1 | AFR=0.020, EUR=0.010, AMR=0.960, EAS=0.010 | 0.404 | 0.639 |
| HLA-B_nov_a7c410d9 | novel | 1 | AFR=0.979, EUR=0.010, AMR=0.010, EAS=0.000 | 0.404 | 0.711 |
| HLA-B_nov_f89fa8c2 | novel | 1 | AFR=0.010, EUR=0.020, AMR=0.960, EAS=0.010 | 0.404 | 0.639 |
| HLA-B_nov_9a4b1870 | novel | 2 | AFR=0.096, EUR=0.000, AMR=0.314, EAS=0.590 | 0.410 | 0.509 |
| 47:03 | known | 3 | AFR=0.650, EUR=0.295, AMR=0.010, EAS=0.045 | 0.422 | 0.063 |
| 51:05 | known | 5 | AFR=0.136, EUR=0.551, AMR=0.074, EAS=0.240 | 0.433 | 0.355 |
| HLA-B_nov_e3fc05f4 | novel | 2 | AFR=0.005, EUR=0.495, AMR=0.000, EAS=0.500 | 0.447 | 0.539 |
| HLA-B_nov_28a66739 | novel | 9 | AFR=0.472, EUR=0.225, AMR=0.288, EAS=0.015 | 0.459 | 0.528 |
| HLA-B_nov_b5feb9c0 | novel | 2 | AFR=0.311, EUR=0.005, AMR=0.684, EAS=0.000 | 0.459 | 0.236 |
| HLA-B_nov_3a379e69 | novel | 2 | AFR=0.011, EUR=0.358, AMR=0.059, EAS=0.572 | 0.467 | 0.539 |
| 07:160 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| 27:01 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| 35:41 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| 41:05 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| 41:09 | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| 44:267N | known | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_6cbd58c1 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_d8ed128a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_8ae66ad7 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_dc21985f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_886a10c6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_23ce36a6 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_5dd52ff0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_8b074c29 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_0e304d54 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_ef0da28e | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_f02db3fe | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_6312be7c | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_088fccbd | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_2755f44b | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_432883d0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_5e419efa | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_d2bfdcde | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_7b58694e | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_ee45d146 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_c368473f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_bef98926 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_5e724eec | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_b473fb45 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_51f79b75 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_4ed134e4 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_c914b066 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_65dd51de | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_77dab031 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_c1c31457 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_09c249df | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_79869e0a | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_92974912 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_8c14ef54 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_a05cc281 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_7373a5dd | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_2df45ea3 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_cc3f5e23 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_20e7ad3f | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_d4c5cfff | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_ab3dd779 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_a8d3b5ae | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_639c4dd0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_67648ed0 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_a0d04fdc | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_248c0970 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.000 |
| HLA-B_nov_0e39d3d2 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_fb611334 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_0628de16 | novel | 1 | AFR=0.000, EUR=1.000, AMR=0.000, EAS=0.000 | 0.484 | 0.705 |
| HLA-B_nov_d4bfb1c6 | novel | 1 | AFR=0.010, EUR=0.020, AMR=0.960, EAS=0.010 | 0.487 | 0.639 |
| HLA-B_nov_59d5eed5 | novel | 1 | AFR=0.000, EUR=0.040, AMR=0.960, EAS=0.000 | 0.488 | 0.639 |
| HLA-B_nov_c272771c | novel | 1 | AFR=0.000, EUR=0.040, AMR=0.960, EAS=0.000 | 0.488 | 0.639 |
| HLA-B_nov_95c48066 | novel | 1 | AFR=0.041, EUR=0.000, AMR=0.959, EAS=0.000 | 0.488 | 0.639 |
| HLA-B_nov_e3452f6d | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.010, EAS=0.020 | 0.491 | 0.711 |
| HLA-B_nov_4fb34829 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.010, EAS=0.020 | 0.491 | 0.711 |
| HLA-B_nov_9fa331de | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.010, EAS=0.020 | 0.491 | 0.711 |
| HLA-B_nov_90aa9c07 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.020, EAS=0.010 | 0.498 | 0.711 |
| HLA-B_nov_9c8d67db | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.020, EAS=0.010 | 0.498 | 0.711 |
| HLA-B_nov_f8139703 | novel | 1 | AFR=0.020, EUR=0.000, AMR=0.950, EAS=0.030 | 0.499 | 0.639 |
| HLA-B_nov_b6e49b3c | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.020, EAS=0.010 | 0.499 | 0.711 |
| 15:52 | known | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_3123b867 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_656b54ac | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_9a725575 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_fb966272 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_1ef08883 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_33724e64 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_e9ce0e06 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_97890c08 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.506 | 0.711 |
| HLA-B_nov_f1d08636 | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.508 | 0.711 |
| HLA-B_nov_8604a1ed | novel | 1 | AFR=0.970, EUR=0.000, AMR=0.030, EAS=0.000 | 0.508 | 0.711 |
| HLA-B_nov_e74d1f15 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031, EAS=0.000 | 0.509 | 0.711 |
| HLA-B_nov_5e2c4f54 | novel | 1 | AFR=0.969, EUR=0.000, AMR=0.031, EAS=0.000 | 0.509 | 0.711 |
| HLA-B_nov_5de4168d | novel | 2 | AFR=0.099, EUR=0.005, AMR=0.396, EAS=0.500 | 0.517 | 0.509 |
| 56:94 | known | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_cdb36446 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_617ca186 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_03a5909a | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_adac2bc6 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_6edc2a7d | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_beb4db1f | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_7fe285e8 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_7da43fcf | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_09d9054d | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_41e2694d | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_ed61f8b6 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.000, EAS=0.010 | 0.526 | 0.705 |
| HLA-B_nov_41b5eac8 | novel | 1 | AFR=0.010, EUR=0.020, AMR=0.949, EAS=0.020 | 0.527 | 0.639 |
| HLA-B_nov_7c4a4f49 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_2f24cf74 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_b947dcb4 | novel | 1 | AFR=0.010, EUR=0.990, AMR=0.000, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_5d00e6a1 | novel | 1 | AFR=0.010, EUR=0.990, AMR=0.000, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_c6d6cdb1 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_8ddb92d2 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_40247a71 | novel | 1 | AFR=0.000, EUR=0.990, AMR=0.010, EAS=0.000 | 0.537 | 0.705 |
| HLA-B_nov_457dac52 | novel | 1 | AFR=0.000, EUR=0.051, AMR=0.949, EAS=0.000 | 0.538 | 0.639 |
| HLA-B_nov_94f7382f | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.020, EAS=0.020 | 0.541 | 0.711 |
| 35:542N | known | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.547 | 0.711 |
| HLA-B_nov_3e66ec37 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.547 | 0.711 |
| HLA-B_nov_a5d0748a | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.547 | 0.711 |
| HLA-B_nov_562d19f4 | novel | 1 | AFR=0.000, EUR=0.020, AMR=0.939, EAS=0.040 | 0.547 | 0.639 |
| HLA-B_nov_85407195 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.547 | 0.711 |
| HLA-B_nov_fd7ba8ef | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.547 | 0.711 |
| HLA-B_nov_4f621ce7 | novel | 1 | AFR=0.958, EUR=0.000, AMR=0.021, EAS=0.021 | 0.548 | 0.711 |
| HLA-B_nov_1a7cc0c5 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.548 | 0.711 |
| HLA-B_nov_415f5cea | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.031, EAS=0.010 | 0.548 | 0.711 |
| HLA-B_nov_af64f753 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.030, EAS=0.010 | 0.548 | 0.711 |
| 08:16 | known | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_cecb3b95 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_5042bbd5 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_f3b4f0b2 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_3ecd044d | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_900e5341 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.556 | 0.711 |
| HLA-B_nov_6de4a489 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.561 | 0.711 |
| HLA-B_nov_89319c90 | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.561 | 0.711 |
| HLA-B_nov_be69758b | novel | 1 | AFR=0.960, EUR=0.000, AMR=0.040, EAS=0.000 | 0.561 | 0.711 |
| HLA-B_nov_709461de | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.041, EAS=0.000 | 0.561 | 0.711 |
| HLA-B_nov_9f3e174f | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.000, EAS=0.020 | 0.562 | 0.705 |
| HLA-B_nov_606353a9 | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.000, EAS=0.020 | 0.562 | 0.705 |
| HLA-B_nov_785793cf | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.000, EAS=0.020 | 0.562 | 0.705 |
| HLA-B_nov_f7a1452f | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.000, EAS=0.020 | 0.562 | 0.705 |
| HLA-B_nov_5678c505 | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.000, EAS=0.020 | 0.562 | 0.705 |
| 15:31 | known | 2 | AFR=0.485, EUR=0.000, AMR=0.515, EAS=0.000 | 0.563 | 0.826 |
| HLA-B_nov_e6e2087c | novel | 1 | AFR=0.959, EUR=0.000, AMR=0.041, EAS=0.000 | 0.563 | 0.711 |
| HLA-B_nov_584448c0 | novel | 1 | AFR=0.030, EUR=0.010, AMR=0.939, EAS=0.020 | 0.563 | 0.639 |
| HLA-B_nov_e05f2aaf | novel | 1 | AFR=0.000, EUR=0.979, AMR=0.000, EAS=0.021 | 0.563 | 0.705 |
| HLA-B_nov_8e8d3c8f | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.010, EAS=0.010 | 0.568 | 0.705 |
| HLA-B_nov_b42ddaac | novel | 1 | AFR=0.020, EUR=0.020, AMR=0.939, EAS=0.020 | 0.568 | 0.639 |
| HLA-B_nov_7d255483 | novel | 1 | AFR=0.957, EUR=0.011, AMR=0.021, EAS=0.011 | 0.568 | 0.711 |
| HLA-B_nov_4c450148 | novel | 1 | AFR=0.051, EUR=0.000, AMR=0.939, EAS=0.010 | 0.568 | 0.639 |
| 44:302 | known | 1 | AFR=0.010, EUR=0.980, AMR=0.010, EAS=0.000 | 0.569 | 0.705 |
| HLA-B_nov_b44619c9 | novel | 1 | AFR=0.000, EUR=0.980, AMR=0.020, EAS=0.000 | 0.569 | 0.705 |
| HLA-B_nov_183b774a | novel | 1 | AFR=0.010, EUR=0.980, AMR=0.010, EAS=0.000 | 0.569 | 0.705 |
| 53:15 | known | 1 | AFR=0.950, EUR=0.000, AMR=0.030, EAS=0.020 | 0.570 | 0.711 |
| HLA-B_nov_a32b50be | novel | 1 | AFR=0.954, EUR=0.023, AMR=0.023, EAS=0.000 | 0.570 | 0.711 |
| HLA-B_nov_8c0e6332 | novel | 1 | AFR=0.950, EUR=0.000, AMR=0.030, EAS=0.020 | 0.570 | 0.711 |
| HLA-B_nov_8345fbf0 | novel | 1 | AFR=0.949, EUR=0.000, AMR=0.020, EAS=0.030 | 0.570 | 0.711 |
| HLA-B_nov_65b3cfad | novel | 1 | AFR=0.000, EUR=0.048, AMR=0.120, EAS=0.831 | 0.570 | 0.208 |
| HLA-B_nov_92dad439 | novel | 1 | AFR=0.950, EUR=0.000, AMR=0.030, EAS=0.020 | 0.570 | 0.711 |
| HLA-B_nov_b4e9b4be | novel | 1 | AFR=0.950, EUR=0.000, AMR=0.040, EAS=0.010 | 0.572 | 0.711 |
| HLA-B_nov_20b30691 | novel | 1 | AFR=0.948, EUR=0.000, AMR=0.031, EAS=0.021 | 0.573 | 0.711 |
| HLA-B_nov_13b01f0f | novel | 1 | AFR=0.950, EUR=0.000, AMR=0.050, EAS=0.000 | 0.574 | 0.711 |
| HLA-B_nov_11794504 | novel | 1 | AFR=0.949, EUR=0.000, AMR=0.051, EAS=0.000 | 0.576 | 0.711 |
| HLA-B_nov_bf1a69a3 | novel | 1 | AFR=0.000, EUR=0.970, AMR=0.010, EAS=0.020 | 0.580 | 0.705 |
| HLA-B_nov_89424f29 | novel | 1 | AFR=0.948, EUR=0.000, AMR=0.052, EAS=0.000 | 0.580 | 0.711 |
| 07:10 | known | 1 | AFR=0.000, EUR=0.970, AMR=0.010, EAS=0.020 | 0.582 | 0.705 |
| HLA-B_nov_4d179b4c | novel | 1 | AFR=0.010, EUR=0.030, AMR=0.929, EAS=0.030 | 0.582 | 0.639 |
| HLA-B_nov_318ba7f2 | novel | 1 | AFR=0.000, EUR=0.970, AMR=0.010, EAS=0.020 | 0.582 | 0.705 |
| HLA-B_nov_dadd5742 | novel | 1 | AFR=0.000, EUR=0.971, AMR=0.029, EAS=0.000 | 0.586 | 0.000 |
| 39:31 | known | 1 | AFR=0.000, EUR=0.967, AMR=0.016, EAS=0.016 | 0.587 | 0.705 |
| HLA-B_nov_aedf1b6f | novel | 1 | AFR=0.944, EUR=0.022, AMR=0.011, EAS=0.022 | 0.587 | 0.711 |
| HLA-B_nov_a5149631 | novel | 1 | AFR=0.031, EUR=0.020, AMR=0.929, EAS=0.020 | 0.587 | 0.639 |
| HLA-B_nov_f07b2d2c | novel | 1 | AFR=0.000, EUR=0.968, AMR=0.032, EAS=0.000 | 0.587 | 0.705 |
| HLA-B_nov_7263fdef | novel | 1 | AFR=0.062, EUR=0.000, AMR=0.927, EAS=0.010 | 0.587 | 0.639 |
| HLA-B_nov_d5349780 | novel | 1 | AFR=0.939, EUR=0.000, AMR=0.040, EAS=0.020 | 0.590 | 0.711 |
| HLA-B_nov_6e966e59 | novel | 1 | AFR=0.000, EUR=0.960, AMR=0.000, EAS=0.040 | 0.590 | 0.705 |
| HLA-B_nov_c6839d69 | novel | 1 | AFR=0.939, EUR=0.000, AMR=0.041, EAS=0.020 | 0.591 | 0.711 |
| HLA-B_nov_442f4a56 | novel | 1 | AFR=0.937, EUR=0.000, AMR=0.032, EAS=0.032 | 0.592 | 0.711 |
| HLA-B_nov_180180d4 | novel | 1 | AFR=0.940, EUR=0.000, AMR=0.060, EAS=0.000 | 0.592 | 0.711 |
| HLA-B_nov_1caccbe8 | novel | 1 | AFR=0.939, EUR=0.010, AMR=0.031, EAS=0.020 | 0.592 | 0.711 |
| HLA-B_nov_d53f646f | novel | 1 | AFR=0.939, EUR=0.000, AMR=0.061, EAS=0.000 | 0.597 | 0.711 |
| HLA-B_nov_8df7809d | novel | 1 | AFR=0.939, EUR=0.010, AMR=0.051, EAS=0.000 | 0.597 | 0.711 |
| HLA-B_nov_09eb855b | novel | 1 | AFR=0.000, EUR=0.960, AMR=0.020, EAS=0.020 | 0.597 | 0.705 |
| HLA-B_nov_1f850e58 | novel | 2 | AFR=0.505, EUR=0.000, AMR=0.490, EAS=0.005 | 0.597 | 0.826 |
| HLA-B_nov_42b8e4ab | novel | 1 | AFR=0.017, EUR=0.050, AMR=0.117, EAS=0.817 | 0.599 | 0.208 |
| HLA-B_nov_3be52441 | novel | 1 | AFR=0.081, EUR=0.000, AMR=0.919, EAS=0.000 | 0.599 | 0.639 |
| HLA-B_nov_a915fad3 | novel | 1 | AFR=0.010, EUR=0.960, AMR=0.020, EAS=0.010 | 0.599 | 0.705 |
| HLA-B_nov_249d46ed | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.030, EAS=0.040 | 0.604 | 0.711 |
| HLA-B_nov_7e151b9e | novel | 1 | AFR=0.930, EUR=0.000, AMR=0.050, EAS=0.020 | 0.607 | 0.711 |
| HLA-B_nov_17512224 | novel | 1 | AFR=0.929, EUR=0.010, AMR=0.020, EAS=0.041 | 0.607 | 0.711 |
| HLA-B_nov_7d2945bd | novel | 1 | AFR=0.930, EUR=0.000, AMR=0.060, EAS=0.010 | 0.609 | 0.711 |
| 51:285 | known | 1 | AFR=0.020, EUR=0.020, AMR=0.909, EAS=0.051 | 0.610 | 0.639 |
| HLA-B_nov_c3e62ffe | novel | 1 | AFR=0.929, EUR=0.010, AMR=0.031, EAS=0.031 | 0.610 | 0.711 |
| HLA-B_nov_726dbc56 | novel | 1 | AFR=0.020, EUR=0.020, AMR=0.909, EAS=0.051 | 0.610 | 0.639 |
| HLA-B_nov_ce029c60 | novel | 1 | AFR=0.000, EUR=0.950, AMR=0.010, EAS=0.040 | 0.610 | 0.705 |
| HLA-B_nov_f12cfa5b | novel | 1 | AFR=0.930, EUR=0.010, AMR=0.050, EAS=0.010 | 0.610 | 0.711 |
| HLA-B_nov_f46b406f | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.071, EAS=0.000 | 0.611 | 0.711 |
| HLA-B_nov_ac5a3f46 | novel | 1 | AFR=0.929, EUR=0.000, AMR=0.061, EAS=0.010 | 0.611 | 0.711 |
| HLA-B_nov_1a6efd18 | novel | 1 | AFR=0.000, EUR=0.949, AMR=0.020, EAS=0.030 | 0.611 | 0.705 |
| 44:21 | known | 1 | AFR=0.000, EUR=0.952, AMR=0.048, EAS=0.000 | 0.613 | 0.000 |
| HLA-B_nov_28d798e5 | novel | 1 | AFR=0.040, EUR=0.020, AMR=0.910, EAS=0.030 | 0.613 | 0.639 |
| HLA-B_nov_e054ddf9 | novel | 1 | AFR=0.010, EUR=0.050, AMR=0.910, EAS=0.030 | 0.617 | 0.639 |
| HLA-B_nov_28799930 | novel | 1 | AFR=0.041, EUR=0.010, AMR=0.907, EAS=0.041 | 0.618 | 0.639 |
| HLA-B_nov_1f3faf7b | novel | 1 | AFR=0.041, EUR=0.010, AMR=0.907, EAS=0.041 | 0.618 | 0.639 |
| HLA-B_nov_c2365afe | novel | 1 | AFR=0.051, EUR=0.010, AMR=0.908, EAS=0.031 | 0.618 | 0.639 |
| HLA-B_nov_8fec16fc | novel | 1 | AFR=0.061, EUR=0.010, AMR=0.909, EAS=0.020 | 0.618 | 0.639 |
| HLA-B_nov_c0671265 | novel | 1 | AFR=0.000, EUR=0.949, AMR=0.030, EAS=0.020 | 0.618 | 0.705 |
| HLA-B_nov_2a7661a6 | novel | 1 | AFR=0.010, EUR=0.950, AMR=0.030, EAS=0.010 | 0.619 | 0.705 |
| HLA-B_nov_e524bf6b | novel | 1 | AFR=0.010, EUR=0.950, AMR=0.030, EAS=0.010 | 0.619 | 0.705 |
| HLA-B_nov_7d8c1f06 | novel | 2 | AFR=0.665, EUR=0.143, AMR=0.081, EAS=0.110 | 0.619 | 0.320 |
| HLA-B_nov_e4455130 | novel | 19 | AFR=0.283, EUR=0.190, AMR=0.275, EAS=0.252 | 0.620 | 0.000 |
| HLA-B_nov_3cdacf93 | novel | 1 | AFR=0.012, EUR=0.950, AMR=0.037, EAS=0.000 | 0.620 | 0.705 |
| HLA-B_nov_c177c874 | novel | 1 | AFR=0.920, EUR=0.000, AMR=0.030, EAS=0.050 | 0.621 | 0.711 |
| HLA-B_nov_f468ba32 | novel | 1 | AFR=0.918, EUR=0.000, AMR=0.020, EAS=0.061 | 0.621 | 0.711 |
| HLA-B_nov_92b7a05e | novel | 1 | AFR=0.031, EUR=0.041, AMR=0.908, EAS=0.020 | 0.621 | 0.639 |
| HLA-B_nov_3893b442 | novel | 1 | AFR=0.062, EUR=0.010, AMR=0.906, EAS=0.021 | 0.621 | 0.639 |
| HLA-B_nov_c60316d0 | novel | 1 | AFR=0.062, EUR=0.010, AMR=0.906, EAS=0.021 | 0.621 | 0.639 |
| HLA-B_nov_3d91afb3 | novel | 1 | AFR=0.920, EUR=0.000, AMR=0.040, EAS=0.040 | 0.621 | 0.711 |
| HLA-B_nov_73714420 | novel | 1 | AFR=0.000, EUR=0.946, AMR=0.054, EAS=0.000 | 0.623 | 0.705 |
| HLA-B_nov_856c84fc | novel | 1 | AFR=0.043, EUR=0.021, AMR=0.904, EAS=0.032 | 0.626 | 0.639 |
| HLA-B_nov_ad8338b7 | novel | 1 | AFR=0.924, EUR=0.043, AMR=0.022, EAS=0.011 | 0.627 | 0.711 |
| HLA-B_nov_eb1b4e4b | novel | 1 | AFR=0.918, EUR=0.000, AMR=0.041, EAS=0.041 | 0.629 | 0.711 |
| HLA-B_nov_f3cd8945 | novel | 1 | AFR=0.920, EUR=0.000, AMR=0.080, EAS=0.000 | 0.629 | 0.711 |
| HLA-B_nov_fa36f7c9 | novel | 1 | AFR=0.920, EUR=0.010, AMR=0.050, EAS=0.020 | 0.629 | 0.711 |
| HLA-B_nov_36457ef8 | novel | 1 | AFR=0.920, EUR=0.000, AMR=0.070, EAS=0.010 | 0.629 | 0.711 |
| HLA-B_nov_17aef3f8 | novel | 1 | AFR=0.000, EUR=0.938, AMR=0.000, EAS=0.062 | 0.629 | 0.000 |
| HLA-B_nov_cbda84f5 | novel | 1 | AFR=0.000, EUR=0.942, AMR=0.035, EAS=0.023 | 0.630 | 0.705 |
| HLA-B_nov_75cec853 | novel | 1 | AFR=0.000, EUR=0.940, AMR=0.020, EAS=0.040 | 0.630 | 0.705 |
| HLA-B_nov_ea5e695b | novel | 1 | AFR=0.919, EUR=0.000, AMR=0.071, EAS=0.010 | 0.630 | 0.711 |
| HLA-B_nov_cb340e43 | novel | 1 | AFR=0.000, EUR=0.940, AMR=0.030, EAS=0.030 | 0.632 | 0.705 |
| HLA-B_nov_6ca27245 | novel | 1 | AFR=0.010, EUR=0.940, AMR=0.020, EAS=0.030 | 0.633 | 0.705 |
| HLA-B_nov_7d4604b5 | novel | 1 | AFR=0.000, EUR=0.939, AMR=0.061, EAS=0.000 | 0.637 | 0.705 |
| HLA-B_nov_f1efcbf8 | novel | 1 | AFR=0.000, EUR=0.939, AMR=0.040, EAS=0.020 | 0.637 | 0.705 |
| 07:14 | known | 1 | AFR=0.010, EUR=0.000, AMR=0.880, EAS=0.110 | 0.638 | 0.639 |
| HLA-B_nov_4f805aaa | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.010, EAS=0.061 | 0.639 | 0.705 |
| HLA-B_nov_f4735665 | novel | 1 | AFR=0.000, EUR=0.930, AMR=0.020, EAS=0.050 | 0.640 | 0.705 |
| HLA-B_nov_cbae2f89 | novel | 1 | AFR=0.000, EUR=0.930, AMR=0.020, EAS=0.050 | 0.640 | 0.705 |
| HLA-B_nov_abfaeb71 | novel | 1 | AFR=0.907, EUR=0.010, AMR=0.031, EAS=0.052 | 0.640 | 0.711 |
| HLA-B_nov_fde1221b | novel | 1 | AFR=0.909, EUR=0.000, AMR=0.081, EAS=0.010 | 0.641 | 0.711 |
| 18:20 | known | 1 | AFR=0.000, EUR=0.928, AMR=0.052, EAS=0.021 | 0.647 | 0.705 |
| 18:34 | known | 1 | AFR=0.000, EUR=0.929, AMR=0.061, EAS=0.010 | 0.647 | 0.705 |
| HLA-B_nov_2868079e | novel | 1 | AFR=0.052, EUR=0.021, AMR=0.887, EAS=0.041 | 0.647 | 0.639 |
| 44:06 | known | 1 | AFR=0.000, EUR=0.926, AMR=0.042, EAS=0.032 | 0.648 | 0.705 |
| HLA-B_nov_85fca092 | novel | 1 | AFR=0.000, EUR=0.929, AMR=0.071, EAS=0.000 | 0.648 | 0.000 |
| HLA-B_nov_bf76c177 | novel | 1 | AFR=0.013, EUR=0.924, AMR=0.025, EAS=0.038 | 0.651 | 0.705 |
| HLA-B_nov_df3dd115 | novel | 1 | AFR=0.000, EUR=0.922, AMR=0.044, EAS=0.033 | 0.651 | 0.705 |
| HLA-B_nov_98197dc1 | novel | 1 | AFR=0.000, EUR=0.000, AMR=0.242, EAS=0.758 | 0.651 | 0.121 |
| 35:516 | known | 1 | AFR=0.000, EUR=0.923, AMR=0.055, EAS=0.022 | 0.654 | 0.705 |
| HLA-B_nov_198ccf04 | novel | 1 | AFR=0.000, EUR=0.923, AMR=0.055, EAS=0.022 | 0.654 | 0.705 |
| HLA-B_nov_59632b99 | novel | 1 | AFR=0.000, EUR=0.923, AMR=0.055, EAS=0.022 | 0.654 | 0.705 |
| HLA-B_nov_a46e00a7 | novel | 1 | AFR=0.903, EUR=0.054, AMR=0.032, EAS=0.011 | 0.656 | 0.711 |
| 08:35 | known | 1 | AFR=0.010, EUR=0.920, AMR=0.030, EAS=0.040 | 0.657 | 0.705 |
| HLA-B_nov_99ead9b9 | novel | 1 | AFR=0.000, EUR=0.922, AMR=0.056, EAS=0.022 | 0.657 | 0.705 |
| HLA-B_nov_30cc397d | novel | 1 | AFR=0.000, EUR=0.923, AMR=0.077, EAS=0.000 | 0.657 | 0.705 |
| HLA-B_nov_9455989c | novel | 1 | AFR=0.000, EUR=0.915, AMR=0.011, EAS=0.074 | 0.657 | 0.705 |
| HLA-B_nov_950fe3aa | novel | 1 | AFR=0.050, EUR=0.030, AMR=0.880, EAS=0.040 | 0.659 | 0.639 |
| HLA-B_nov_633f823f | novel | 1 | AFR=0.000, EUR=0.918, AMR=0.035, EAS=0.047 | 0.659 | 0.705 |
| HLA-B_nov_ff49c3b2 | novel | 1 | AFR=0.000, EUR=0.918, AMR=0.035, EAS=0.047 | 0.659 | 0.705 |
| 27:08 | known | 1 | AFR=0.000, EUR=0.120, AMR=0.880, EAS=0.000 | 0.663 | 0.639 |
| HLA-B_nov_69ee4ede | novel | 1 | AFR=0.030, EUR=0.920, AMR=0.030, EAS=0.020 | 0.663 | 0.705 |
| HLA-B_nov_5f30e7ec | novel | 1 | AFR=0.000, EUR=0.912, AMR=0.033, EAS=0.055 | 0.666 | 0.705 |
| HLA-B_nov_90596be6 | novel | 1 | AFR=0.889, EUR=0.000, AMR=0.091, EAS=0.020 | 0.668 | 0.711 |
| HLA-B_nov_60277aac | novel | 1 | AFR=0.031, EUR=0.072, AMR=0.876, EAS=0.021 | 0.668 | 0.639 |
| 14:85N | known | 1 | AFR=0.013, EUR=0.911, AMR=0.051, EAS=0.025 | 0.674 | 0.705 |
| HLA-B_nov_672c5e49 | novel | 1 | AFR=0.000, EUR=0.899, AMR=0.010, EAS=0.091 | 0.674 | 0.705 |
| HLA-B_nov_3444d893 | novel | 1 | AFR=0.111, EUR=0.010, AMR=0.869, EAS=0.010 | 0.674 | 0.639 |
| HLA-B_nov_9a3a8c00 | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.069, EAS=0.023 | 0.677 | 0.705 |
| HLA-B_nov_81b45220 | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.080, EAS=0.011 | 0.677 | 0.705 |
| HLA-B_nov_abbf46cb | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.080, EAS=0.011 | 0.677 | 0.705 |
| HLA-B_nov_d661c981 | novel | 1 | AFR=0.000, EUR=0.908, AMR=0.069, EAS=0.023 | 0.677 | 0.705 |
| HLA-B_nov_0453be48 | novel | 1 | AFR=0.000, EUR=0.901, AMR=0.025, EAS=0.074 | 0.677 | 0.705 |
| HLA-B_nov_2d37ddd9 | novel | 1 | AFR=0.000, EUR=0.905, AMR=0.054, EAS=0.041 | 0.678 | 0.705 |
| HLA-B_nov_3995636d | novel | 1 | AFR=0.000, EUR=0.898, AMR=0.020, EAS=0.082 | 0.679 | 0.705 |
| HLA-B_nov_274e126e | novel | 1 | AFR=0.073, EUR=0.021, AMR=0.865, EAS=0.042 | 0.680 | 0.639 |
| HLA-B_nov_fd593752 | novel | 1 | AFR=0.000, EUR=0.902, AMR=0.049, EAS=0.049 | 0.681 | 0.705 |
| HLA-B_nov_2ec3ff0e | novel | 1 | AFR=0.000, EUR=0.904, AMR=0.085, EAS=0.011 | 0.687 | 0.705 |
| HLA-B_nov_fdc9e59d | novel | 1 | AFR=0.879, EUR=0.000, AMR=0.121, EAS=0.000 | 0.687 | 0.711 |
| HLA-B_nov_4c946edc | novel | 1 | AFR=0.000, EUR=0.897, AMR=0.034, EAS=0.069 | 0.688 | 0.000 |
| HLA-B_nov_67886062 | novel | 1 | AFR=0.000, EUR=0.900, AMR=0.067, EAS=0.033 | 0.689 | 0.705 |
| 44:04 | known | 2 | AFR=0.510, EUR=0.020, AMR=0.455, EAS=0.015 | 0.691 | 0.826 |
| HLA-B_nov_7302f76d | novel | 1 | AFR=0.874, EUR=0.000, AMR=0.126, EAS=0.000 | 0.691 | 0.711 |
| HLA-B_nov_a8e4136a | novel | 1 | AFR=0.000, EUR=0.895, AMR=0.058, EAS=0.047 | 0.691 | 0.705 |
| HLA-B_nov_54cd9fa3 | novel | 1 | AFR=0.000, EUR=0.897, AMR=0.093, EAS=0.010 | 0.692 | 0.705 |
| HLA-B_nov_f50f5dda | novel | 1 | AFR=0.020, EUR=0.890, AMR=0.010, EAS=0.080 | 0.692 | 0.705 |
| HLA-B_nov_d687a2c5 | novel | 1 | AFR=0.870, EUR=0.000, AMR=0.130, EAS=0.000 | 0.694 | 0.711 |
| HLA-B_nov_d9ef93f6 | novel | 1 | AFR=0.137, EUR=0.000, AMR=0.853, EAS=0.011 | 0.694 | 0.639 |
| HLA-B_nov_0ae312f6 | novel | 1 | AFR=0.112, EUR=0.020, AMR=0.857, EAS=0.010 | 0.694 | 0.639 |
| HLA-B_nov_2424728f | novel | 1 | AFR=0.041, EUR=0.072, AMR=0.856, EAS=0.031 | 0.698 | 0.639 |
| HLA-B_nov_997b815e | novel | 1 | AFR=0.084, EUR=0.032, AMR=0.853, EAS=0.032 | 0.699 | 0.639 |
| HLA-B_nov_149a7a9c | novel | 1 | AFR=0.871, EUR=0.065, AMR=0.043, EAS=0.022 | 0.699 | 0.711 |
| HLA-B_nov_c4f7dbd5 | novel | 1 | AFR=0.122, EUR=0.000, AMR=0.847, EAS=0.031 | 0.699 | 0.639 |
| HLA-B_nov_96a3d8c0 | novel | 1 | AFR=0.000, EUR=0.875, AMR=0.000, EAS=0.125 | 0.699 | 0.000 |
| HLA-B_nov_bad7b9fc | novel | 1 | AFR=0.000, EUR=0.889, AMR=0.056, EAS=0.056 | 0.699 | 0.705 |
| HLA-B_nov_e8b0f3b7 | novel | 1 | AFR=0.010, EUR=0.880, AMR=0.010, EAS=0.100 | 0.699 | 0.705 |
| 27:17 | known | 1 | AFR=0.040, EUR=0.890, AMR=0.040, EAS=0.030 | 0.701 | 0.705 |
| HLA-B_nov_c3a4bc1a | novel | 1 | AFR=0.012, EUR=0.884, AMR=0.035, EAS=0.070 | 0.701 | 0.705 |
| HLA-B_nov_7ae61930 | novel | 1 | AFR=0.000, EUR=0.886, AMR=0.068, EAS=0.045 | 0.701 | 0.705 |
| HLA-B_nov_0786781c | novel | 1 | AFR=0.040, EUR=0.890, AMR=0.040, EAS=0.030 | 0.701 | 0.705 |
| HLA-B_nov_d9f13efa | novel | 1 | AFR=0.040, EUR=0.080, AMR=0.850, EAS=0.030 | 0.701 | 0.639 |
| 35:55 | known | 1 | AFR=0.000, EUR=0.884, AMR=0.070, EAS=0.047 | 0.702 | 0.705 |
| HLA-B_nov_e14db612 | novel | 1 | AFR=0.111, EUR=0.030, AMR=0.848, EAS=0.010 | 0.702 | 0.639 |
| HLA-B_nov_e2b93788 | novel | 1 | AFR=0.000, EUR=0.882, AMR=0.059, EAS=0.059 | 0.702 | 0.705 |
| HLA-B_nov_4951f11c | novel | 1 | AFR=0.000, EUR=0.885, AMR=0.090, EAS=0.026 | 0.703 | 0.705 |
| HLA-B_nov_8101def2 | novel | 1 | AFR=0.062, EUR=0.062, AMR=0.845, EAS=0.031 | 0.709 | 0.639 |
| HLA-B_nov_d3ddc74c | novel | 1 | AFR=0.856, EUR=0.033, AMR=0.044, EAS=0.067 | 0.709 | 0.711 |
| HLA-B_nov_bf61ed92 | novel | 1 | AFR=0.000, EUR=0.880, AMR=0.084, EAS=0.036 | 0.709 | 0.705 |
| HLA-B_nov_b0ad9c35 | novel | 1 | AFR=0.010, EUR=0.870, AMR=0.020, EAS=0.100 | 0.709 | 0.705 |
| 18:05 | known | 1 | AFR=0.000, EUR=0.130, AMR=0.840, EAS=0.030 | 0.711 | 0.639 |
| HLA-B_nov_d12e5487 | novel | 1 | AFR=0.030, EUR=0.880, AMR=0.050, EAS=0.040 | 0.712 | 0.705 |
| HLA-B_nov_7dc0e55c | novel | 1 | AFR=0.000, EUR=0.876, AMR=0.082, EAS=0.041 | 0.713 | 0.705 |
| HLA-B_nov_bfb6333b | novel | 1 | AFR=0.000, EUR=0.877, AMR=0.096, EAS=0.027 | 0.714 | 0.705 |
| HLA-B_nov_c8196e09 | novel | 1 | AFR=0.000, EUR=0.874, AMR=0.069, EAS=0.057 | 0.716 | 0.705 |
| HLA-B_nov_c2d85fff | novel | 1 | AFR=0.020, EUR=0.101, AMR=0.838, EAS=0.040 | 0.716 | 0.639 |
| HLA-B_nov_88cce230 | novel | 1 | AFR=0.000, EUR=0.875, AMR=0.104, EAS=0.021 | 0.718 | 0.705 |
| HLA-B_nov_8738690d | novel | 1 | AFR=0.000, EUR=0.872, AMR=0.081, EAS=0.047 | 0.719 | 0.705 |
| HLA-B_nov_0213941a | novel | 1 | AFR=0.000, EUR=0.872, AMR=0.085, EAS=0.043 | 0.720 | 0.705 |
| HLA-B_nov_64107ce0 | novel | 1 | AFR=0.000, EUR=0.872, AMR=0.090, EAS=0.038 | 0.721 | 0.705 |
| HLA-B_nov_c31c522f | novel | 1 | AFR=0.011, EUR=0.875, AMR=0.091, EAS=0.023 | 0.721 | 0.705 |
| HLA-B_nov_970e6c4a | novel | 1 | AFR=0.000, EUR=0.869, AMR=0.107, EAS=0.024 | 0.723 | 0.705 |
| HLA-B_nov_dfac340c | novel | 1 | AFR=0.012, EUR=0.867, AMR=0.072, EAS=0.048 | 0.724 | 0.705 |
| HLA-B_nov_8a010d86 | novel | 1 | AFR=0.013, EUR=0.867, AMR=0.067, EAS=0.053 | 0.726 | 0.705 |
| HLA-B_nov_d63ce2a9 | novel | 1 | AFR=0.017, EUR=0.867, AMR=0.067, EAS=0.050 | 0.726 | 0.705 |
| HLA-B_nov_a97e900d | novel | 1 | AFR=0.000, EUR=0.862, AMR=0.080, EAS=0.057 | 0.728 | 0.705 |
| HLA-B_nov_5971e459 | novel | 1 | AFR=0.000, EUR=0.862, AMR=0.088, EAS=0.050 | 0.731 | 0.705 |
| HLA-B_nov_c30cf4d2 | novel | 2 | AFR=0.151, EUR=0.108, AMR=0.613, EAS=0.129 | 0.732 | 0.410 |
| HLA-B_nov_fa264b30 | novel | 1 | AFR=0.000, EUR=0.859, AMR=0.082, EAS=0.059 | 0.736 | 0.705 |
| HLA-B_nov_67b2eeab | novel | 1 | AFR=0.062, EUR=0.072, AMR=0.825, EAS=0.041 | 0.737 | 0.639 |
| HLA-B_nov_160154e6 | novel | 1 | AFR=0.012, EUR=0.859, AMR=0.071, EAS=0.059 | 0.737 | 0.705 |
| HLA-B_nov_d73ebc0c | novel | 1 | AFR=0.163, EUR=0.010, AMR=0.816, EAS=0.010 | 0.739 | 0.639 |
| HLA-B_nov_69b8d710 | novel | 1 | AFR=0.020, EUR=0.110, AMR=0.820, EAS=0.050 | 0.739 | 0.639 |
| 27:09 | known | 1 | AFR=0.000, EUR=0.857, AMR=0.143, EAS=0.000 | 0.742 | 0.000 |
| HLA-B_nov_f46291ae | novel | 1 | AFR=0.830, EUR=0.000, AMR=0.160, EAS=0.011 | 0.742 | 0.711 |
| HLA-B_nov_89cf2b08 | novel | 1 | AFR=0.000, EUR=0.856, AMR=0.089, EAS=0.056 | 0.742 | 0.705 |
| HLA-B_nov_98c711e5 | novel | 1 | AFR=0.833, EUR=0.031, AMR=0.094, EAS=0.042 | 0.744 | 0.711 |
| HLA-B_nov_d48371e1 | novel | 1 | AFR=0.827, EUR=0.010, AMR=0.163, EAS=0.000 | 0.744 | 0.711 |
| HLA-B_nov_470c7b6a | novel | 1 | AFR=0.010, EUR=0.153, AMR=0.816, EAS=0.020 | 0.744 | 0.639 |
| HLA-B_nov_f718667b | novel | 1 | AFR=0.000, EUR=0.846, AMR=0.077, EAS=0.077 | 0.748 | 0.000 |
| HLA-B_nov_ebe632a1 | novel | 1 | AFR=0.013, EUR=0.848, AMR=0.089, EAS=0.051 | 0.751 | 0.705 |
| HLA-B_nov_903af24a | novel | 1 | AFR=0.000, EUR=0.848, AMR=0.109, EAS=0.043 | 0.751 | 0.000 |
| HLA-B_nov_9914427e | novel | 1 | AFR=0.010, EUR=0.160, AMR=0.810, EAS=0.020 | 0.751 | 0.639 |
| HLA-B_nov_6209ee87 | novel | 1 | AFR=0.010, EUR=0.160, AMR=0.810, EAS=0.020 | 0.751 | 0.639 |
| HLA-B_nov_560d3c37 | novel | 1 | AFR=0.000, EUR=0.845, AMR=0.113, EAS=0.041 | 0.752 | 0.705 |
| HLA-B_nov_4f658252 | novel | 1 | AFR=0.143, EUR=0.020, AMR=0.806, EAS=0.031 | 0.752 | 0.639 |
| HLA-B_nov_abeace88 | novel | 1 | AFR=0.823, EUR=0.021, AMR=0.115, EAS=0.042 | 0.752 | 0.711 |
| HLA-B_nov_3cf50b0e | novel | 1 | AFR=0.818, EUR=0.000, AMR=0.162, EAS=0.020 | 0.752 | 0.711 |
| HLA-B_nov_893848f1 | novel | 1 | AFR=0.013, EUR=0.842, AMR=0.079, EAS=0.066 | 0.753 | 0.705 |
| HLA-B_nov_257a6409 | novel | 1 | AFR=0.020, EUR=0.162, AMR=0.808, EAS=0.010 | 0.753 | 0.639 |
| HLA-B_nov_54da8a2c | novel | 1 | AFR=0.012, EUR=0.841, AMR=0.073, EAS=0.073 | 0.753 | 0.705 |
| HLA-B_nov_94a057b6 | novel | 1 | AFR=0.014, EUR=0.843, AMR=0.086, EAS=0.057 | 0.756 | 0.705 |
| HLA-B_nov_d8da2570 | novel | 1 | AFR=0.194, EUR=0.010, AMR=0.796, EAS=0.000 | 0.758 | 0.639 |
| HLA-B_nov_2423c642 | novel | 1 | AFR=0.000, EUR=0.837, AMR=0.105, EAS=0.058 | 0.761 | 0.705 |
| HLA-B_nov_a023356a | novel | 1 | AFR=0.012, EUR=0.841, AMR=0.122, EAS=0.024 | 0.761 | 0.705 |
| HLA-B_nov_35f7e815 | novel | 1 | AFR=0.000, EUR=0.839, AMR=0.126, EAS=0.034 | 0.761 | 0.705 |
| HLA-B_nov_70d6526e | novel | 2 | AFR=0.030, EUR=0.625, AMR=0.230, EAS=0.115 | 0.762 | 0.823 |
| HLA-B_nov_009d823b | novel | 1 | AFR=0.011, EUR=0.839, AMR=0.103, EAS=0.046 | 0.763 | 0.705 |
| HLA-B_nov_708e0c17 | novel | 1 | AFR=0.814, EUR=0.041, AMR=0.082, EAS=0.062 | 0.766 | 0.711 |
| HLA-B_nov_9479734a | novel | 1 | AFR=0.143, EUR=0.031, AMR=0.796, EAS=0.031 | 0.771 | 0.639 |
| HLA-B_nov_f4e609fe | novel | 1 | AFR=0.000, EUR=0.830, AMR=0.102, EAS=0.068 | 0.772 | 0.705 |
| HLA-B_nov_ca434506 | novel | 1 | AFR=0.000, EUR=0.831, AMR=0.130, EAS=0.039 | 0.772 | 0.705 |
| HLA-B_nov_e258e8f9 | novel | 1 | AFR=0.000, EUR=0.830, AMR=0.114, EAS=0.057 | 0.773 | 0.705 |
| 38:91 | known | 1 | AFR=0.013, EUR=0.831, AMR=0.104, EAS=0.052 | 0.776 | 0.705 |
| HLA-B_nov_ce78bc9f | novel | 1 | AFR=0.000, EUR=0.827, AMR=0.123, EAS=0.049 | 0.777 | 0.705 |
| HLA-B_nov_cfb772eb | novel | 1 | AFR=0.011, EUR=0.830, AMR=0.125, EAS=0.034 | 0.778 | 0.705 |
| HLA-B_nov_a0e01a96 | novel | 1 | AFR=0.811, EUR=0.056, AMR=0.111, EAS=0.022 | 0.782 | 0.711 |
| HLA-B_nov_a80fc10b | novel | 1 | AFR=0.000, EUR=0.821, AMR=0.095, EAS=0.083 | 0.786 | 0.705 |
| HLA-B_nov_9fbe8078 | novel | 1 | AFR=0.010, EUR=0.806, AMR=0.031, EAS=0.153 | 0.787 | 0.705 |
| HLA-B_nov_3c32a88a | novel | 1 | AFR=0.214, EUR=0.010, AMR=0.776, EAS=0.000 | 0.792 | 0.639 |
| HLA-B_nov_6f459aa7 | novel | 1 | AFR=0.000, EUR=0.820, AMR=0.135, EAS=0.045 | 0.793 | 0.705 |
| HLA-B_nov_bbbd239d | novel | 1 | AFR=0.000, EUR=0.820, AMR=0.135, EAS=0.045 | 0.793 | 0.705 |
| HLA-B_nov_ff663f66 | novel | 1 | AFR=0.796, EUR=0.022, AMR=0.108, EAS=0.075 | 0.793 | 0.711 |
| HLA-B_nov_2d7e9819 | novel | 1 | AFR=0.786, EUR=0.000, AMR=0.214, EAS=0.000 | 0.794 | 0.711 |
| HLA-B_nov_05c6957c | novel | 1 | AFR=0.216, EUR=0.010, AMR=0.773, EAS=0.000 | 0.794 | 0.639 |
| HLA-B_nov_7ab154bd | novel | 1 | AFR=0.000, EUR=0.818, AMR=0.136, EAS=0.045 | 0.794 | 0.705 |
| HLA-B_nov_2920fc70 | novel | 1 | AFR=0.000, EUR=0.814, AMR=0.116, EAS=0.070 | 0.796 | 0.705 |
| HLA-B_nov_eef5bdfe | novel | 1 | AFR=0.000, EUR=0.815, AMR=0.123, EAS=0.062 | 0.796 | 0.705 |
| HLA-B_nov_10d36ffc | novel | 1 | AFR=0.162, EUR=0.030, AMR=0.778, EAS=0.030 | 0.797 | 0.639 |
| HLA-B_nov_75d764e9 | novel | 1 | AFR=0.010, EUR=0.814, AMR=0.155, EAS=0.021 | 0.800 | 0.705 |
| HLA-B_nov_f8a83d5c | novel | 1 | AFR=0.173, EUR=0.000, AMR=0.765, EAS=0.061 | 0.800 | 0.639 |
| HLA-B_nov_d4b14b9a | novel | 1 | AFR=0.153, EUR=0.031, AMR=0.776, EAS=0.041 | 0.800 | 0.639 |
| HLA-B_nov_dd140996 | novel | 1 | AFR=0.012, EUR=0.812, AMR=0.141, EAS=0.035 | 0.802 | 0.705 |
| HLA-B_nov_75c0c1a4 | novel | 1 | AFR=0.030, EUR=0.810, AMR=0.080, EAS=0.080 | 0.802 | 0.705 |
| HLA-B_nov_b6619d8d | novel | 1 | AFR=0.012, EUR=0.810, AMR=0.119, EAS=0.060 | 0.803 | 0.705 |
| HLA-B_nov_2783493a | novel | 1 | AFR=0.000, EUR=0.807, AMR=0.148, EAS=0.045 | 0.804 | 0.705 |
| HLA-B_nov_f0ecc128 | novel | 1 | AFR=0.770, EUR=0.000, AMR=0.230, EAS=0.000 | 0.806 | 0.711 |
| HLA-B_nov_0a3ca459 | novel | 1 | AFR=0.082, EUR=0.082, AMR=0.773, EAS=0.062 | 0.806 | 0.639 |
| HLA-B_nov_c567a698 | novel | 1 | AFR=0.000, EUR=0.750, AMR=0.000, EAS=0.250 | 0.806 | 0.000 |
| HLA-B_nov_01bdf7b0 | novel | 1 | AFR=0.000, EUR=0.750, AMR=0.000, EAS=0.250 | 0.806 | 0.000 |
| HLA-B_nov_80c40461 | novel | 1 | AFR=0.790, EUR=0.086, AMR=0.086, EAS=0.037 | 0.810 | 0.711 |
| HLA-B_nov_3f16cac8 | novel | 1 | AFR=0.000, EUR=0.795, AMR=0.120, EAS=0.084 | 0.811 | 0.705 |
| HLA-B_nov_f05a082c | novel | 1 | AFR=0.000, EUR=0.798, AMR=0.131, EAS=0.071 | 0.811 | 0.705 |
| HLA-B_nov_5f4e9a68 | novel | 4 | AFR=0.258, EUR=0.495, AMR=0.216, EAS=0.031 | 0.814 | 0.842 |
| HLA-B_nov_a5806dd7 | novel | 1 | AFR=0.200, EUR=0.022, AMR=0.756, EAS=0.022 | 0.814 | 0.639 |
| HLA-B_nov_41499f6d | novel | 1 | AFR=0.000, EUR=0.795, AMR=0.157, EAS=0.048 | 0.816 | 0.705 |
| HLA-B_nov_f29be7f4 | novel | 1 | AFR=0.000, EUR=0.789, AMR=0.132, EAS=0.079 | 0.826 | 0.705 |
| HLA-B_nov_d5cd981a | novel | 2 | AFR=0.145, EUR=0.626, AMR=0.212, EAS=0.016 | 0.828 | 0.000 |
| HLA-B_nov_fb8dafaa | novel | 1 | AFR=0.270, EUR=0.000, AMR=0.730, EAS=0.000 | 0.829 | 0.639 |
| HLA-B_nov_5a5869cb | novel | 1 | AFR=0.753, EUR=0.010, AMR=0.237, EAS=0.000 | 0.829 | 0.711 |
| HLA-B_nov_62fbec54 | novel | 1 | AFR=0.206, EUR=0.010, AMR=0.742, EAS=0.041 | 0.829 | 0.639 |
| HLA-B_nov_983d840a | novel | 1 | AFR=0.012, EUR=0.238, AMR=0.107, EAS=0.643 | 0.829 | 0.208 |
| HLA-B_nov_f83c8686 | novel | 1 | AFR=0.163, EUR=0.010, AMR=0.745, EAS=0.082 | 0.829 | 0.639 |
| HLA-B_nov_f02e6ac0 | novel | 1 | AFR=0.000, EUR=0.218, AMR=0.141, EAS=0.641 | 0.829 | 0.208 |
| HLA-B_nov_e7abe413 | novel | 1 | AFR=0.753, EUR=0.010, AMR=0.227, EAS=0.010 | 0.830 | 0.711 |
| HLA-B_nov_38e03e62 | novel | 1 | AFR=0.012, EUR=0.787, AMR=0.163, EAS=0.037 | 0.830 | 0.705 |
| HLA-B_nov_03cd7f7f | novel | 1 | AFR=0.000, EUR=0.780, AMR=0.134, EAS=0.085 | 0.834 | 0.705 |
| HLA-B_nov_2fdaced0 | novel | 1 | AFR=0.755, EUR=0.043, AMR=0.085, EAS=0.117 | 0.837 | 0.711 |
| HLA-B_nov_3afdde33 | novel | 1 | AFR=0.283, EUR=0.000, AMR=0.717, EAS=0.000 | 0.837 | 0.639 |
| HLA-B_nov_a612b678 | novel | 1 | AFR=0.760, EUR=0.062, AMR=0.115, EAS=0.062 | 0.837 | 0.711 |
| HLA-B_nov_d25ed3df | novel | 1 | AFR=0.760, EUR=0.062, AMR=0.115, EAS=0.062 | 0.837 | 0.711 |
| HLA-B_nov_bd6c67aa | novel | 1 | AFR=0.755, EUR=0.117, AMR=0.032, EAS=0.096 | 0.839 | 0.711 |
| HLA-B_nov_9ee0b8f8 | novel | 2 | AFR=0.490, EUR=0.455, AMR=0.030, EAS=0.025 | 0.840 | 0.856 |
| HLA-B_nov_6de6f05c | novel | 1 | AFR=0.010, EUR=0.768, AMR=0.212, EAS=0.010 | 0.840 | 0.705 |
| HLA-B_nov_f62b3ed6 | novel | 1 | AFR=0.021, EUR=0.216, AMR=0.732, EAS=0.031 | 0.841 | 0.639 |
| HLA-B_nov_451add33 | novel | 1 | AFR=0.000, EUR=0.765, AMR=0.176, EAS=0.059 | 0.842 | 0.000 |
| HLA-B_nov_e81ee869 | novel | 1 | AFR=0.013, EUR=0.763, AMR=0.125, EAS=0.100 | 0.846 | 0.705 |
| HLA-B_nov_7409cfcb | novel | 1 | AFR=0.013, EUR=0.766, AMR=0.156, EAS=0.065 | 0.846 | 0.705 |
| HLA-B_nov_27431add | novel | 1 | AFR=0.000, EUR=0.500, AMR=0.000, EAS=0.500 | 0.849 | 0.000 |
| HLA-B_nov_587a4cc2 | novel | 1 | AFR=0.253, EUR=0.010, AMR=0.707, EAS=0.030 | 0.851 | 0.639 |
| HLA-B_nov_d39c0c1c | novel | 1 | AFR=0.043, EUR=0.159, AMR=0.174, EAS=0.623 | 0.851 | 0.208 |
| HLA-B_nov_3ac46116 | novel | 1 | AFR=0.000, EUR=0.742, AMR=0.097, EAS=0.161 | 0.851 | 0.705 |
| HLA-B_nov_c734a2f0 | novel | 1 | AFR=0.244, EUR=0.000, AMR=0.700, EAS=0.056 | 0.852 | 0.639 |
| HLA-B_nov_5ca36281 | novel | 1 | AFR=0.097, EUR=0.075, AMR=0.720, EAS=0.108 | 0.852 | 0.639 |
| HLA-B_nov_c4efb5c4 | novel | 1 | AFR=0.219, EUR=0.021, AMR=0.708, EAS=0.052 | 0.852 | 0.639 |
| HLA-B_nov_1745f037 | novel | 1 | AFR=0.030, EUR=0.758, AMR=0.182, EAS=0.030 | 0.852 | 0.121 |
| HLA-B_nov_52cc5c74 | novel | 1 | AFR=0.000, EUR=0.295, AMR=0.115, EAS=0.590 | 0.854 | 0.208 |
| HLA-B_nov_fb8ee4e5 | novel | 1 | AFR=0.684, EUR=0.000, AMR=0.316, EAS=0.000 | 0.858 | 0.711 |
| HLA-B_nov_a4143085 | novel | 1 | AFR=0.684, EUR=0.000, AMR=0.316, EAS=0.000 | 0.858 | 0.711 |
| HLA-B_nov_21f78371 | novel | 1 | AFR=0.287, EUR=0.011, AMR=0.681, EAS=0.021 | 0.859 | 0.639 |
| HLA-B_nov_88b5e881 | novel | 1 | AFR=0.316, EUR=0.010, AMR=0.673, EAS=0.000 | 0.859 | 0.639 |
| HLA-B_nov_f924859c | novel | 1 | AFR=0.095, EUR=0.116, AMR=0.705, EAS=0.084 | 0.861 | 0.639 |
| HLA-B_nov_5600bf31 | novel | 1 | AFR=0.711, EUR=0.041, AMR=0.216, EAS=0.031 | 0.861 | 0.711 |
| HLA-B_nov_4fa3d83a | novel | 1 | AFR=0.083, EUR=0.750, AMR=0.167, EAS=0.000 | 0.861 | 0.000 |
| HLA-B_nov_736a8631 | novel | 1 | AFR=0.041, EUR=0.206, AMR=0.701, EAS=0.052 | 0.862 | 0.639 |
| HLA-B_nov_048d0d89 | novel | 1 | AFR=0.691, EUR=0.021, AMR=0.255, EAS=0.032 | 0.864 | 0.711 |
| 08:09 | known | 1 | AFR=0.330, EUR=0.000, AMR=0.649, EAS=0.021 | 0.866 | 0.639 |
| HLA-B_nov_29d1f9f7 | novel | 1 | AFR=0.650, EUR=0.000, AMR=0.350, EAS=0.000 | 0.866 | 0.711 |
| HLA-B_nov_4d67a1e9 | novel | 1 | AFR=0.012, EUR=0.312, AMR=0.112, EAS=0.562 | 0.868 | 0.208 |
| HLA-B_nov_97d4c44c | novel | 1 | AFR=0.000, EUR=0.700, AMR=0.280, EAS=0.020 | 0.869 | 0.705 |
| HLA-B_nov_978b3f67 | novel | 1 | AFR=0.640, EUR=0.000, AMR=0.360, EAS=0.000 | 0.869 | 0.711 |
| HLA-B_nov_6fe23a07 | novel | 1 | AFR=0.012, EUR=0.310, AMR=0.119, EAS=0.560 | 0.870 | 0.208 |
| HLA-B_nov_beff239d | novel | 1 | AFR=0.042, EUR=0.663, AMR=0.032, EAS=0.263 | 0.870 | 0.705 |
| HLA-B_nov_f8b8518b | novel | 1 | AFR=0.713, EUR=0.138, AMR=0.113, EAS=0.037 | 0.870 | 0.711 |
| 55:97N | known | 1 | AFR=0.090, EUR=0.180, AMR=0.690, EAS=0.040 | 0.872 | 0.639 |
| HLA-B_nov_88f0f010 | novel | 1 | AFR=0.102, EUR=0.143, AMR=0.684, EAS=0.071 | 0.877 | 0.639 |
| HLA-B_nov_10aabe32 | novel | 1 | AFR=0.168, EUR=0.063, AMR=0.674, EAS=0.095 | 0.877 | 0.639 |
| HLA-B_nov_2ef23553 | novel | 1 | AFR=0.000, EUR=0.667, AMR=0.333, EAS=0.000 | 0.877 | 0.000 |
| HLA-B_nov_1db89393 | novel | 1 | AFR=0.384, EUR=0.010, AMR=0.606, EAS=0.000 | 0.882 | 0.639 |
| HLA-B_nov_80d1d7b0 | novel | 1 | AFR=0.385, EUR=0.010, AMR=0.604, EAS=0.000 | 0.882 | 0.639 |
| HLA-B_nov_b757e4ca | novel | 1 | AFR=0.021, EUR=0.632, AMR=0.053, EAS=0.295 | 0.882 | 0.705 |
| HLA-B_nov_723315d8 | novel | 1 | AFR=0.420, EUR=0.000, AMR=0.580, EAS=0.000 | 0.884 | 0.639 |
| HLA-B_nov_fc9e0671 | novel | 1 | AFR=0.596, EUR=0.000, AMR=0.404, EAS=0.000 | 0.884 | 0.711 |
| 18:13 | known | 1 | AFR=0.155, EUR=0.062, AMR=0.660, EAS=0.124 | 0.887 | 0.639 |
| HLA-B_nov_bf18667a | novel | 1 | AFR=0.449, EUR=0.000, AMR=0.551, EAS=0.000 | 0.889 | 0.639 |
| 35:49 | known | 1 | AFR=0.455, EUR=0.000, AMR=0.545, EAS=0.000 | 0.890 | 0.639 |
| HLA-B_nov_b1992d3b | novel | 1 | AFR=0.021, EUR=0.436, AMR=0.064, EAS=0.479 | 0.890 | 0.208 |
| HLA-B_nov_2035a20d | novel | 1 | AFR=0.052, EUR=0.660, AMR=0.062, EAS=0.227 | 0.890 | 0.705 |
| HLA-B_nov_608ea697 | novel | 1 | AFR=0.421, EUR=0.000, AMR=0.568, EAS=0.011 | 0.890 | 0.639 |
| HLA-B_nov_0087a08e | novel | 1 | AFR=0.042, EUR=0.084, AMR=0.579, EAS=0.295 | 0.893 | 0.639 |
| HLA-B_nov_8e6f67ac | novel | 1 | AFR=0.510, EUR=0.000, AMR=0.490, EAS=0.000 | 0.894 | 0.711 |
| HLA-B_nov_aa4fcd4b | novel | 1 | AFR=0.202, EUR=0.053, AMR=0.649, EAS=0.096 | 0.894 | 0.639 |
| HLA-B_nov_95cbb061 | novel | 1 | AFR=0.033, EUR=0.370, AMR=0.087, EAS=0.511 | 0.894 | 0.208 |
| HLA-B_nov_8cc7819b | novel | 1 | AFR=0.663, EUR=0.043, AMR=0.196, EAS=0.098 | 0.894 | 0.711 |
| HLA-B_nov_530c73c3 | novel | 1 | AFR=0.022, EUR=0.333, AMR=0.122, EAS=0.522 | 0.896 | 0.208 |
| HLA-B_nov_d30b9a84 | novel | 1 | AFR=0.124, EUR=0.093, AMR=0.639, EAS=0.144 | 0.898 | 0.639 |
| HLA-B_nov_42efe47f | novel | 1 | AFR=0.000, EUR=0.439, AMR=0.561, EAS=0.000 | 0.898 | 0.639 |
| 38:09 | known | 2 | AFR=0.465, EUR=0.448, AMR=0.065, EAS=0.022 | 0.899 | 0.856 |
| HLA-B_nov_3b05e11b | novel | 1 | AFR=0.663, EUR=0.109, AMR=0.185, EAS=0.043 | 0.904 | 0.711 |
| 53:05 | known | 1 | AFR=0.000, EUR=0.556, AMR=0.444, EAS=0.000 | 0.907 | 0.000 |
| HLA-B_nov_6e4db735 | novel | 1 | AFR=0.128, EUR=0.074, AMR=0.617, EAS=0.181 | 0.907 | 0.639 |
| HLA-B_nov_68a9e512 | novel | 1 | AFR=0.021, EUR=0.371, AMR=0.588, EAS=0.021 | 0.910 | 0.639 |
| HLA-B_nov_2516c759 | novel | 1 | AFR=0.098, EUR=0.174, AMR=0.641, EAS=0.087 | 0.910 | 0.639 |
| HLA-B_nov_9080fd27 | novel | 1 | AFR=0.011, EUR=0.479, AMR=0.096, EAS=0.415 | 0.910 | 0.705 |
| HLA-B_nov_e609b2a7 | novel | 1 | AFR=0.011, EUR=0.479, AMR=0.096, EAS=0.415 | 0.910 | 0.705 |
| HLA-B_nov_8b9103a4 | novel | 1 | AFR=0.000, EUR=0.571, AMR=0.129, EAS=0.300 | 0.916 | 0.705 |
| HLA-B_nov_fd264f91 | novel | 1 | AFR=0.371, EUR=0.021, AMR=0.557, EAS=0.052 | 0.917 | 0.639 |
| HLA-B_nov_5e92b57a | novel | 1 | AFR=0.062, EUR=0.646, AMR=0.156, EAS=0.135 | 0.919 | 0.705 |
| HLA-B_nov_feaa531e | novel | 1 | AFR=0.022, EUR=0.630, AMR=0.217, EAS=0.130 | 0.919 | 0.121 |
| HLA-B_nov_4379b659 | novel | 1 | AFR=0.091, EUR=0.202, AMR=0.616, EAS=0.091 | 0.924 | 0.639 |
| HLA-B_nov_cb221017 | novel | 1 | AFR=0.010, EUR=0.485, AMR=0.485, EAS=0.020 | 0.924 | 0.639 |
| HLA-B_nov_5daeb7bc | novel | 1 | AFR=0.032, EUR=0.226, AMR=0.258, EAS=0.484 | 0.924 | 0.121 |
| HLA-B_nov_915a65b7 | novel | 1 | AFR=0.091, EUR=0.636, AMR=0.273, EAS=0.000 | 0.926 | 0.000 |
| HLA-B_nov_3974baa8 | novel | 1 | AFR=0.630, EUR=0.111, AMR=0.185, EAS=0.074 | 0.927 | 0.711 |
| HLA-B_nov_121201a0 | novel | 1 | AFR=0.117, EUR=0.181, AMR=0.617, EAS=0.085 | 0.927 | 0.639 |
| HLA-B_nov_e6962c16 | novel | 1 | AFR=0.032, EUR=0.258, AMR=0.226, EAS=0.484 | 0.931 | 0.121 |
| HLA-B_nov_4a3fccf1 | novel | 1 | AFR=0.630, EUR=0.174, AMR=0.109, EAS=0.087 | 0.932 | 0.711 |
| HLA-B_nov_d082bee6 | novel | 1 | AFR=0.074, EUR=0.632, AMR=0.158, EAS=0.137 | 0.937 | 0.705 |
| HLA-B_nov_3e25c4fd | novel | 1 | AFR=0.032, EUR=0.442, AMR=0.126, EAS=0.400 | 0.940 | 0.705 |
| HLA-B_nov_32ecffa1 | novel | 1 | AFR=0.383, EUR=0.025, AMR=0.185, EAS=0.407 | 0.943 | 0.208 |
| HLA-B_nov_e8dacf34 | novel | 1 | AFR=0.609, EUR=0.161, AMR=0.103, EAS=0.126 | 0.944 | 0.711 |
| 51:14 | known | 1 | AFR=0.082, EUR=0.235, AMR=0.582, EAS=0.102 | 0.946 | 0.639 |
| HLA-B_nov_72abb4d5 | novel | 1 | AFR=0.531, EUR=0.062, AMR=0.375, EAS=0.031 | 0.946 | 0.711 |
| HLA-B_nov_128bc54d | novel | 1 | AFR=0.596, EUR=0.135, AMR=0.112, EAS=0.157 | 0.946 | 0.711 |
| HLA-B_nov_9ad05e51 | novel | 1 | AFR=0.030, EUR=0.530, AMR=0.400, EAS=0.040 | 0.946 | 0.705 |
| 40:139 | known | 1 | AFR=0.600, EUR=0.189, AMR=0.111, EAS=0.100 | 0.949 | 0.711 |
| HLA-B_nov_d9fdfca4 | novel | 1 | AFR=0.600, EUR=0.189, AMR=0.111, EAS=0.100 | 0.949 | 0.711 |
| HLA-B_nov_6d041ba3 | novel | 1 | AFR=0.402, EUR=0.052, AMR=0.485, EAS=0.062 | 0.950 | 0.639 |
| HLA-B_nov_e368bf67 | novel | 1 | AFR=0.589, EUR=0.167, AMR=0.156, EAS=0.089 | 0.951 | 0.711 |
| HLA-B_nov_3a91837f | novel | 1 | AFR=0.141, EUR=0.172, AMR=0.566, EAS=0.121 | 0.951 | 0.639 |
| 15:34 | known | 2 | AFR=0.096, EUR=0.490, AMR=0.348, EAS=0.065 | 0.952 | 0.823 |
| HLA-B_nov_f1cd89ab | novel | 1 | AFR=0.043, EUR=0.362, AMR=0.191, EAS=0.404 | 0.952 | 0.121 |
| HLA-B_nov_3654bb38 | novel | 1 | AFR=0.012, EUR=0.488, AMR=0.221, EAS=0.279 | 0.953 | 0.705 |
| HLA-B_nov_124cb7d8 | novel | 1 | AFR=0.465, EUR=0.070, AMR=0.163, EAS=0.302 | 0.958 | 0.121 |
| HLA-B_nov_13b4e364 | novel | 1 | AFR=0.040, EUR=0.444, AMR=0.444, EAS=0.071 | 0.958 | 0.639 |
| HLA-B_nov_3d00e477 | novel | 1 | AFR=0.040, EUR=0.444, AMR=0.444, EAS=0.071 | 0.958 | 0.639 |
| 27:14 | known | 1 | AFR=0.111, EUR=0.283, AMR=0.535, EAS=0.071 | 0.960 | 0.639 |
| HLA-B_nov_5947d4af | novel | 1 | AFR=0.556, EUR=0.133, AMR=0.211, EAS=0.100 | 0.960 | 0.711 |
| 27:10 | known | 1 | AFR=0.061, EUR=0.323, AMR=0.505, EAS=0.111 | 0.961 | 0.639 |
| HLA-B_nov_7cdbfe7c | novel | 1 | AFR=0.104, EUR=0.552, AMR=0.146, EAS=0.198 | 0.961 | 0.705 |
| HLA-B_nov_4e2d3568 | novel | 1 | AFR=0.562, EUR=0.180, AMR=0.169, EAS=0.090 | 0.962 | 0.711 |
| HLA-B_nov_34f5e896 | novel | 1 | AFR=0.562, EUR=0.180, AMR=0.169, EAS=0.090 | 0.962 | 0.711 |
| HLA-B_nov_b6c41ddd | novel | 1 | AFR=0.320, EUR=0.040, AMR=0.320, EAS=0.320 | 0.964 | 0.121 |
| HLA-B_nov_ae080ac9 | novel | 1 | AFR=0.060, EUR=0.360, AMR=0.480, EAS=0.100 | 0.964 | 0.639 |
| HLA-B_nov_ae76889a | novel | 1 | AFR=0.070, EUR=0.400, AMR=0.450, EAS=0.080 | 0.966 | 0.639 |
| HLA-B_nov_8c596552 | novel | 1 | AFR=0.131, EUR=0.182, AMR=0.495, EAS=0.192 | 0.967 | 0.639 |
| HLA-B_nov_1c72f545 | novel | 1 | AFR=0.083, EUR=0.500, AMR=0.208, EAS=0.208 | 0.969 | 0.000 |
| HLA-B_nov_ceb2f866 | novel | 1 | AFR=0.442, EUR=0.140, AMR=0.174, EAS=0.244 | 0.976 | 0.711 |
| HLA-B_nov_bf505b36 | novel | 1 | AFR=0.110, EUR=0.407, AMR=0.187, EAS=0.297 | 0.977 | 0.705 |
| HLA-B_nov_1db2a33a | novel | 1 | AFR=0.146, EUR=0.250, AMR=0.417, EAS=0.188 | 0.982 | 0.639 |
| HLA-B_nov_087e3604 | novel | 1 | AFR=0.427, EUR=0.247, AMR=0.124, EAS=0.202 | 0.982 | 0.711 |
| HLA-B_nov_164a08f8 | novel | 1 | AFR=0.114, EUR=0.352, AMR=0.250, EAS=0.284 | 0.982 | 0.705 |
| HLA-B_nov_4fac5a8f | novel | 1 | AFR=0.449, EUR=0.180, AMR=0.180, EAS=0.191 | 0.982 | 0.711 |
| 35:205 | known | 1 | AFR=0.354, EUR=0.125, AMR=0.312, EAS=0.208 | 0.983 | 0.121 |
| HLA-B_nov_cc11ebef | novel | 1 | AFR=0.218, EUR=0.322, AMR=0.149, EAS=0.310 | 0.983 | 0.705 |
| HLA-B_nov_6b37665d | novel | 1 | AFR=0.414, EUR=0.276, AMR=0.126, EAS=0.184 | 0.983 | 0.711 |
| HLA-B_nov_81116cfa | novel | 1 | AFR=0.302, EUR=0.140, AMR=0.349, EAS=0.209 | 0.984 | 0.639 |
| HLA-B_nov_0f47e2cd | novel | 1 | AFR=0.302, EUR=0.140, AMR=0.349, EAS=0.209 | 0.984 | 0.639 |
| HLA-B_nov_7d093379 | novel | 1 | AFR=0.267, EUR=0.444, AMR=0.178, EAS=0.111 | 0.987 | 0.705 |
| HLA-B_nov_8877fa30 | novel | 1 | AFR=0.398, EUR=0.273, AMR=0.148, EAS=0.182 | 0.987 | 0.711 |
| HLA-B_nov_b7d1f5da | novel | 1 | AFR=0.370, EUR=0.321, AMR=0.136, EAS=0.173 | 0.989 | 0.711 |
| HLA-B_nov_877845eb | novel | 1 | AFR=0.370, EUR=0.348, AMR=0.141, EAS=0.141 | 0.989 | 0.711 |
| HLA-B_nov_e3a8b2e1 | novel | 1 | AFR=0.170, EUR=0.394, AMR=0.266, EAS=0.170 | 0.990 | 0.705 |
| HLA-B_nov_35865647 | novel | 1 | AFR=0.227, EUR=0.239, AMR=0.250, EAS=0.284 | 0.990 | 0.208 |
| HLA-B_nov_c012ebe3 | novel | 1 | AFR=0.385, EUR=0.297, AMR=0.209, EAS=0.110 | 0.990 | 0.711 |
| HLA-B_nov_6c14ecd9 | novel | 1 | AFR=0.174, EUR=0.348, AMR=0.250, EAS=0.228 | 0.990 | 0.705 |
| HLA-B_nov_bfa755fa | novel | 1 | AFR=0.221, EUR=0.338, AMR=0.208, EAS=0.234 | 0.993 | 0.705 |
| HLA-B_nov_c0d643d8 | novel | 1 | AFR=0.300, EUR=0.356, AMR=0.189, EAS=0.156 | 0.997 | 0.705 |
| HLA-B_nov_eaa1f00a | novel | 1 | AFR=0.264, EUR=0.333, AMR=0.195, EAS=0.207 | 0.998 | 0.705 |
| HLA-B_nov_a58e54e7 | novel | 1 | AFR=0.347, EUR=0.278, AMR=0.264, EAS=0.111 | 0.998 | 0.711 |
| 15:113 | known | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| 15:75 | known | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| 40:23 | known | 4 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.000 |
| 48:04 | known | 3 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.000 |
| 50:04 | known | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| 52:04 | known | 5 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.000 |
| 52:07 | known | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| 57:61 | known | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_82017cd4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_a6fe07fd | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_99e3a911 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_4cd85998 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_e30cee8e | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_ede46040 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_63aa1ae4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_8976677d | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_f98f02f4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_3f8f404d | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_6a5cc862 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_bfd7e2ea | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_db3f6940 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_15670a99 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_9a1bd411 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d9bb7679 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d857b794 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_4b101839 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_80e14459 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_bcbd6621 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_c0369ebd | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d2dafa0e | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_9a779eb0 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_cbd8313f | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_beb4d124 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_026afb40 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_296c34c7 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_6a8fd8e2 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_a0acfc2b | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_17621e32 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_9f0c5e12 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_6afb2f6a | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d284e465 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_6beaadc6 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_6f36e254 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.000 |
| HLA-B_nov_009a400c | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_ab73e528 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_faaa0ea8 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_44d61bb0 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_f010c0ce | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_027777b6 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_348960dc | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_f18f7081 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_5d23f27a | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_bdf27ef9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_be336e5a | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_021f48ad | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_964805a7 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d02f84ab | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_387c65d5 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_a3f62ec9 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_49b0379f | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_b0c6c4a8 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_663f338e | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_8e0e0943 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_af30735e | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_41eb4588 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_746946c1 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_9e1df2cd | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_5ea37b02 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_f374e762 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_4a14d990 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_091f0220 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_16472b9f | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_2284d26e | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_d1996f44 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_67f61a63 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_85ee38b4 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |
| HLA-B_nov_001d03b1 | novel | 1 | AFR=nan, EUR=nan, AMR=nan, EAS=nan | NA | 0.121 |