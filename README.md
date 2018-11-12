# Assessing spectral estimation methods for Electric Network Frequency (ENF) Extraction.

Georgios Karantaidis and Constantine Kotropoulos  
PCI 2018  
© 2018 Association for Computing Machinery


# Notes

Code used in [1] for ENF extraction is provided including:

1) Signal filtering with a very sharp bandpass filter (main.m).

2) ENF extraction procedure using overplapping windows and a Quadratic Interpolation (main.m). Very good results achieved due to proper parameter tuning both in the pre-processing step and in the main procedure.

3) Various estimation techniques are employed:
  STFT, Capon, Fast Capon, Fast Iterative Adaptive Approach (treat F-IAA as individual main file), Daniell, Blackman-Tukey, Welch, ESPRIT and MUSIC.
  
  # Important
  
If you use this software you should cite the following in any resulting publication:

[1] Georgios Karantaidis and Constantine Kotropoulos. 2018. **"ASSESSING SPECTRAL ESTIMATION METHODS FOR ELECTRIC NETWORK FREQUENCY EXTRACTION"**,In 22nd Pan-Hellenic Conference on Informatics (PCI ’18), November 29-December 1, 2018, Athens, Greece. ACM, New York, NY, USA, 6 pages. https://doi.org/10.1145/3291533.3291538

The dataset was created by Spectral Analysis Lab, Dept. of ECE, University of Florida and ALL COPYRIGHTS remain with them : http://www.sal.ufl.edu/download.html
