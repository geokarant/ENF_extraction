# ENF_extraction
Assessing spectral estimation methods for Electric Network Frequency (ENF) Extraction using overplapping windows and Quadratic Interpolation.

A code for ENF extraction is provided including:

1) Signal filtering with a very sharp bandpass filter.

2) ENF extraction procedure using overplapping windows and a Quadratic Interpolation. Very good results achieved due to proper parameter tuning both in the pre-processing step and in the main procedure.

3) Various estimation techniques are employed:
  STFT, Capon, Fast Capon, Fast Iterative Adaptive Approach, Daniell, Blackman-Tukey, Welch, ESPRIT and MUSIC.
  
  
We hope you find our code useful and we'd appreciate it if you cite our work: 

[1] Georgios Karantaidis and Constantine Kotropoulos. 2018. **"ASSESSING SPECTRAL ESTIMATION METHODS FOR ELECTRIC NETWORK FREQUENCY EXTRACTION"**,In 22nd Pan-Hellenic Conference on Informatics (PCI ’18), November 29-December 1, 2018, Athens, Greece. ACM, New York, NY, USA, 6 pages. https://doi.org/10.1145/3291533.3291538

The dataset was created by Spectral Analysis Lab, Dept. of ECE, University of Florida: http://www.sal.ufl.edu/download.html
