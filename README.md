# Evaluating Mechanistic Whole-Brain Models of Empirical MEG Data

This project focuses on reimplementing, extending, and critically evaluating mechanistic whole-brain models for empirical MEG data analysis, with a particular emphasis on Deco's multiple-frequency brain model and its integration with the Wilson-Cowan framework.

## Project Overview

We recreate and analyze Deco's multiple-frequency brain model, extend it using the Wilson-Cowan framework, and provide critical insights into model assumptions and parameterization. Our work aims to improve the accuracy and reproducibility of computational models in simulating complex brain dynamics.

## Key Features

1. **Model Reimplementation and Analysis**
   - Recreation of Deco's multiple-frequency brain model
   - Application of bifurcation analysis, nullcline plotting, and Euler-Maruyama method
   - Investigation of replicability challenges and parameter sensitivity

2. **Enhanced Model Integration and Comparison**
   - Extension of Stuart-Landau model with Wilson-Cowan framework
   - Implementation of signal processing techniques (band-pass filtering, envelope functional connectivity analysis)
   - Comparative analysis of model performance against empirical MEG data

3. **Critical Evaluation and Insights**
   - Identification of oversights in noise parameterization and model assumptions
   - Recommendations for refining computational models in neuroscience

## Results

Our analysis reveals:

- Challenges in replicability and parameter sensitivity of Deco's model
- Superior alignment of the Wilson-Cowan model with empirical MEG data
- Critical insights into noise parameterization and model assumptions
For a detailed presentation of our findings, please refer to our [full report](final_report.pdf).

## References

- [AAL90](https://nilearn.github.io/dev/modules/description/aal_SPM12.html): result of an automated anatomical parcellation of the spatially normalized single-subject high-resolution T1 volume provided by the Montreal Neurological Institute (MNI) (Collins et al.).


