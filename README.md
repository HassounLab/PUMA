# PUMA

**Overview** <br>
Untargeted metabolomics comprehensively characterizes small molecules and elucidates activities of biochemical pathways within a biological sample. Despite computational advances, interpreting collected measurements and determining their biological role remains a chal-lenge.

To interpret measurements, we present an inference based approach, termed Probabilistic modeling for Untargeted Metabolomics Analysis (PUMA). Our approach captures measurements and known information about the sample under study in a generative model and uses stochastic sampling to compute posterior probability distributions. PUMA predicts the likelihood of pathways being active, and then derives a probabilistic annotation, which assigns chemical identities to the measurements. PUMA is validated on synthetic datasets. When applied to test cases, the resulting pathway activities are biologically meaningful and distinctly different from those obtained using statistical pathway enrichment techniques. Annotation results are in agreement to those obtained using other tools that utilize additional information in the form of spectral signatures. Importantly, PUMA annotates many additional measurements. 

# Getting Started
`Python 3.5.6` is used for development. We recommend installing packages using Anaconda as follows:<br>
`conda create --name PUMA --file enviroment.yml`<br>
`conda activate PUMA`<br>

# Dataset
Two examples are provided: <br>
[CHO_cell:](data/CHO_cell) Chinese Hampster Ovary Cell example provided in the paper <br>
[human_urine:](data/human_urine) A human urine example, see reference in the paper <br>

# How to run the code?
 Start executing the code by `python run_puma.py`

# Authors
This software is written by Ramtin Hosseini, Neda Hassanpour, Liping Liu, and Soha Hassoun (Soha.Hassoun@tufts.edu).<br>
**Paper title:** Pathway Activity Analysis and Metabolite Annotation for Untargeted Metabolomics using Probabilistic Modeling
