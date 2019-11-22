# PUMA

Overview
Untargeted metabolomics comprehensively characterizes small molecules and elu-cidates activities of biochemical pathways within a biological sample. Despite computational ad-vances, interpreting collected measurements and determining their biological role remains a chal-lenge.

To interpret measurements, we present an inference-based approach, termed Probabil-istic modeling for Untargeted Metabolomics Analysis (PUMA). Our approach captures measure-ments and known information about the sample under study in a generative model and uses sto-chastic sampling to compute posterior probability distributions. PUMA predicts the likelihood of pathways being active, and then derives a probabilistic annotation, which assigns chemical iden-tities to the measurements. PUMA is validated on synthetic datasets. When applied to test cases, the resulting pathway activities are biologically meaningful and distinctly different from those obtained using statistical pathway enrichment techniques. Annotation results are in agreement to those obtained using other tools that utilize additional information in the form of spectral signa-tures. Importantly, PUMA annotates many additional measurements. 

Getting Started

There are two examples: 

data/CHO_cell: Chinese Hampster Ovary Cell example provided in the paper

data/human_urine: A human urine example, see reference in the paper

Start running quickly by executing run_puma.py

