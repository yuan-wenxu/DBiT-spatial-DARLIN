# DBiT-spatial-DARLIN Technical Documentation

This document provides a detailed technical description of the DBiT-spatial-DARLIN pipeline.

- Step scripts:
  - `dbit_mrna.sh` 
  - `image.sh`
  - `dbit_amplicon.sh`
  - `plot_cell_filtered.sh`

---

## 0. Pipeline Overview

This project is used for quality control of sequencing data from DBiT.

There are four shell scripts to process data from three modalities.

`plot_cell_filtered.sh` is used to integrate data from three modalities.

---

## 1. Transcriptome

1. Extract spatial barcode and UMI based on sequence information.
  
    This is the structure of barcode and UMI.

    ![Barcode and UMI structure](image/barcode.png)


2. Using STARsolo to align transcriptome to genome.
3. Cluster using Leiden algorithm and plot according to spatial location.

## 2. Image

1. Extract the sampling area based on the transcriptome plot results.
2. Cut the image according to the actual sampling area and use stardist to predict the number of cells.

## 3. Amplicon

1. Extract the DARLIN sequence. (optional)
2. Extract spatial barcode and UMI based on sequence information.
3. Correct DARLIN sequence.
