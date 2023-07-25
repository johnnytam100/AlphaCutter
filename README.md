# AlphaCutter

This repository contains the AlphaCutter.py for the removal of non-globular regions from predicted protein structures.

Tam, C., & Iwasaki, W. (2023). AlphaCutter: Efficient removal of non‐globular regions from predicted protein structures. In PROTEOMICS. Wiley. https://doi.org/10.1002/pmic.202300176

![20220816_one-on-one](https://user-images.githubusercontent.com/51283097/212842122-2a05502f-2672-473f-8efe-90fdeb165090.png)

<p align="center">
  <img src="https://github.com/johnnytam100/AlphaCutter/assets/51283097/4178cf66-383a-4ab7-a76b-d3b78e2ec796" width="600">
</p>

---

# How to run AlphaCutter

## Step 1 : Clone, create environment.
````
git clone https://github.com/johnnytam100/AlphaCutter.git

cd AlphaCutter

conda env create -f environment.yml

conda activate AlphaCutter
````

## Step 2 : Copy all PDB files to current directory.
````
cp (path to your PDB files) ./
````

## Step 3 : Run AlphaCutter.

**(1) Default Settings**

````
python AlphaCutter_v111.py \
  --loop_min 20 \
  --helix_min 30 \
  --fragment_min 5 \
  --domain_min 0 \
  --pLDDT_min 0 \
  --local_contact_range 5 \
  --single_out \
  --domain_out
````

**(2) Remove shorter loops and helices**

````
python AlphaCutter_v111.py \
  --loop_min 15 \
  --helix_min 20 \
  --fragment_min 5 \
  --domain_min 0 \
  --pLDDT_min 0 \
  --local_contact_range 5 \
  --single_out \
  --domain_out
````

---

# Options

* `--loop_min`        – Min. loop length to be considered as non-globular region. (Default: `20`. Decrease it to remove shorter loops.) Note: This option is not compulsory if `helix_min` is specified.

* `--helix_min`       – Min. helix length to be considered as non-globular region. (Default: `30`. Decrease it to remove shorter helices.) Note: This option is not compulsory if `loop_min` is specified.

* `--fragment_min`    – Min. domain fragment length to be considered as domain region. (Default: `5`. Increase it to exclude larger domain fragments.)

* `--domain_min`      – Min. domain length to be output as a domain. (Default: `0`. Increase it to exclude larger domains in output.)

* `--pLDDT_min`       – Min. average pLDDT of domain residues to be output as a domain  (Default: `0`. Increase it to exclude domains predicted with low pLDDT.)

* `--local_contact_range`   – Distinguish non-local contact and local contact used in defining non-globular loops and helices. Loops and helices forming any non-local contacts will be redirected as domain fragments. (Default: `5`. Not recommended to change it. If you do, also test with `--loop_min` and `--helix_min`.)

* `--single_out`   – Output all domains as a single PDB file.

* `--domain_out`   – Output every non-contacting domain as a separate PDB file.

---

# AlphaCutter-cleaned SwissProt Protein Structures

Step 1: Download all parts files from https://doi.org/10.5281/zenodo.7944483

Step 2: Run the following
````
cat 20230510_AlphaCutter_cleaned_SwissProt_parts-a? > 20230510_AlphaCutter_cleaned_SwissProt.tar.gz
tar zxvf 20230510_AlphaCutter_cleaned_SwissProt.tar.gz
````
