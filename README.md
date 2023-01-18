# AlphaCutter

This repository contains the AlphaCutter.py for the removal of disordered regions from predicted protein structures.

![20220816_one-on-one](https://user-images.githubusercontent.com/51283097/212842122-2a05502f-2672-473f-8efe-90fdeb165090.png)

(link to paper)

(major results)

---

# How to run AlphaCutter

## Step 1 : Copy all PDB files to current directory.
````
cp (path to your PDB files) ./
````

## Step 2 : Run AlphaCutter.

**(1) Default Settings**

````
python AlphaCutter_v111.py 
  --loop_min 20 
  --helix_min 30 
  --fragment_min 0 
  --domain_min 0 
  --pLDDT_min 0 
  --local_contact_range 5 
  --single_out 
  --domain_out
````

**(2) Remove shorter loops and helices**

````
python AlphaCutter_v111.py 
  --loop_min 15 
  --helix_min 20 
  --fragment_min 0 
  --domain_min 0 
  --pLDDT_min 0 
  --local_contact_range 5 
  --single_out 
  --domain_out
````

---

# Options

* `--loop_min`        – Min. loop length to be consider as disordered region. (Default: `20`. Decrease it to remove shorter loops.)

* `--helix_min`       – Min. helix length to be consider as disordered region. (Default: `30`. Decrease it to remove shorter helices.)

* `--fragment_min`    – Min. domain fragment length to be consider as domain region. (Default: `0`. Increase it to exclude small domain fragments.)

* `--domain_min`      – Min. domain length to be output as a domain. (Default: `0`. Increase it to exclude small domains in output.)

* `--pLDDT_min`       – Min. average pLDDT of domain residues to be output as a domain  (Default: `0`. Increase it to exclude domains predicted with low pLDDT.)

* `--local_contact_range`   – Distinguish non-local contact and local contact used in defining disordered loops and helices. Loops and helices forming any non-local contacts will be redirected as domain fragments. (Default: `5`. Not recommended to change it. If you do, also test with `--loop_min` and `--helix_min`.)

* `--single_out`   – Output every domain as a separate PDB file.

* `--domain_out`   – Output all domains as a single PDB file.
