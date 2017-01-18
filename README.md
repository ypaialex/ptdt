#pTDT (Polygenic Transmission Disequilibrium Test) `v1.0.0`

`pTDT` is a python-based command line tool for conducting significance tests on the over-transmission of polygenic risk from parents to offspring.

This `pTDT` wiki describes: 

1. The file formats accepted by `pTDT` 
2. The basic usage of `pTDT` 
3. Additional analytic flags available through `pTDT`  
4. A brief `pTDT` tutorial

For a detailed description of the `pTDT`, please consult [Weiner et al. 2017] (http://biorxiv.org/content/early/2016/11/23/089342). If you publish using `pTDT`, please cite that paper.

## System requirements and `pTDT` download
1. `Python 3.x`
2. `argparse`
4. `numpy`
5. `pandas`
6. `scipy`

These may be accessed XXX

`pTDT` may be downloaded using XXX

## 1) Basic file formats

`pTDT` requires two files to run:

  1. Polygenic Risk Score file (from here: `PRS file`)
    * Contents: Mapping between individuals and their polygenic risk scores
    * Format: Text file with four columns: 1) Family ID 2) Individual ID 3) Case/control status 4) polygenic risk score. 
    * Number of rows in `PRS file` = number of individuals in cohort
  2. Family Structure file (from here: `Structure file`)
    * Contents: File identifying which individual belongs in which family
    * Format: Text file with four or five columns: 1) Family ID 2) Proband Individual ID 3) Father Individual ID 4) Mother Individual ID   5) Sibling IID (*optional*)
    * Number of rows in `Structure file` = number of families in cohort ~ (length of `PRS file` / n) (n = 3 for trio families, n = 4 for quads)

Individuals missing from the `Structure file` should be marked as "NA"

**Important note about unique IDs**

Individuals must be assigned unique IDs. It is this unique ID that allows `ptdt` to map individuals from their "Individual ID" in the `PRS file` to their appropriate "Individual ID" in the `Structure file`.

## 2) Basic usage 

All `pTDT` commands must begin by calling 1) Python and 2) `pTDT` as follows:

```
python ptdt.py
```



If successfully run, the following output will print:

```
usage: ptdt.py [-h] [--prs FILENAME X Y [FILENAME X Y ...]] --structure
                     FILENAME [--subset FILENAME] [--quad] [--print]
                     --out OUT
ptdt.py: error: the following arguments are required: --structure, --out
```

A basic `pTDT` analysis must include 1) `PRS file` 2) `Structure file` 3) out name, as follows:
```
python ptdt.py --prs [PRS file] --structure [Structure file] --out [outname]

```
If successfully run, the following output will print:

```
Writing log file to [outname].ptdt.log
Options invoked:
	--PRS [PRS file]
	--structure [Structure file]

n families loaded from structure file.
Creating pTDT matrix... done.
QC pass.
m probands used in pTDT analysis (n-m skipped due to missingness).

------------------------------------------
Proband analysis
pTDT mean: X SD
pTDT SE: Y
pTDT pvalue: Z
------------------------------------------

--output: Results written to [outname].ptdt.result.
```
Brief description of basic output
* `QC pass` If correlation between average parent PRS and offspring PRS > 0.2 (flags data scramble)
* `pTDT mean` Average of the pTDT deviation distribution
* `pTDT SE` Standard error of the pTDT deviation distribution
* `pTDT pvalue` pvalue of one-sample t-test of the pTDT deviation distribution against null of pTDT mean = 0
* `missingness` If any family members in the `Struructure file` are not found in the `PRS file`

## 3) Additional flags

The following are flags that can be added to the `pTDT` basic usage. All available flags are also available with the `-help` flag as follows:

```
python ptdt.py --help
```

`--subset [subset file]` 
* Contents: Allows `pTDT` to be performed on a subset of the families in the `Structure file`
* Format: `subset file` is a text file containing one column of Family IDs (no header)

`--quad`
* Contents: If Sibling column present in `Structure file`, invoking `--quad` performs pTDT `pTDT` analysis on these siblings as well
* Format: No file required, flag sufficient

`--print`
* Contents: Outputs a table to working directory that contains intermediate values in the `pTDT` calculation. The table resembles a `Structure file` where the individual IDs have been replaced by their corresponding PRS from the `PRS file` and the pTDT values calculated in the far right columns
* Format: Text file with 6 columns: 1) Family ID 2) proband PRS 3) father PRS 4) mother PRS 5) average parent PRS 6) proband pTDT value. If `--quad` invoked, additional sibling PRS and sibling PRS columns added.  

`--prs` X Y 
* Contents: The default `PRS file` format is set to match the output from [PRS scoring in Plink] (http://pngu.mgh.harvard.edu/~purcell/plink/profile.shtml), with the Individual ID in the 2nd column and the PRS in the 4th column. This modification to the `--prs` flag allows `pTDT` to accept files with different column ordering
* Format: `--prs X Y` where X is the integer number of the column containing the Individual ID and Y is the integer number of the column containing the PRS value. Default is X = 2 and Y = 4.

## 4) Tutorial

You can download the three files for the tutorial [here] (some extension). Those files are:
* `PRS file` called *demo_prs_file* which contains individual-PRS mapping for 7,780 individuals (1,945 families X 4 members per family)
* `Structure file` called *demo_structure_file* which contains mapping between family ID and individual ID for 1,945 families
* `Subset file` called *demo_subset_file* which contains a subset of 1,410 families from the *demo_structure_file*

**Basic pTDT on proband and unaffected sibling**

The following flags would be used:
```
python ptdt.py --prs demo_prs_upload --structure demo_structure_upload --quad --out demo_prs_quad
```
The following output will print:

```
pTDT Script v1.0.3 (10 Jan 2017)
Created by Alex Pai, Daniel Weiner and Elise Robinson

Writing log file to demo_prs_quad.ptdt.log
Options invoked:
	--PRS demo_prs_file
	--structure demo_structure_file
	--quad

1945 families loaded from structure file.
Creating pTDT matrix... done.
QC pass.
1945 probands used in pTDT analysis (0 skipped due to missingness).
1945 siblings used in pTDT analysis (0 skipped due to missingness).

------------------------------------------
Proband analysis
pTDT mean: 3.841E-02 SD
pTDT SE: 1.690E-02
pTDT pvalue: 2.310E-02

Sibling analysis
pTDT mean: 8.762E-03 SD
pTDT SE: 1.674E-02
pTDT pvalue: 6.007E-01
------------------------------------------

--output: Results written to demo_prs_quad.ptdt.result.
```
The commands invoked and the results are written to `demo_prs_quad.ptdt.log`and `demo_prs_quad.ptdt.result` respectively.

**pTDT on subset of cohort, without unaffected siblings**

The following flags would be used:
```
python ptdt.py --prs demo_prs_upload --structure demo_structure_upload --subset demo_subset_file --out demo_prs_subset

```
The following output will print:

```
pTDT Script v1.0.3 (10 Jan 2017)
Created by Alex Pai, Daniel Weiner and Elise Robinson

Writing log file to demo_prs_subset.ptdt.log
Options invoked:
	--PRS demo_prs_file
	--structure demo_structure_file
	--subset demo_subset_file

1945 families loaded from structure file.
1410 families loaded from subset file (535 families excluded from analysis).
Creating pTDT matrix... done.
QC pass.
1410 probands used in pTDT analysis (0 skipped due to missingness).

------------------------------------------
Proband analysis
pTDT mean: 3.122E-02 SD
pTDT SE: 2.034E-02
pTDT pvalue: 1.251E-01
------------------------------------------

--output: Results written to demo_prs_subset.ptdt.result.
```
The commands invoked and the results are written to `demo_prs_subset.ptdt.log`and `demo_prs_subset.ptdt.result` respectively.

## Updating pTDT

Download of newest version of `pTDT` available through `git` XXXXXX

## Authors

Alex Pai, Daniel Weiner and Elise Robinson (Stanley Center for Psychiatric Research at Broad Institute)
