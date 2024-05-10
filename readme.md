#  RNAfcg: RNA flexibility prediction based on topological centrality and global features



RNAfcg is a method for predicting RNA flexibility.


## Preparation
Before running the program, the solvent-accessible surface area of each nucleotide residue in RNA needs to be calculated using NACCESS software and the secondary structure information of RNA needs to be extracted using the server (http://web.x3dna.org/).In addition, some python packages need to be installed to ensure that the program runs properly.

RNAfcg uses the following dependencies:
- python 3.8
- numpy
- pandas
- pickle
- networkx 2.8.2
- math


## Example
An RNA PDB file, '4LNT.cif', is used as an example to show the process.


## How to run
**Step 1:** Put the '4LNT.cif' file and RNA_fcg.py in the same directory with other codes.

**Step 2:** Put the file '4LNT_RA.rsa' calculated by the NACCESS software and the secondary structure file 'SS_4LNT_RA.csv' extracted by X3DNA's server into the same folder. Run script 'feature_extraction.py' to extract the features and store the file as '4LNT_RA_all_feature.csv'.

**Step 3:** Run script 'RNA_fcg.py' to predict the B-factor value of the corresponding RNA and store the output file as 'predict_b_factor.csv' file.



The B factor values predicted by RNA_fcg are stored in the predict_b_factor.csv file.

---

## Help
For any questions, please contact us by chunhuali@bjut.edu.cn or chunhuali@bjut.edu.cn