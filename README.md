# ProteinDataBank-Reading

## Description
This program takes in a pdb file (eg:- name.pdb) and parses it to give the following output in .txt file (eg:- name_output.txt):
- Title/Name of the protein.
- Total length of the protein/ number of residues in the given pdb
- No of chains present in protein and their names(Ascending order if chains are named
numerically followed by alphabetical order)
- All aminoacid ratios present in the protein in alphabetical order.
ex: Ratio(Leu) = no of leucine present/total length of protein.
- Are there any unknown aminoacids present ? If so, mention total count.
- Are there any ligand molecules other than water ? If so, mention their names.
- Calculate all possible phi, psi, omega angles for the given pdb.

## Running the Program
- Make sure the .pdb file is in the same directory as the program
- `python3 script.py`
- Enter File Name - (eg:- 2wsc.pdb)
- A file 2wsc_output.txt will be generated having the result.

## Requirements
- `pip3 install numpy`
