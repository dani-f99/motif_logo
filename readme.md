--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman
                                      2025
--------------------------------------------------------------------------------


# PROJECT: Clonal Selection Pressure Analysis

## 1. OVERVIEW
This program aim is to visualize the amino acid usage bais across 
defiened BCR heavy chain variable section. If metadata exists the program
can output multiple sub-plots defined by the metadata in order to comapre 
motifs across different conditions.


## 2. PREREQUISITES
Please ensure the following python modules are installed:

- Pandas 
- NumPy 
- Matplotlib 
- logomaker


## 3. USAGE GUIDE
1. Place the sequences dataset in the `input` folder. 
2. Follow the steps illustrated in the 'motif_logo.ipynb' notebool. 
3. Access the output figures via the 'output' folder.


## 4. DIRECTORY STRUCTURE
The program uses the following folder structure: 

- input/  : Input folder (for the sequences data). 
- output/ : Result - output motif figures. 

		   
## 5. RESOUCES
- Logomaker documentation: https://logomaker.readthedocs.io/en/latest/
