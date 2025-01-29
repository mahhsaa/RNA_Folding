# RNA Folding Objective Function

This project implements an RNA folding energy estimator using interatomic distance distributions from PDB files. It includes scripts to train a scoring model, visualize interaction profiles, and evaluate RNA structures.

## 📋 Project Description
The goal is to create an objective function that estimates the Gibbs free energy of RNA conformations by:
1. **Training** on known 3D structures (PDB files).
2. **Plotting** interaction profiles (text-based).
3. **Scoring** predicted RNA structures.

## 🛠 Installation
1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/rna-folding.git
   cd rna-folding


## Requirements:
Python 3.6+ (no external packages required).

## Usage
1. **Training the Model**
Process PDB files to generate scoring profiles: python train.py data/pdb_files/ scores/
Input: PDB files in data/pdb_files/.
Output: 10 score files (e.g., AA.txt, AU.txt) in scores/.

2. **Plotting Interaction Profiles**
Generate terminal-friendly profiles: python plot.py scores/

3. **Scoring an RNA Structure**
Calculate the Gibbs free energy of a PDB structure: python score.py examples/target.pdb scores/
