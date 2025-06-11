# Protein-pocket-
Automated pocket detection and HETATM summary from PDB file of 1tup protein using fpocket and Python.

## 🧪 Features

- Run `fpocket` on a given PDB file to detect potential binding pockets.
- Extract and summarize residues and HETATM (ligands, ions, etc.) present in the structure.
- Export results to `.csv` for downstream analysis.
- Visualization of pockets using PyMOL.

---

## 📁 Project Structure

- Step 1 : Parse Structure
- Step 2: Generate a csv for residuue data
- Step 3: Run fpocket
- Step 4: Parse fpocket output
- Step 5: Write to CSV
- Step 6: Open PyMol for visualisation

---

## 🚀 Getting Started

### 📦 Prerequisites

- Python 3.8+
- [`fpocket`](https://github.com/Discngine/fpocket) installed and in your `PATH`
- PyMOL installed and optionally available via command line

### 🔧 Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/protein-pocket-analysis.git
cd protein-pocket-analysis
```
### Dependencies

Biopython : Can be installed with "pip install biopython"

---

## 🧬 Usage 

To run the pocket detection, extraction of data, generating csv files and visualising protein run the following in terminal/command prompt : python protein_pocket.py

---

## 📄 Output

- residue_data.csv – Summary of each residue, including chain ID, residue name, HETATM flag, and sequence number.

- fpocket output directory (e.g., 1tup_out/) with .pml for visualization and pocket structure files. 

---

## 🧠 Motivation

fpocket efficiently detects pockets, but doesn't summarize HETATMs. This script combines both tasks—automating pocket detection and producing a structured CSV of relevant residues, helping streamline druggability and ligand-binding analysis.

---

## 🛠️ Future Improvements

I will be improvinf upon this project by: 

- Extending support for multiple PDBs (beyond 1tup)
- Pocket ranking based on custom metrics

---

### Note 

This code is written in python and run on the macOS. Adjust paths or PyMOL calls if using Windows or Linux.


