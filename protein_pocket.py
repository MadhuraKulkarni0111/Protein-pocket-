import os
import csv
import pandas as pd
from Bio.PDB import PDBParser 
import subprocess

# ---------- Step 1: Parse Structure ----------
def extract_residues(pdb_file):
    parser = PDBParser(QUIET=True)
    residue_data = [] # a list to store all the information 
    structure = parser.get_structure("protein", pdb_file) # parsers through the file and treats it like a gene named p53
    for model in structure: # the for loop iterates through the model eg. X-ray, NMR structure
        for chain in model: # the model comprises of chains (distinct polypeptide chain), this loop iterates through thr model
            print(f"Chain {chain.id}") # prints the chain id
            for residue in chain: 
                resname = residue.get_resname() # returns a 3-letter amino acid code  
                het_flag, resid, _ = residue.get_id() # residue ID tuple (HETATM flag, residue sequence number, insertion code)
                residue_data.append((resname, het_flag, resid))
    return residue_data

# ---------- Step 2: Generate a csv for residuue data ----------
def write_residue_data(residue_data, output_file_residues="residue_data.csv"):
    with open(output_file_residues, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Amino acid", "HETATM Flag", "Residue Sequence Number"])
        for resname, het_flag, resid in residue_data:
            writer.writerow([resname, het_flag, resid])
    print(f"CSV summary saved to {output_file_residues}")


# ---------- Step 3: Run fpocket ----------
def run_fpocket(pdb_file):
    print("Running fpocket...")
    os.system(f"./fpocket/bin/fpocket -f {pdb_file}")

# ---------- Step 4: Parse fpocket output ----------
def parse_fpocket_output(pdb_file):

    base_name = os.path.splitext(pdb_file)[0]
    pocket_dir = os.path.join(base_name + "_out", "pockets" )
    pockets_info = []

   
    for pocket_id in range(1, 100):  
            #pocket_file = f"{pocket_dir}/pocket{pocket_id}_atm.pdb"
            pocket_file = os.path.join(pocket_dir, f"pocket{pocket_id}_atm.pdb")
            if not os.path.exists(pocket_file):
                break  # No more pockets
            with open(pocket_file) as f:
                residues = set()
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        resname = line[17:20].strip()
                        resid = line[22:26].strip()
                        chain = line[21].strip()
                        residues.add((chain, resname, resid))
                pockets_info.append({
                    "Pocket_ID": pocket_id,
                    "Residues": list(residues),
                    "Residue_Count": len(residues)
                })
    return pockets_info

# ---------- Step 5: Write to CSV ----------
def write_pockets_to_csv(pockets_info, output_file="pockets_summary.csv"):
    with open(output_file, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Pocket_ID", "Residue_Count", "Residues"])
        for pocket in pockets_info:
            residue_str = "; ".join([f"{c}:{r}{i}" for (c, r, i) in pocket["Residues"]])
            writer.writerow([pocket["Pocket_ID"], pocket["Residue_Count"], residue_str])
    print(f"CSV summary saved to {output_file}")

# ---------- Step 6: Open PyMol for visualisation ----------
def visualize_with_pymol(pocket_dir):
    pymol_script = os.path.join(pocket_dir, "1tup.pml")  
    if os.path.exists(pymol_script):
        print(f"Opening PyMOL session: {pymol_script}")
        subprocess.run(["open", "-a", "PyMOL", pymol_script])
    else:
        print(f"pymol script not found at {pymol_script}")

# ---------- Run Full Pipeline ----------
if __name__ == "__main__":
    pdb_file = "1tup.pdb"  

    # Step 1: Parse and extract residue info
    residues = extract_residues(pdb_file)
    print("Extracted residues:", residues[:5], "...")

    # Step 2
    write_residue_data(residues)

    # Step 3
    run_fpocket(pdb_file)

    # Step 4
    pocket_info = parse_fpocket_output(pdb_file)

    # Step 5
    write_pockets_to_csv(pocket_info)

    # Step 6
    output_dir = "1tup_out"
    visualize_with_pymol("1tup_out") 