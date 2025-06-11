import os
import pandas as pd
from Bio.PDB import PDBParser

# Step 1: Parse full protein structure (optional, for other analyses)
def extract_residues(pdb_file):
    parser = PDBParser(QUIET=True)
    residue_data = []
    structure = parser.get_structure("protein", pdb_file)
    for model in structure:
        for chain in model:
            print(f"Chain {chain.id}")
            for residue in chain:
                resname = residue.get_resname()
                resid = residue.get_id()[1]
                residue_data.append((chain.id, resname, resid))
    return residue_data

# Step 2: Run P2Rank prediction
def run_p2rank(pdb_file):
    print("Running P2Rank...")
    os.system(f"./p2rank_2.5/prank predict -f {pdb_file}")

# Step 3: Read P2Rank's pockets.csv output
def parse_p2rank_csv(pdb_file):
    base_name = os.path.splitext(pdb_file)[0]
    csv_path = os.path.join(f"{base_name}_out", "pockets.csv")
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"{csv_path} not found. Ensure P2Rank ran correctly.")
    
    df = pd.read_csv(csv_path)
    return df

# Step 4: Optionally write a summarized CSV
def summarize_pockets_csv(df, output_file="pockets_summary.csv"):
    summary = []
    for _, row in df.iterrows():
        pocket_id = row['pocket_number']
        residue_list = row['pocket_residues']  # Format: "A:SER:45;A:THR:46"
        residue_count = len(residue_list.split(";")) if pd.notna(residue_list) else 0
        summary.append([pocket_id, residue_count, residue_list])
    
    summary_df = pd.DataFrame(summary, columns=["Pocket_ID", "Residue_Count", "Residues"])
    summary_df.to_csv(output_file, index=False)
    print(f"CSV summary saved to {output_file}")

# ---------- Run the pipeline ----------
if __name__ == "__main__":
    pdb_file = "1tup.pdb"

    # (Optional) Step 1: Print available residues in the structure
    residues = extract_residues(pdb_file)
    print("Extracted residues (first 5):", residues[:5], "...")

    # Step 2: Run P2Rank
    run_p2rank(pdb_file)

    # Step 3: Read P2Rank's CSV output
    pockets_df = parse_p2rank_csv(pdb_file)
    print(f"\nDetected {len(pockets_df)} pockets.")
    print(pockets_df[["pocket_number", "score", "pocket_residues"]].head())

    # Step 4: Create summarized CSV
    summarize_pockets_csv(pockets_df)
