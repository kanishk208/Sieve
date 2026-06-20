"""
clinvar_transformer.py — Transform ClinVar JSON to standard format
Handles the mismatch between ClinVar's complex nested structure and 
our simplified {mutation, resnum, known_pathogenicity} format.

Usage:
    from clinvar_transformer import transform_clinvar_json
    variants = transform_clinvar_json("data/VCP/VCP_clinvar_variants.json")
"""

import re
from pathlib import Path


def extract_mutations_from_protein_change(protein_change: str) -> list[dict]:
    """
    Parse protein_change field: "A115S, M158I, D640V" → list of mutations
    
    Returns list of dicts with 'mutation' and 'resnum' keys.
    Example: [{"mutation": "A115S", "resnum": 115}, ...]
    """
    if not protein_change or not isinstance(protein_change, str):
        return []
    
    mutations = []
    # Split by comma and strip whitespace
    parts = [p.strip() for p in protein_change.split(",")]
    
    for part in parts:
        # Pattern: Single letter + digits + single letter (e.g., "A115S", "M158I")
        match = re.match(r"([A-Z])(\d+)([A-Z])", part)
        if match:
            wt_aa, resnum, mut_aa = match.groups()
            mutation = f"{wt_aa}{resnum}{mut_aa}"
            mutations.append({
                "mutation": mutation,
                "resnum": int(resnum)
            })
    
    return mutations


def transform_clinvar_json(json_path: str) -> list[dict]:
    """
    Transform ClinVar JSON format to standard variant format.
    
    Input (ClinVar format):
    {
        "variation_id": "4812857",
        "protein_change": "A115S, A160S",
        "search_class": "pathogenic",
        ...
    }
    
    Output (standard format):
    {
        "mutation": "A115S",
        "resnum": 115,
        "known_pathogenicity": "pathogenic",
        "variation_id": "4812857",
        ...
    }
    """
    import json
    
    try:
        with open(json_path, "r") as f:
            clinvar_records = json.load(f)
    except Exception as e:
        print(f"[!] Failed to load {json_path}: {e}")
        return []
    
    transformed = []
    
    for record in clinvar_records:
        # Extract mutations from protein_change field
        protein_change = record.get("protein_change", "")
        mutations = extract_mutations_from_protein_change(protein_change)
        
        if not mutations:
            continue  # Skip if no valid mutations found
        
        # Map search_class to known_pathogenicity
        search_class = record.get("search_class", "unknown")
        # Normalize: "pathogenic" → "pathogenic", anything else → "benign"
        known_pathogenicity = "pathogenic" if search_class.lower() == "pathogenic" else "benign"
        
        # Create variant record for each mutation
        for mut_info in mutations:
            variant = {
                "mutation": mut_info["mutation"],
                "resnum": mut_info["resnum"],
                "known_pathogenicity": known_pathogenicity,
                "variation_id": record.get("variation_id", ""),
                "title": record.get("title", ""),
                "significance": record.get("significance", ""),
                "protein_change": protein_change,
            }
            transformed.append(variant)
    
    return transformed


def save_transformed_json(variants: list[dict], output_path: str):
    """Save transformed variants to JSON file."""
    import json
    try:
        with open(output_path, "w") as f:
            json.dump(variants, f, indent=2)
        print(f"[✓] Saved transformed JSON: {output_path}")
    except Exception as e:
        print(f"[!] Failed to save {output_path}: {e}")


def main():
    """Transform all ClinVar JSON files in data directories."""
    import json
    from pathlib import Path
    
    data_dir = Path(__file__).parent / "data"
    
    for gene_dir in sorted(data_dir.glob("*")):
        if not gene_dir.is_dir():
            continue
        
        json_file = gene_dir / f"{gene_dir.name}_clinvar_variants.json"
        
        if not json_file.exists():
            print(f"[⚠] Skipping {gene_dir.name} - JSON not found")
            continue
        
        print(f"\n[Transform] Processing {gene_dir.name}...")
        variants = transform_clinvar_json(str(json_file))
        
        if variants:
            # Save transformed version
            output_file = gene_dir / f"{gene_dir.name}_clinvar_variants_transformed.json"
            save_transformed_json(variants, str(output_file))
            print(f"  Transformed {len(variants)} mutations from {len(set([v.get('variation_id') for v in variants]))} variants")
        else:
            print(f"  [!] No valid mutations found")


if __name__ == "__main__":
    main()
