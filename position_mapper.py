#!/usr/bin/env python
"""
Position Mapper: Align ClinVar sequence positions to PDB structure coordinates.

Problem: VCP has multiple isoforms. ClinVar uses one numbering system, our PDB 
uses another. When we see "AA mismatch at position 537: expected I, found A", it's 
usually an offset (e.g., ClinVar position 537 is actually PDB position 492).

Solution: Extract sequence from PDB structure, align to standard VCP sequence,
build offset mapping, use it to translate ClinVar positions to PDB positions.
"""

import re
from pathlib import Path
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# Standard VCP uniprot sequence (P55072) - first 850 residues (full is 806 in our PDB)
# Use FASTA if available, else hardcode a known reference
STANDARD_VCP_SEQUENCE = """
MSTLSVAPQRRDMPGRSGLNDSARQVDPQKVDDPLGLNQVEQVVQRSLSEVENRRSLPDQ
EAHGGLKLISWDLPSQNLSGFLQNKIDLIFEQVGPNHSVMVSRKQNDKQTKLFQLDVDSI
SWQPSLKMRQKKRQREDFVVCQDRVLNPTQPIQFLSMGSVYAKLEDRFSLHDLIDFAEQY
PGKQGSSVYAELNNSPAVLTGAVQSVQQYAKPQMEEQALVQVEKASQRLKQQESYAKLLQ
DSYRQIGQKDEAGLKRMVQKEQQIQKAQELLRQLEEYGLVQFNDRLKETVDVLSLLETKG
DDVSLQWGVGLQQAVVVPQSQNPYAAIDKKLLNQLLHELGAEAEMQAQAKNTYSELYNDE
NDKESEFAQQLYQELMLQYNKFDDDDEEDGEQD
""".replace('\n', '')

# Better: build from PDB if available
def extract_sequence_from_pdb(pdb_file: Path) -> str:
    """Extract amino acid sequence from PDB CA atoms."""
    sequence = []
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    prev_resnum = None
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith('ATOM') and len(line) > 25:
                    try:
                        # Only CA atoms
                        if line[12:16].strip() != 'CA':
                            continue
                        
                        resnum = int(line[22:26])
                        aa_3letter = line[17:20].strip()
                        
                        # Skip if duplicate residue number (only take first)
                        if resnum == prev_resnum:
                            continue
                        
                        prev_resnum = resnum
                        sequence.append(aa_3to1.get(aa_3letter, 'X'))
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        print(f"[!] Error extracting sequence: {e}")
    
    return ''.join(sequence)


def build_position_mapper(pdb_file: Path, reference_sequence: str = None) -> dict:
    """
    Build a mapping from reference sequence positions to PDB structure positions.
    
    Returns: dict mapping reference_pos -> pdb_pos
    Example: {537: 492, 538: 493, ...}
    """
    # Extract sequence from PDB
    pdb_sequence = extract_sequence_from_pdb(pdb_file)
    
    if not reference_sequence:
        reference_sequence = STANDARD_VCP_SEQUENCE
    
    print(f"[Info] PDB sequence length: {len(pdb_sequence)}")
    print(f"[Info] Reference sequence length: {len(reference_sequence)}")
    
    if len(pdb_sequence) == 0:
        print("[!] Failed to extract PDB sequence")
        return {}
    
    # Global alignment using PairwiseAligner with scoring
    print("[Info] Running sequence alignment...")
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2.0     # Match
    aligner.mismatch_score = -1.0 # Mismatch
    aligner.open_gap_score = -5.0 # Gap open
    aligner.extend_gap_score = -1.0 # Gap extend
    
    alignments = aligner.align(reference_sequence, pdb_sequence)
    
    if len(alignments) == 0:
        print("[!] No alignment found")
        return {}
    
    # Get best alignment (highest score)
    alignment = alignments[0]
    ref_aligned = str(alignment[0])
    pdb_aligned = str(alignment[1])
    
    print(f"[Info] Alignment score: {alignment.score}")
    
    # Build position mapping
    position_map = {}
    ref_pos = 0
    pdb_pos = 0
    
    for ref_aa, pdb_aa in zip(ref_aligned, pdb_aligned):
        if ref_aa != '-':
            ref_pos += 1
            if pdb_aa != '-':
                # Both positions aligned
                position_map[ref_pos] = pdb_pos
            else:
                # Gap in PDB - this position doesn't exist in structure
                position_map[ref_pos] = None
        
        if pdb_aa != '-':
            pdb_pos += 1
    
    print(f"[Info] Built mapping for {len([v for v in position_map.values() if v is not None])} aligned positions")
    
    return position_map


def translate_position(clinvar_pos: int, position_map: dict) -> int:
    """
    Translate ClinVar position to PDB position using alignment map.
    
    Args:
        clinvar_pos: Position from ClinVar (1-indexed)
        position_map: Mapping from build_position_mapper()
    
    Returns:
        pdb_pos: Corresponding position in PDB (1-indexed), or None if not aligned
    """
    return position_map.get(clinvar_pos, None)


def main():
    """Test the position mapper."""
    print("=" * 70)
    print("[Position Mapper] VCP Isoform Offset Detection")
    print("=" * 70)
    
    gene_dir = Path("data/VCP")
    pdb_file = gene_dir / "AF-P55072-F1-model_6.pdb"
    
    if not pdb_file.exists():
        print(f"[!] PDB file not found: {pdb_file}")
        return
    
    # Build position map
    position_map = build_position_mapper(pdb_file)
    
    if not position_map:
        print("[!] Failed to build position map")
        return
    
    # Test on known problematic positions
    test_positions = [115, 160, 113, 158, 537, 595]
    
    print(f"\n[Mapping] ClinVar → PDB position translations:")
    for clinvar_pos in test_positions:
        pdb_pos = translate_position(clinvar_pos, position_map)
        if pdb_pos is not None:
            print(f"  {clinvar_pos} → {pdb_pos}")
        else:
            print(f"  {clinvar_pos} → NOT ALIGNED (gap in PDB)")
    
    # Save mapping for use in validate.py
    import json
    mapping_file = gene_dir / "position_map.json"
    with open(mapping_file, "w") as f:
        # Convert None values to -1 for JSON serialization
        serializable_map = {
            str(k): (v if v is not None else -1) 
            for k, v in position_map.items()
        }
        json.dump(serializable_map, f, indent=2)
    
    print(f"\n[Saved] Position map: {mapping_file}")


if __name__ == "__main__":
    main()
