"""
grantham.py — Amino Acid Chemical Distance (Grantham 1974)

Provides a fast, zero-dependency lookup for the chemical distance between
any pair of standard amino acids.  Used by the validation pipeline as a
"chemical safety valve" to override false-positive topology signals for
conservative substitutions.

Reference:
    Grantham R. "Amino acid difference formula to help explain protein
    evolution." Science 185, 862-864 (1974).
"""
from __future__ import annotations

# ─── Full Grantham Distance Matrix ──────────────────────────────────────────
# Upper-triangular lookup keyed by sorted 2-tuple of 1-letter amino acid codes.
# Values are integer Grantham distances.  Identical pairs → 0.

_GRANTHAM_RAW: dict[tuple[str, str], int] = {
    ("A", "R"): 112, ("A", "N"): 111, ("A", "D"): 126, ("A", "C"):  195,
    ("A", "Q"): 91,  ("A", "E"): 107, ("A", "G"):  60, ("A", "H"):  86,
    ("A", "I"):  94, ("A", "L"):  96, ("A", "K"): 106, ("A", "M"):  84,
    ("A", "F"): 113, ("A", "P"):  27, ("A", "S"):  99, ("A", "T"):  58,
    ("A", "W"): 148, ("A", "Y"): 112, ("A", "V"):  64,

    ("R", "N"):  86, ("R", "D"):  96, ("R", "C"): 180, ("R", "Q"):  43,
    ("R", "E"):  54, ("R", "G"): 125, ("R", "H"):  29, ("R", "I"):  97,
    ("R", "L"): 102, ("R", "K"):  26, ("R", "M"):  91, ("R", "F"): 97,
    ("R", "P"): 103, ("R", "S"): 110, ("R", "T"):  71, ("R", "W"): 101,
    ("R", "Y"):  77, ("R", "V"):  96,

    ("N", "D"):  23, ("N", "C"): 139, ("N", "Q"):  46, ("N", "E"):  42,
    ("N", "G"):  80, ("N", "H"):  68, ("N", "I"): 149, ("N", "L"): 153,
    ("N", "K"):  94, ("N", "M"): 142, ("N", "F"): 158, ("N", "P"):  91,
    ("N", "S"):  46, ("N", "T"):  65, ("N", "W"): 174, ("N", "Y"): 143,
    ("N", "V"): 133,

    ("D", "C"): 154, ("D", "Q"):  61, ("D", "E"):  45, ("D", "G"):  94,
    ("D", "H"):  81, ("D", "I"): 168, ("D", "L"): 172, ("D", "K"): 101,
    ("D", "M"): 160, ("D", "F"): 177, ("D", "P"): 108, ("D", "S"):  65,
    ("D", "T"):  85, ("D", "W"): 181, ("D", "Y"): 160, ("D", "V"): 152,

    ("C", "Q"): 154, ("C", "E"): 170, ("C", "G"): 159, ("C", "H"): 174,
    ("C", "I"):  198, ("C", "L"): 198, ("C", "K"): 202, ("C", "M"): 196,
    ("C", "F"): 205, ("C", "P"): 169, ("C", "S"): 112, ("C", "T"): 149,
    ("C", "W"): 215, ("C", "Y"): 194, ("C", "V"): 192,

    ("Q", "E"):  29, ("Q", "G"):  87, ("Q", "H"):  24, ("Q", "I"): 109,
    ("Q", "L"): 113, ("Q", "K"):  53, ("Q", "M"):  101, ("Q", "F"): 116,
    ("Q", "P"):  76, ("Q", "S"):  68, ("Q", "T"):  42, ("Q", "W"): 130,
    ("Q", "Y"):  99, ("Q", "V"):  96,

    ("E", "G"):  98, ("E", "H"):  40, ("E", "I"): 134, ("E", "L"): 138,
    ("E", "K"):  56, ("E", "M"): 126, ("E", "F"): 140, ("E", "P"):  93,
    ("E", "S"):  80, ("E", "T"):  65, ("E", "W"): 152, ("E", "Y"): 122,
    ("E", "V"): 121,

    ("G", "H"):  98, ("G", "I"): 135, ("G", "L"): 138, ("G", "K"): 127,
    ("G", "M"): 127, ("G", "F"): 153, ("G", "P"):  42, ("G", "S"):  56,
    ("G", "T"):  59, ("G", "W"): 184, ("G", "Y"): 147, ("G", "V"): 109,

    ("H", "I"):  94, ("H", "L"):  99, ("H", "K"):  32, ("H", "M"):  87,
    ("H", "F"):  100, ("H", "P"):  77, ("H", "S"):  89, ("H", "T"):  47,
    ("H", "W"): 115, ("H", "Y"):  83, ("H", "V"):  84,

    ("I", "L"):   5, ("I", "K"): 102, ("I", "M"):  10, ("I", "F"):  21,
    ("I", "P"):  95, ("I", "S"): 142, ("I", "T"):  89, ("I", "W"): 61,
    ("I", "Y"):  33, ("I", "V"):  29,

    ("L", "K"): 107, ("L", "M"):  15, ("L", "F"):  22, ("L", "P"):  98,
    ("L", "S"): 145, ("L", "T"):  92, ("L", "W"):  61, ("L", "Y"):  36,
    ("L", "V"):  32,

    ("K", "M"):  95, ("K", "F"): 102, ("K", "P"): 103, ("K", "S"): 121,
    ("K", "T"):  78, ("K", "W"): 110, ("K", "Y"):  85, ("K", "V"):  97,

    ("M", "F"):  28, ("M", "P"):  87, ("M", "S"): 135, ("M", "T"):  81,
    ("M", "W"):  67, ("M", "Y"):  36, ("M", "V"):  21,

    ("F", "P"):  114, ("F", "S"): 155, ("F", "T"): 103, ("F", "W"):  40,
    ("F", "Y"):  22, ("F", "V"):  50,

    ("P", "S"):  74, ("P", "T"):  38, ("P", "W"): 147, ("P", "Y"): 110,
    ("P", "V"):  68,

    ("S", "T"):  58, ("S", "W"): 177, ("S", "Y"): 144, ("S", "V"): 124,

    ("T", "W"): 128, ("T", "Y"):  92, ("T", "V"):  69,

    ("W", "Y"):  37, ("W", "V"):  88,

    ("Y", "V"):  55,
}

# Build a symmetric O(1) lookup dict
GRANTHAM: dict[tuple[str, str], int] = {}
for (a, b), d in _GRANTHAM_RAW.items():
    GRANTHAM[(a, b)] = d
    GRANTHAM[(b, a)] = d

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def grantham_distance(aa1: str, aa2: str) -> int:
    """
    Return the Grantham distance between two amino acids (1-letter codes).
    Returns 0 for identical residues, -1 for non-standard amino acids.
    """
    a, b = aa1.upper(), aa2.upper()
    if a == b:
        return 0
    if a not in STANDARD_AA or b not in STANDARD_AA:
        return -1
    return GRANTHAM.get((a, b), -1)


def is_chemically_conservative(wt: str, mutant: str, threshold: int = 50) -> bool:
    """
    Return True if the amino acid substitution is chemically conservative
    (Grantham distance < threshold).

    Default threshold = 50 (standard conservative boundary per Grantham 1974).
    Examples of conservative swaps (distance < 50):
        I → V  (29)   I → L  (5)    I → M  (10)    F → Y  (22)
        V → M  (21)   L → M  (15)   N → D  (23)    K → R  (26)
    """
    d = grantham_distance(wt, mutant)
    if d < 0:
        return False  # unknown → not conservative
    return d < threshold


# ── Group-based conservative classification (for quick display) ──────────
CONSERVATIVE_GROUPS = [
    frozenset({"I", "L", "V", "M"}),        # Aliphatic / Hydrophobic
    frozenset({"F", "Y", "W"}),              # Aromatic
    frozenset({"K", "R", "H"}),              # Positively Charged
    frozenset({"D", "E"}),                   # Negatively Charged
    frozenset({"S", "T", "N", "Q"}),         # Polar Uncharged
    frozenset({"A", "G"}),                   # Small
    frozenset({"C"}),                        # Special — Cysteine
    frozenset({"P"}),                        # Special — Proline
]


def chemical_group(aa: str) -> str:
    """Return the chemical group name for an amino acid."""
    aa = aa.upper()
    groups = {
        "I": "Aliphatic", "L": "Aliphatic", "V": "Aliphatic", "M": "Aliphatic",
        "F": "Aromatic", "Y": "Aromatic", "W": "Aromatic",
        "K": "Positive", "R": "Positive", "H": "Positive",
        "D": "Negative", "E": "Negative",
        "S": "Polar", "T": "Polar", "N": "Polar", "Q": "Polar",
        "A": "Small", "G": "Small",
        "C": "Cysteine", "P": "Proline",
    }
    return groups.get(aa, "Unknown")


if __name__ == "__main__":
    # Quick sanity checks
    test_pairs = [
        ("V", "M", 21), ("I", "V", 29), ("F", "I", 21),
        ("R", "C", 180), ("D", "E", 45), ("I", "L", 5),
    ]
    print("Grantham Distance Sanity Checks:")
    for a, b, expected in test_pairs:
        d = grantham_distance(a, b)
        ok = "✓" if d == expected else f"✗ (got {d})"
        print(f"  {a} → {b}: {d:>4}  conservative={is_chemically_conservative(a, b)}  {ok}")
