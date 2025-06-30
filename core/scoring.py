# core/scoring.py

import sys
import os
from typing import List, Dict, Any, Union, Optional

# --- Fix for relative import when running this file directly ---
# This block ensures that the project's root directory is added to sys.path
# BEFORE any relative imports are attempted. This allows Python to correctly
# locate modules within the 'core' package when this script is run directly.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Import editor specifications from the editor_specs module
from core.editor_specs import get_editor_spec

# Define a type alias for clarity for the sgRNA candidate dictionary structure.
SgRNACandidate = Dict[str, Any]

def calculate_on_target_efficiency(sgRNA_sequence: str, editor_name: str) -> float:
    """
    Calculates a placeholder on-target efficiency score for an sgRNA.
    Future improvements could integrate models like Doench '16, Azimuth, or deep learning models.

    Args:
        sgRNA_sequence (str): The 20bp protospacer sequence.
        editor_name (str): The name of the base editor (e.g., "BE4max", "ABE8e").
                           Can be used to apply editor-specific scoring rules in advanced versions.

    Returns:
        float: A dummy on-target efficiency score (e.g., based on GC content or a fixed value).
               Returns 0.0 if the sequence is empty, otherwise a value between 0.0 and 1.0.
    """
    # Placeholder implementation: A very simple score based on GC content.
    # This is NOT a biologically accurate prediction but serves as a functional placeholder.
    if not sgRNA_sequence:
        return 0.0

    gc_count = sgRNA_sequence.count('G') + sgRNA_sequence.count('C')
    # Normalize GC content to a 0-1 scale, then apply a simple transformation.

    gc_content = gc_count / len(sgRNA_sequence)
    
    # A simple, arbitrary scoring logic for demonstration:
    # Higher GC content might (hypothetically) lead to higher efficiency in some contexts.
    # This is a very rough approximation.
    score = 0.5 + (gc_content - 0.5) * 0.8 # Scale around 0.5, max 0.9, min 0.1
    return round(max(0.1, min(0.9, score)), 2) # Ensure score is between 0.1 and 0.9

def identify_bystander_edits(
    sgRNA_info: SgRNACandidate,
    editor_spec: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    Identifies potential bystander edits within the editable window of an sgRNA.
    Bystander edits occur when other bases of the target type (e.g., 'C' for BE4max, 'A' for ABE8e)
    are present in the editable window besides the primary target base.

    Args:
        sgRNA_info (SgRNACandidate): A dictionary containing sgRNA candidate information,
                                     including 'sequence', 'target_base_relative_pos',
                                     'genomic_start', and 'strand'.
        editor_spec (Dict[str, Any]): A dictionary containing the base editor's specifications,
                                      including 'editable_window_start', 'editable_window_end', and 'edit_type'.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, where each dictionary describes a potential
                              bystander edit. Each dictionary contains:
                              - 'base': The base that could be edited (e.g., 'C' or 'A').
                              - 'relative_pos': The 0-indexed position within the protospacer.
                              - 'genomic_pos': The 0-indexed genomic position on the original sequence.
                              Returns an empty list if no bystander edits are found.
    """
    protospacer_sequence = sgRNA_info["sequence"]
    primary_target_relative_pos = sgRNA_info["target_base_relative_pos"]
    genomic_protospacer_start = sgRNA_info["genomic_start"]
    strand = sgRNA_info["strand"]

    window_start = editor_spec["editable_window_start"]
    window_end = editor_spec["editable_window_end"]
    editor_edit_type = editor_spec["edit_type"]

    target_base_type = ''
    if editor_edit_type == "C_TO_T":
        target_base_type = 'C'
    elif editor_edit_type == "A_TO_G":
        target_base_type = 'A'
    else:
        # Handle unknown edit types if necessary, though current specs only include C_TO_T and A_TO_G
        return []

    bystander_edits_found: List[Dict[str, Any]] = []

    # Iterate through the editable window of the protospacer sequence
    for relative_pos in range(window_start, window_end + 1):
        # Ensure the relative position is within the bounds of the protospacer sequence
        if 0 <= relative_pos < len(protospacer_sequence):
            # Check if this position is not the primary target base
            if relative_pos != primary_target_relative_pos:
                base_at_pos = protospacer_sequence[relative_pos].upper()
                
                # If the base matches the editor's target type, it's a potential bystander edit
                if base_at_pos == target_base_type:
                    # Calculate the genomic position of the bystander edit
                    genomic_pos = genomic_protospacer_start + relative_pos
                    
                    bystander_edits_found.append({
                        "base": base_at_pos,
                        "relative_pos": relative_pos,
                        "genomic_pos": genomic_pos,
                        "strand": strand # Include strand for clarity
                    })
    return bystander_edits_found

def check_pam_blocking_mutation(
    sgRNA_info: SgRNACandidate,
    editor_spec: Dict[str, Any]
) -> bool:
    """
    Checks if the desired target edit, if it occurs, would fall within or modify
    the PAM sequence, potentially blocking Cas enzyme recognition.

    Args:
        sgRNA_info (SgRNACandidate): A dictionary containing sgRNA candidate information,
                                     including 'target_base_relative_pos', 'genomic_start',
                                     'pam_sequence', 'pam_start', and 'pam_end'.
        editor_spec (Dict[str, Any]): A dictionary containing the base editor's specifications.
                                      This is primarily used to confirm the PAM sequence.

    Returns:
        bool: True if the target edit position overlaps with the PAM sequence, False otherwise.
    """
    target_base_relative_pos = sgRNA_info["target_base_relative_pos"]
    protospacer_sequence_length = len(sgRNA_info["sequence"]) # Should be 20 for typical sgRNAs
    
    # The PAM is located immediately 3' to the protospacer.
    # If the protospacer is 0-indexed from 0 to 19, the PAM starts at index 20 (relative to the guide).
    # This assumes the guide is 20bp.
    pam_relative_start = protospacer_sequence_length
    pam_relative_end = protospacer_sequence_length + len(sgRNA_info["pam_sequence"]) - 1

    # Check if the target base's relative position falls within the PAM's relative positions
    # If the target base is being edited, and that position is part of the PAM, it's a blocking mutation.
    if pam_relative_start <= target_base_relative_pos <= pam_relative_end:
        return True
    
    return False

if __name__ == '__main__':
    # --- Example Usage for scoring.py ---
    print("--- Testing scoring.py ---")

    # Dummy editor spec for testing (should match editor_specs.py)
    test_editor_specs = {
        "be4max": {
            "description": "BE4max (Cytosine Base Editor)",
            "editable_window_start": 4,
            "editable_window_end": 8,
            "pam_sequence": "NGG",
            "edit_type": "C_TO_T"
        },
        "abe8e": {
            "description": "ABE8e (Adenine Base Editor)",
            "editable_window_start": 5,
            "editable_window_end": 9,
            "pam_sequence": "NGG",
            "edit_type": "A_TO_G"
        },
        "spcas9-ng": {
            "description": "SpCas9-NG (Broader PAM CBE - Hypothetical)",
            "editable_window_start": 4,
            "editable_window_end": 8,
            "pam_sequence": "NG",
            "edit_type": "C_TO_T"
        },
    }

    # Helper function to get editor spec for tests
    def get_test_editor_spec(name):
        return test_editor_specs.get(name.lower())

    # --- Test Case 1: calculate_on_target_efficiency ---
    print("\n--- Test Case 1: calculate_on_target_efficiency ---")
    seq_high_gc = "GCGCGCGCGCGCGCGCGCGC" # 20bp, 100% GC
    seq_low_gc = "ATATATATATATATATATAT"   # 20bp, 0% GC
    seq_mixed_gc = "ATGCATGCATGCATGCATGC" # 20bp, 50% GC

    score_high = calculate_on_target_efficiency(seq_high_gc, "BE4max")
    score_low = calculate_on_target_efficiency(seq_low_gc, "ABE8e")
    score_mixed = calculate_on_target_efficiency(seq_mixed_gc, "BE4max")

    print(f"  Score for '{seq_high_gc}': {score_high}")
    print(f"  Score for '{seq_low_gc}': {score_low}")
    print(f"  Score for '{seq_mixed_gc}': {score_mixed}")
    print("  (Note: Scores are placeholders and based on simple GC content for demonstration)")

    # --- Test Case 2: identify_bystander_edits (BE4max, C->T) ---
    print("\n--- Test Case 2: identify_bystander_edits (BE4max, C->T) ---")
    # Protospacer: ATGCATGCACGTGCATGCATGC (C at 0-idx 9 is target)
    # Editable window (4-8): C at 4, C at 6, G at 7, C at 8
    # Target is C at 9. So, C at 4, C at 6, C at 8 are potential bystanders.
    sgRNA_info_be4max: SgRNACandidate = {
        "sequence": "ATGCATGCACGTGCATGCATGC", # Length 22, for example, assuming longer seq for illustration
        "genomic_start": 100,
        "genomic_end": 119,
        "pam_sequence": "NGG",
        "pam_start": 120,
        "pam_end": 122,
        "target_base_relative_pos": 9, # C at 0-indexed 9
        "strand": 1,
        "original_target_base": "C"
    }
    # Adjust sequence to be 20bp for consistency with typical sgRNA_length
    sgRNA_info_be4max["sequence"] = "ATGCATGCACGTGCATGCAT" # 20bp
    sgRNA_info_be4max["target_base_relative_pos"] = 9 # C at 0-indexed 9

    editor_spec_be4max = get_test_editor_spec("BE4max")
    bystanders_be4max = identify_bystander_edits(sgRNA_info_be4max, editor_spec_be4max)

    print(f"  Protospacer: {sgRNA_info_be4max['sequence']}")
    print(f"  Target C at relative pos: {sgRNA_info_be4max['target_base_relative_pos']}")
    print(f"  BE4max window: {editor_spec_be4max['editable_window_start']}-{editor_spec_be4max['editable_window_end']}")
    print(f"  Bystander edits found ({len(bystanders_be4max)}):")
    if bystanders_be4max:
        for b in bystanders_be4max:
            print(f"    - Base: {b['base']}, Relative Pos: {b['relative_pos']}, Genomic Pos: {b['genomic_pos']}")
    else:
        print("    None.")

    # --- Test Case 3: identify_bystander_edits (ABE8e, A->G) ---
    print("\n--- Test Case 3: identify_bystander_edits (ABE8e, A->G) ---")
    # Protospacer: GGGGGAGGGGAGGGGAGGGG (A at 0-idx 5 is target)
    # Editable window (5-9): A at 5, G at 6, G at 7, G at 8, G at 9
    # Target is A at 5. No other A's in window (5-9).
    sgRNA_info_abe8e: SgRNACandidate = {
        "sequence": "GGGGGAGGGGAGGGGAGGGG", # 20bp
        "genomic_start": 200,
        "genomic_end": 219,
        "pam_sequence": "NGG",
        "pam_start": 220,
        "pam_end": 222,
        "target_base_relative_pos": 5, # A at 0-indexed 5
        "strand": 1,
        "original_target_base": "A"
    }
    editor_spec_abe8e = get_test_editor_spec("ABE8e")
    bystanders_abe8e = identify_bystander_edits(sgRNA_info_abe8e, editor_spec_abe8e)

    print(f"  Protospacer: {sgRNA_info_abe8e['sequence']}")
    print(f"  Target A at relative pos: {sgRNA_info_abe8e['target_base_relative_pos']}")
    print(f"  ABE8e window: {editor_spec_abe8e['editable_window_start']}-{editor_spec_abe8e['editable_window_end']}")
    print(f"  Bystander edits found ({len(bystanders_abe8e)}):")
    if bystanders_abe8e:
        for b in bystanders_abe8e:
            print(f"    - Base: {b['base']}, Relative Pos: {b['relative_pos']}, Genomic Pos: {b['genomic_pos']}")
    else:
        print("    None.")

    # --- Test Case 4: check_pam_blocking_mutation (No blocking) ---
    print("\n--- Test Case 4: check_pam_blocking_mutation (No blocking) ---")
    # sgRNA length 20. PAM starts at relative pos 20.
    # Target at relative pos 9 (well before PAM).
    sgRNA_info_no_block: SgRNACandidate = {
        "sequence": "ATGCATGCACGTGCATGCAT", # 20bp
        "genomic_start": 100,
        "genomic_end": 119,
        "pam_sequence": "NGG",
        "pam_start": 120,
        "pam_end": 122,
        "target_base_relative_pos": 9, # C at 0-indexed 9
        "strand": 1,
        "original_target_base": "C"
    }
    editor_spec_be4max_for_pam = get_test_editor_spec("BE4max")
    is_blocking_no = check_pam_blocking_mutation(sgRNA_info_no_block, editor_spec_be4max_for_pam)
    print(f"  Target at relative pos {sgRNA_info_no_block['target_base_relative_pos']}. PAM starts at relative pos {len(sgRNA_info_no_block['sequence'])}.")
    print(f"  Is PAM blocking mutation? {is_blocking_no} (Expected: False)")

    # --- Test Case 5: check_pam_blocking_mutation (Blocking) ---
    print("\n--- Test Case 5: check_pam_blocking_mutation (Blocking) ---")
    # Target at relative pos 20 (which is the start of PAM for 20bp guide).
    sgRNA_info_block: SgRNACandidate = {
        "sequence": "ATGCATGCATGCATGCATGC", # 20bp
        "genomic_start": 300,
        "genomic_end": 319,
        "pam_sequence": "NGG",
        "pam_start": 320,
        "pam_end": 322,
        "target_base_relative_pos": 20, # Hypothetically, if target falls on 1st base of PAM
        "strand": 1,
        "original_target_base": "G" # Example, could be any base in PAM
    }
    editor_spec_be4max_for_pam_block = get_test_editor_spec("BE4max")
    is_blocking_yes = check_pam_blocking_mutation(sgRNA_info_block, editor_spec_be4max_for_pam_block)
    print(f"  Target at relative pos {sgRNA_info_block['target_base_relative_pos']}. PAM starts at relative pos {len(sgRNA_info_block['sequence'])}.")
    print(f"  Is PAM blocking mutation? {is_blocking_yes} (Expected: True)")
