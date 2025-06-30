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
    are present in the editable window besides the primary target_base_relative_pos.
    The primary target itself is NOT considered a bystander edit.

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
            # A bystander edit is an editable base *other than* the primary target in the window.
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

