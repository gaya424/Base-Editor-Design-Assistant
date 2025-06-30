import sys
import os

# --- Fix for relative import when running this file directly ---
# This block ensures that the project's root directory is added to sys.path
# BEFORE any relative imports are attempted. This allows Python to correctly
# locate modules within the 'core' package when this script is run directly.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Now, import necessary modules. The relative import for editor_specs should now work.
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Any, Union, Optional
import re

# Import editor specifications from the editor_specs module
# This relative import is correct when sgRNA_designer.py is imported as part of the 'core' package.
# With the sys.path fix above, it will also work when run directly.
from core.editor_specs import get_editor_spec # Changed to absolute import path for direct execution

# Define a type alias for clarity for the sgRNA candidate dictionary structure.
# This helps in maintaining consistency across functions that handle sgRNA data.
SgRNACandidate = Dict[str, Any]

def find_potential_sgRNAs(
    sequence: Union[Seq, SeqRecord],
    target_edit_pos: int, # 1-based genomic coordinate of the desired edit
    pam_constraint: str = "NGG",
    sgRNA_length: int = 20
) -> List[SgRNACandidate]:
    """
    Identifies all potential sgRNA candidates (protospacer + PAM) around a target
    edit position on both forward and reverse strands of the provided sequence.

    Args:
        sequence (Union[Seq, SeqRecord]): The genomic DNA sequence to search within.
        target_edit_pos (int): The 1-based genomic coordinate of the desired base edit.
                               This position must contain the base that the editor targets (e.g., 'C' for BE4max).
        pam_constraint (str): The PAM sequence pattern (e.g., "NGG", "NG"). 'N' acts as a wildcard, matching any base.
        sgRNA_length (int): The length of the protospacer (the region of the target DNA that the guide RNA binds to).
                            Typically 20 base pairs.

    Returns:
        List[SgRNACandidate]: A list of dictionaries, each representing a potential sgRNA candidate.
                              Each dictionary contains detailed information about the candidate,
                              including its sequence, genomic coordinates, and placeholders for
                              scores and flags to be filled by subsequent modules.
    Raises:
        TypeError: If the input 'sequence' is not a Biopython Seq or SeqRecord object.
        ValueError: If the 'target_edit_pos' is out of bounds for the given sequence.
    """
    # Extract the raw sequence object from SeqRecord if provided
    if isinstance(sequence, SeqRecord):
        dna_sequence = sequence.seq
    elif isinstance(sequence, Seq):
        dna_sequence = sequence
    else:
        raise TypeError("Input 'sequence' must be a Biopython Seq or SeqRecord object.")

    # Convert the 1-based target_edit_pos to a 0-based index for Python string/list indexing
    target_edit_idx = target_edit_pos - 1
    
    # Validate that the target edit position is within the bounds of the DNA sequence
    if not (0 <= target_edit_idx < len(dna_sequence)):
        raise ValueError(f"Target edit position {target_edit_pos} is out of bounds for sequence length {len(dna_sequence)}.")

    potential_sgRNAs: List[SgRNACandidate] = []
    pam_length = len(pam_constraint)
    
    # Convert the PAM constraint string (e.g., "NGG") into a regular expression pattern
    # '.' in regex matches any character, effectively handling 'N' as a wildcard.
    pam_regex = pam_constraint.replace('N', '.')

    # Iterate over both the forward (+1) and reverse complement (-1) strands
    for strand_val, current_dna_seq in [(1, dna_sequence), (-1, dna_sequence.reverse_complement())]:
        # Iterate through all possible starting positions for the PAM on the current strand.
        # The PAM is typically located immediately 3' to the protospacer.
        # So, if the protospacer is 20bp, and the PAM is 3bp, the region is 20bp guide + 3bp PAM.
        # The loop iterates through all possible PAM start positions.
        for pam_start_on_current_strand in range(len(current_dna_seq) - pam_length + 1):
            # Extract the potential PAM sequence slice
            pam_seq_slice = current_dna_seq[pam_start_on_current_strand : pam_start_on_current_strand + pam_length]

            # Check if the extracted PAM sequence matches the defined PAM constraint pattern
            if re.match(pam_regex, str(pam_seq_slice).upper()):
                # Calculate the start and end positions of the protospacer relative to the current strand.
                # The protospacer is sgRNA_length bases immediately 5' to the PAM.
                protospacer_start_on_current_strand = pam_start_on_current_strand - sgRNA_length
                protospacer_end_on_current_strand = pam_start_on_current_strand - 1 # Inclusive end

                # Ensure the calculated protospacer region is within the bounds of the current strand sequence
                if protospacer_start_on_current_strand >= 0:
                    # Extract the protospacer sequence
                    protospacer_sequence = current_dna_seq[protospacer_start_on_current_strand : protospacer_start_on_current_strand + sgRNA_length]

                    # Map the coordinates back to the original genomic sequence (0-based)
                    if strand_val == 1: # Forward strand
                        genomic_protospacer_start = protospacer_start_on_current_strand
                        genomic_protospacer_end = protospacer_end_on_current_strand
                        genomic_pam_start = pam_start_on_current_strand
                        genomic_pam_end = pam_start_on_current_strand + pam_length - 1
                    else: # Reverse strand
                        # For the reverse complement strand, positions are inverted relative to the original sequence.
                        # A position 'k' on the reverse complement corresponds to 'len(original_seq) - 1 - k' on the original.
                        genomic_protospacer_start = len(dna_sequence) - 1 - protospacer_end_on_current_strand
                        genomic_protospacer_end = len(dna_sequence) - 1 - protospacer_start_on_current_strand
                        genomic_pam_start = len(dna_sequence) - 1 - (pam_start_on_current_strand + pam_length - 1)
                        genomic_pam_end = len(dna_sequence) - 1 - pam_start_on_current_strand

                    # Check if the desired target edit position falls within the genomic coordinates
                    # of this protospacer on the original sequence.
                    if genomic_protospacer_start <= target_edit_idx <= genomic_protospacer_end:
                        # Calculate the relative position of the target base within the 20bp protospacer sequence (0-based)
                        target_base_relative_pos = target_edit_idx - genomic_protospacer_start
                        # Get the original base at the target position from the main DNA sequence
                        original_target_base = dna_sequence[target_edit_idx].upper()

                        # Construct the sgRNA candidate dictionary
                        candidate: SgRNACandidate = {
                            "sequence": str(protospacer_sequence).upper(), # The 20bp protospacer sequence (target DNA strand)
                            "genomic_start": genomic_protospacer_start,   # 0-based start of protospacer on original genome
                            "genomic_end": genomic_protospacer_end,       # 0-based end of protospacer on original genome
                            "pam_sequence": str(pam_seq_slice).upper(),   # The PAM sequence found
                            "pam_start": genomic_pam_start,               # 0-based start of PAM on original genome
                            "pam_end": genomic_pam_end,                   # 0-based end of PAM on original genome
                            "target_base_relative_pos": target_base_relative_pos, # 0-based position of target base within protospacer
                            "strand": strand_val,                         # Strand targeted (+1 for forward, -1 for reverse)
                            "original_target_base": original_target_base, # The base at target_edit_pos BEFORE editing
                            # Placeholders for results from subsequent scoring/analysis modules
                            "is_compatible": False,       # To be updated by check_editable_window
                            "on_target_efficiency": None, # To be filled by scoring.py
                            "off_target_warnings": None,  # To be filled by scoring.py or external tool
                            "bystander_edits": [],        # To be filled by scoring.py
                            "pam_blocking_mutation": False # To be filled by scoring.py
                        }
                        potential_sgRNAs.append(candidate)
    return potential_sgRNAs

def check_editable_window(sgRNA_info: SgRNACandidate, editor_spec: Dict[str, Any]) -> bool:
    """
    Determines if the target base falls within the editable window of the selected base editor.
    Also verifies if the original target base type matches the editor's specified edit type.

    Args:
        sgRNA_info (SgRNACandidate): A dictionary containing sgRNA candidate information,
                                     specifically 'target_base_relative_pos' and 'original_target_base'.
        editor_spec (Dict[str, Any]): A dictionary containing the base editor's specifications,
                                      including 'editable_window_start', 'editable_window_end', and 'edit_type'.

    Returns:
        bool: True if the target base's relative position is within the editor's window AND
              its original base type is compatible with the editor's function (e.g., 'C' for C->T editor),
              False otherwise.
    """
    target_relative_pos = sgRNA_info["target_base_relative_pos"]
    original_target_base = sgRNA_info["original_target_base"]
    
    window_start = editor_spec["editable_window_start"]
    window_end = editor_spec["editable_window_end"]
    editor_edit_type = editor_spec["edit_type"]

    # First, check if the original base at the target position is the type that the editor can modify.
    # For a C->T editor, the target base must be 'C'. For an A->G editor, it must be 'A'.
    if editor_edit_type == "C_TO_T" and original_target_base != 'C':
        return False
    if editor_edit_type == "A_TO_G" and original_target_base != 'A':
        return False

    # Second, check if the target base's relative position within the protospacer
    # falls within the editor's defined editable window.
    return window_start <= target_relative_pos <= window_end

def design_sgRNAs(
    sequence: Union[Seq, SeqRecord],
    desired_edit: str, # e.g., "C123T" or "A256G"
    editor_name: str,
    pam_constraint: str = "NGG",
    sgRNA_length: int = 20
) -> List[SgRNACandidate]:
    """
    The main orchestration function for this module. It designs and filters sgRNAs
    based on a desired genomic edit and a selected base editor.

    Args:
        sequence (Union[Seq, SeqRecord]): The genomic DNA sequence to analyze.
        desired_edit (str): A string specifying the desired base change.
                            Format: <original_base><1-based_position><new_base> (e.g., "C123T", "A256G").
        editor_name (str): The name of the base editor to use (e.g., "BE4max", "ABE8e").
        pam_constraint (str): The PAM sequence pattern (e.g., "NGG").

    Returns:
        List[SgRNACandidate]: A list of compatible sgRNA candidates that meet all criteria:
                              - They contain the target edit position within their protospacer.
                              - The target base is of the correct type for the chosen editor.
                              - The target base falls within the editor's specific editable window.

    Raises:
        ValueError: If the 'desired_edit' format is invalid, the 'editor_name' is unknown,
                    the desired edit type is incompatible with the editor, or the specified
                    original base in 'desired_edit' does not match the actual base in the sequence.
    """
    # 1. Parse the desired_edit string (e.g., "C123T" -> original_base='C', target_pos=123, new_base='T')
    # Use regex to validate the format: Starts with A,C,G,T,N; followed by digits; ends with A,C,G,T,N.
    if not re.match(r"^[ACGTN]\d+[ACGTN]$", desired_edit, re.IGNORECASE):
        raise ValueError(f"Invalid desired_edit format: '{desired_edit}'. Expected format like 'C123T' (e.g., original_base + 1-based_position + new_base).")

    original_base_from_edit_str = desired_edit[0].upper()
    target_pos_1_based = int(desired_edit[1:-1]) # Extract the 1-based position
    new_base_from_edit_str = desired_edit[-1].upper()

    # 2. Retrieve the specifications for the chosen base editor
    editor_spec = get_editor_spec(editor_name)
    if not editor_spec:
        raise ValueError(f"Unknown base editor: '{editor_name}'. Please ensure it is defined in editor_specs.py.")

    # 3. Validate if the desired edit type (e.g., C->T) matches the editor's capability
    expected_editor_edit_type = editor_spec["edit_type"]
    if expected_editor_edit_type == "C_TO_T" and not (original_base_from_edit_str == 'C' and new_base_from_edit_str == 'T'):
        raise ValueError(f"Editor '{editor_name}' is a Cytosine Base Editor (C->T). Desired edit '{desired_edit}' is incompatible with its function.")
    if expected_editor_edit_type == "A_TO_G" and not (original_base_from_edit_str == 'A' and new_base_from_edit_str == 'G'):
        raise ValueError(f"Editor '{editor_name}' is an Adenine Base Editor (A->G). Desired edit '{desired_edit}' is incompatible with its function.")

    # 4. Get the actual base at the target position from the provided genomic sequence
    # This ensures that the user's specified original base actually exists at that position.
    try:
        # Access sequence content, handling both SeqRecord and Seq objects
        if isinstance(sequence, SeqRecord):
            actual_base_at_pos = sequence.seq[target_pos_1_based - 1].upper()
        else:
            actual_base_at_pos = sequence[target_pos_1_based - 1].upper()
    except IndexError:
        raise ValueError(f"Target position {target_pos_1_based} is out of bounds for the provided sequence (length {len(sequence)}).")

    if original_base_from_edit_str != actual_base_at_pos:
        raise ValueError(f"The original base specified in desired_edit ('{original_base_from_edit_str}') does not match "
                         f"the base found at position {target_pos_1_based} in the sequence ('{actual_base_at_pos}').")

    # 5. Find all potential sgRNAs that could target the region around the desired edit
    all_potential_sgRNAs = find_potential_sgRNAs(sequence, target_pos_1_based, pam_constraint, sgRNA_length)

    # 6. Filter the potential sgRNAs based on their compatibility with the editor's editable window
    compatible_sgRNAs: List[SgRNACandidate] = []
    for sgRNA_candidate in all_potential_sgRNAs:
        if check_editable_window(sgRNA_candidate, editor_spec):
            sgRNA_candidate["is_compatible"] = True # Mark as compatible
            compatible_sgRNAs.append(sgRNA_candidate)
            
    return compatible_sgRNAs

if __name__ == '__main__':
    # --- Example Usage for sgRNA_designer.py ---
    # This block demonstrates how the functions in this module can be used.
    # It will only execute when this script is run directly (python sgRNA_designer.py).
    # It is not part of the main application logic when imported as a module.
    print("--- Testing sgRNA_designer.py ---")

    # Create dummy DNA sequences for testing purposes.
    # These sequences are specifically crafted so that:
    # 1. The 'desired_edit' matches the actual base at the target position.
    # 2. There are NGG/NG PAMs in positions that allow the target base to fall
    #    within the respective editor's editable window.

    # Test sequence for C->T edits (e.g., C10T)
    # Target C at 1-based pos 10 (0-indexed 9).
    # BE4max editable window: 4-8 (relative to guide).
    # We need a protospacer where C at genomic 9 falls into relative positions 4-8.
    # If relative pos is 4, protospacer starts at 9-4=5. Protospacer 5-24. PAM at 25.
    # Sequence: ...C...[15 bases]...NGG...
    # Example: ATGCATGCACGCATGCATGCATGCATGGATGCATGCATGCATGC (C at 9, TGG at 25)
    test_seq_for_ct_edit = Seq("ATGCATGCACGCATGCATGCATGCATGGATGCATGCATGCATGC") # C at index 9, TGG (NGG) at index 25

    # Test sequence for A->G edits (e.g., A15G)
    # Target A at 1-based pos 15 (0-indexed 14).
    # ABE8e editable window: 5-9 (relative to guide).
    # We need a protospacer where A at genomic 14 falls into relative positions 5-9.
    # If relative pos is 5, protospacer starts at 14-5=9. Protospacer 9-28. PAM at 29.
    # Sequence: ...A...[14 bases]...NGG...
    # Example: ATGCATGCATGCATGCATGCATGCATGCAGGATGCATGCATGC (A at 14, AGG (NGG) at 29)
    test_seq_for_ag_edit = Seq("ATGCATGCATGCATGCATGCATGCATGCAGGATGCATGCATGC") # A at index 14, AGG (NGG) at index 29

    # --- Test Case 1: BE4max (C->T) on custom sequence with C10T ---
    print("\n--- Test Case 1: BE4max (C->T) on custom sequence with C10T ---")
    desired_edit_1 = "C10T"
    editor_name_1 = "BE4max"

    try:
        compatible_sgRNAs_1 = design_sgRNAs(test_seq_for_ct_edit, desired_edit_1, editor_name_1)
        print(f"  Found {len(compatible_sgRNAs_1)} compatible sgRNAs for {desired_edit_1} with {editor_name_1}:")
        for i, sgRNA in enumerate(compatible_sgRNAs_1):
            print(f"    sgRNA {i+1}:")
            print(f"      Protospacer Seq: {sgRNA['sequence']}")
            print(f"      Genomic Coords (0-idx): {sgRNA['genomic_start']}-{sgRNA['genomic_end']}")
            print(f"      PAM: {sgRNA['pam_sequence']} (Genomic {sgRNA['pam_start']}-{sgRNA['pam_end']})")
            print(f"      Target Base: {sgRNA['original_target_base']} at relative pos {sgRNA['target_base_relative_pos']}")
            print(f"      Strand: {'Forward' if sgRNA['strand'] == 1 else 'Reverse'}")
            print(f"      Compatible: {sgRNA['is_compatible']}")
    except ValueError as e:
        print(f"  Error: {e}")

    # --- Test Case 2: ABE8e (A->G) on custom sequence with A17G ---
    print("\n--- Test Case 2: ABE8e (A->G) on custom sequence with A17G ---")
    desired_edit_2 = "A17G"
    editor_name_2 = "ABE8e"

    try:
        compatible_sgRNAs_2 = design_sgRNAs(test_seq_for_ag_edit, desired_edit_2, editor_name_2)
        print(f"  Found {len(compatible_sgRNAs_2)} compatible sgRNAs for {desired_edit_2} with {editor_name_2}:")
        for i, sgRNA in enumerate(compatible_sgRNAs_2):
            print(f"    sgRNA {i+1}:")
            print(f"      Protospacer Seq: {sgRNA['sequence']}")
            print(f"      Genomic Coords (0-idx): {sgRNA['genomic_start']}-{sgRNA['genomic_end']}")
            print(f"      PAM: {sgRNA['pam_sequence']} (Genomic {sgRNA['pam_start']}-{sgRNA['pam_end']})")
            print(f"      Target Base: {sgRNA['original_target_base']} at relative pos {sgRNA['target_base_relative_pos']}")
            print(f"      Strand: {'Forward' if sgRNA['strand'] == 1 else 'Reverse'}")
            print(f"      Compatible: {sgRNA['is_compatible']}")
    except ValueError as e:
        print(f"  Error: {e}")

    # --- Test Case 3: Incompatible edit type (C10T with ABE8e) ---
    print("\n--- Test Case 3: Incompatible edit type (C10T with ABE8e) ---")
    desired_edit_3 = "C10T" # ABE8e is A->G, not C->T
    editor_name_3 = "ABE8e"
    try:
        design_sgRNAs(test_seq_for_ct_edit, desired_edit_3, editor_name_3)
    except ValueError as e:
        print(f"  Caught expected error: {e}")

    # --- Test Case 4: Original base mismatch (A10G but sequence has C at 10) ---
    # This test now explicitly uses a sequence where the base *does not* match the desired_edit
    print("\n--- Test Case 4: Original base mismatch (A10G but sequence has C at 10) ---")
    # Using test_seq_for_ct_edit which has 'C' at pos 10 (0-indexed 9)
    # We are asking for 'A10G', which should cause a mismatch.
    desired_edit_4 = "A10G"
    editor_name_4 = "ABE8e" # Editor type doesn't matter for this specific mismatch test
    try:
        design_sgRNAs(test_seq_for_ct_edit, desired_edit_4, editor_name_4)
    except ValueError as e:
        print(f"  Caught expected error: {e}")