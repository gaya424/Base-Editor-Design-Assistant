from typing import Dict, Any, Optional

# Define the specifications for various base editors.
# Each editor is a dictionary containing its specific properties:
# - 'editable_window_start': The 0-indexed start position of the editable window
#                            relative to the 5' end of the 20bp guide sequence.
# - 'editable_window_end': The 0-indexed end position (inclusive) of the editable window.
# - 'pam_sequence': The Protospacer Adjacent Motif (PAM) sequence recognized by the
#                   Cas enzyme associated with this editor (e.g., 'NGG').
# - 'edit_type': The type of base conversion this editor performs ('C_TO_T' or 'A_TO_G').
BASE_EDITOR_SPECS: Dict[str, Dict[str, Any]] = {
    "be4max": {
        "description": "BE4max (Cytosine Base Editor)",
        "editable_window_start": 4,  # e.g., C at position 4 (0-indexed) of the guide
        "editable_window_end": 8,    # e.g., C at position 8 (0-indexed) of the guide (inclusive)
        "pam_sequence": "NGG",
        "edit_type": "C_TO_T"
    },
    "abe8e": {
        "description": "ABE8e (Adenine Base Editor)",
        "editable_window_start": 5,  # e.g., A at position 5 (0-indexed) of the guide
        "editable_window_end": 9,    # e.g., A at position 9 (0-indexed) of the guide (inclusive)
        "pam_sequence": "NGG",
        "edit_type": "A_TO_G"
    }
}

def get_editor_spec(editor_name: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves the specifications for a given base editor.

    Args:
        editor_name (str): The name of the base editor (case-insensitive, e.g., "BE4max", "abe8e").

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing the editor's specifications if found,
                                  otherwise None.
    """
    return BASE_EDITOR_SPECS.get(editor_name.lower())

if __name__ == '__main__':
    # --- Example Usage ---
    print("--- Testing editor_specs.py ---")

    # Test 1: Retrieve specifications for BE4max
    editor_name_1 = "BE4max"
    print(f"\n1. Retrieving specs for '{editor_name_1}':")
    be4max_spec = get_editor_spec(editor_name_1)
    if be4max_spec:
        print(f"   Description: {be4max_spec.get('description')}")
        print(f"   Editable Window: {be4max_spec.get('editable_window_start')} to {be4max_spec.get('editable_window_end')} (0-indexed)")
        print(f"   PAM Sequence: {be4max_spec.get('pam_sequence')}")
        print(f"   Edit Type: {be4max_spec.get('edit_type')}")
    else:
        print(f"   Editor '{editor_name_1}' not found.")

    # Test 2: Retrieve specifications for ABE8e (case-insensitive check)
    editor_name_2 = "abe8e"
    print(f"\n2. Retrieving specs for '{editor_name_2}':")
    abe8e_spec = get_editor_spec(editor_name_2)
    if abe8e_spec:
        print(f"   Description: {abe8e_spec.get('description')}")
        print(f"   Editable Window: {abe8e_spec.get('editable_window_start')} to {abe8e_spec.get('editable_window_end')} (0-indexed)")
        print(f"   PAM Sequence: {abe8e_spec.get('pam_sequence')}")
        print(f"   Edit Type: {abe8e_spec.get('edit_type')}")
    else:
        print(f"   Editor '{editor_name_2}' not found.")

    # Test 3: Retrieve specifications for a non-existent editor
    editor_name_3 = "NonExistentEditor"
    print(f"\n3. Retrieving specs for '{editor_name_3}':")
    non_existent_spec = get_editor_spec(editor_name_3)
    if non_existent_spec:
        print(f"   Found specs: {non_existent_spec}")
    else:
        print(f"   Editor '{editor_name_3}' not found (as expected).")

    # Test 4: Accessing a specific property directly (after checking existence)
    print("\n4. Directly accessing a property:")
    be4max_spec_direct = get_editor_spec("BE4max")
    if be4max_spec_direct:
        print(f"   BE4max editable window start: {be4max_spec_direct['editable_window_start']}")
    else:
        print("   BE4max spec not found for direct access test.")