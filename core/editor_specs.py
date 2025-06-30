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

