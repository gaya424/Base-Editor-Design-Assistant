# core/sequence_parser.py

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from typing import Union, Optional
import sys # Import sys for debugging print

def load_sequence(filepath_or_text: str, file_format: Optional[str] = None) -> SeqRecord:
    """
    Loads a DNA sequence from a file (FASTA, GenBank) or a plain text string.

    Args:
        filepath_or_text (str): The file path to the sequence file or a plain DNA sequence string.
        file_format (str, optional): The format of the file ('fasta' or 'genbank').
                                     If None, the function attempts to guess the format from the file extension.
                                     If the input is a plain text sequence, this should be left as None.

    Returns:
        SeqRecord: A Biopython SeqRecord object containing the sequence and metadata.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file format is not supported or the input is invalid.
        IOError: If there is a problem reading the file.
    """
    print(f"DEBUG: load_sequence received input: '{filepath_or_text}'")
    print(f"DEBUG: Current working directory: {os.getcwd()}")
    
    # Determine if the input string is likely a file path based on its content (extension)
    # This helps distinguish between "my_sequence.fasta" (intended file) and "ATGC" (plain sequence)
    looks_like_file_by_extension = False
    ext = os.path.splitext(filepath_or_text)[1].lower()
    if ext in ('.fasta', '.fna', '.fa', '.gb', '.gbk'):
        looks_like_file_by_extension = True

    # 1. Check if the input is an existing file path
    if os.path.exists(filepath_or_text):
        print(f"DEBUG: '{filepath_or_text}' exists. Attempting to load as a file.")
        try:
            # If no format is provided, try to guess from the extension
            if file_format is None:
                # If it exists and has a common sequence file extension, use that.
                if looks_like_file_by_extension:
                    file_format = ext.lstrip('.') # Remove leading dot
                    # Special handling for .fna, .fa which are fasta
                    if file_format in ('fna', 'fa'):
                        file_format = 'fasta'
                else:
                    # If it's an existing file but doesn't have a known sequence extension,
                    # we can't auto-detect. Raise an error or try a default like 'fasta'.
                    # For now, let's assume it should have a recognizable extension if it's a file.
                    raise ValueError(f"Could not determine file format for existing file '{filepath_or_text}'. Please specify format ('fasta' or 'genbank').")

            # Use Biopython's SeqIO to parse the file
            with open(filepath_or_text, "r") as handle:
                record = next(SeqIO.parse(handle, file_format))
                print(f"DEBUG: Successfully parsed file as '{file_format}'.")
                return record
        except StopIteration:
            raise ValueError(f"Error: No sequence records found in file '{filepath_or_text}' with format '{file_format}'. Is the file empty or malformed?")
        except ValueError as e:
            # Catch Biopython's parsing errors specifically
            raise ValueError(f"Error parsing file '{filepath_or_text}' as '{file_format}'. It may be malformed or the format is incorrect. Details: {e}")
        except Exception as e:
            raise IOError(f"An unexpected error occurred while reading the file: {e}")
            
    # 2. If it's not an *existing* file path, decide if it was *intended* to be a file.
    elif file_format is not None or looks_like_file_by_extension:
        # If a format was specified, or it looked like a file by extension,
        # but os.path.exists returned False, then it's a FileNotFoundError.
        raise FileNotFoundError(f"Error: Sequence file not found at '{filepath_or_text}'.")
        
    # 3. If it's not an existing file AND not explicitly intended as a file,
    # then treat it as a plain sequence string.
    elif filepath_or_text.strip(): # Check if the string is not just whitespace
        print(f"DEBUG: '{filepath_or_text}' does not exist as a file and was not explicitly specified as one. Treating as plain text sequence.")
        # Ensure it only contains valid DNA characters (A, T, C, G, N)
        if not all(base.upper() in 'ATCGN' for base in filepath_or_text):
            raise ValueError("Input string contains characters that are not valid DNA bases (A, T, C, G, N).")
        
        # Create a SeqRecord from the plain sequence
        seq_obj = Seq(filepath_or_text)
        record = SeqRecord(seq_obj, id="plain_sequence", name="plain_sequence")
        return record
        
    # 4. If input is empty or invalid
    else:
        raise ValueError("Input is neither a valid file path nor a non-empty sequence string.")

def get_genomic_slice(sequence: Union[Seq, SeqRecord], start: int, end: int) -> Seq:
    """
    Extracts a slice of a DNA sequence from a Seq or SeqRecord object.

    Args:
        sequence (Union[Seq, SeqRecord]): The input sequence object.
        start (int): The 1-based start position of the slice (inclusive).
        end (int): The 1-based end position of the slice (inclusive).

    Returns:
        Seq: A Biopython Seq object representing the extracted slice.

    Raises:
        ValueError: If the start or end coordinates are invalid (e.g., out of bounds).
    """
    # 1. Get the sequence object from SeqRecord if needed
    if isinstance(sequence, SeqRecord):
        seq_obj = sequence.seq
    elif isinstance(sequence, Seq):
        seq_obj = sequence
    else:
        raise TypeError("Input 'sequence' must be a Biopython Seq or SeqRecord object.")
        
    # 2. Validate the coordinates
    # Convert 1-based coordinates to 0-based Python slicing
    slice_start = start - 1
    slice_end = end
    
    if not (1 <= start <= len(seq_obj) and 1 <= end <= len(seq_obj) and start <= end):
        raise ValueError(
            f"Invalid coordinates: start={start}, end={end}. "
            f"Coordinates must be within the sequence length (1 to {len(seq_obj)}) and start must be <= end."
        )
        
    # 3. Return the slice
    return seq_obj[slice_start:slice_end]

