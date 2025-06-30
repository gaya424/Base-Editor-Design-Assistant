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
    
    # Check if the input is a file path
    if os.path.exists(filepath_or_text):
        print(f"DEBUG: '{filepath_or_text}' exists. Attempting to load as a file.")
        try:
            # If no format is provided, try to guess from the extension
            if file_format is None:
                ext = os.path.splitext(filepath_or_text)[1].lower()
                if ext in ('.fasta', '.fna', '.fa'):
                    file_format = 'fasta'
                elif ext in ('.gb', '.gbk'):
                    file_format = 'genbank'
                else:
                    # If it's a file but extension is unknown, try fasta as a common fallback
                    print(f"DEBUG: Unknown file extension '{ext}'. Trying 'fasta' as default format.")
                    file_format = 'fasta' # Try fasta as a default if extension is ambiguous but file exists

            # Use Biopython's SeqIO to parse the file
            with open(filepath_or_text, "r") as handle:
                # SeqIO.parse returns an iterator, we take the first record
                # This is suitable for single-record files.
                # If a file has multiple records, you'd need to decide how to handle them.
                record = next(SeqIO.parse(handle, file_format))
                print(f"DEBUG: Successfully parsed file as '{file_format}'.")
                return record
        except FileNotFoundError as e: # This should ideally be caught by os.path.exists, but good to keep
            raise FileNotFoundError(f"Error: Sequence file not found at '{filepath_or_text}'.") from e
        except StopIteration:
            raise ValueError(f"Error: No sequence records found in file '{filepath_or_text}' with format '{file_format}'. Is the file empty or malformed?")
        except ValueError as e:
            raise ValueError(f"Error parsing file '{filepath_or_text}' as '{file_format}'. It may be malformed or the format is incorrect. Details: {e}")
        except Exception as e:
            raise IOError(f"An unexpected error occurred while reading the file: {e}")
            
    # If it's not a file path, treat it as a plain sequence string
    elif filepath_or_text.strip(): # Check if the string is not just whitespace
        print(f"DEBUG: '{filepath_or_text}' does not exist as a file. Treating as plain text sequence.")
        # Ensure it only contains valid DNA characters (A, T, C, G, N)
        if not all(base.upper() in 'ATCGN' for base in filepath_or_text):
            raise ValueError("Input string contains characters that are not valid DNA bases (A, T, C, G, N).")
        
        # Create a SeqRecord from the plain sequence
        seq_obj = Seq(filepath_or_text)
        record = SeqRecord(seq_obj, id="plain_sequence", name="plain_sequence")
        return record
        
    # If input is empty or invalid
    else:
        raise ValueError("Input is neither a valid file path nor a non-empty sequence string.")

def get_genomic_slice(sequence: Union[Seq, SeqRecord], start: int, end: int) -> Seq:
    """
    Extracts a slice of a DNA sequence from a Seq or SeqRecord object.
    ... (rest of docstring) ...
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

if __name__ == '__main__':
    # --- Example Usage ---
    print("--- Testing sequence_parser.py ---")

    # Example 1: Load from a plain text sequence
    plain_seq_text = "ATGCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    print(f"\n1. Loading from plain text: '{plain_seq_text[:15]}...'")
    try:
        seq_record_plain = load_sequence(plain_seq_text)
        print(f"   Success! ID: {seq_record_plain.id}, Length: {len(seq_record_plain.seq)}")
        print(f"   Sequence: {seq_record_plain.seq}")
    except Exception as e:
        print(f"   Failed: {e}")

    # Example 2: Create a dummy FASTA file for testing
    dummy_fasta_content = ">test_gene_A\nATGCATGCATGCATGC\n>test_gene_B\nCGATCGATCGATCGAT"
    dummy_fasta_path = "dummy_test.fasta"
    with open(dummy_fasta_path, "w") as f:
        f.write(dummy_fasta_content)
    
    print(f"\n2. Loading from a FASTA file: '{dummy_fasta_path}'")
    try:
        # Note: SeqIO.parse will only give the first record for a multi-record file
        seq_record_fasta = load_sequence(dummy_fasta_path)
        print(f"   Success! ID: {seq_record_fasta.id}, Description: {seq_record_fasta.description}")
        print(f"   Sequence: {seq_record_fasta.seq}")
    except Exception as e:
        print(f"   Failed: {e}")
        
    # Example 3: Test with a non-existent file
    print("\n3. Testing with a non-existent file: 'non_existent.fasta'")
    try:
        load_sequence("non_existent.fasta")
    except FileNotFoundError as e:
        print(f"   Caught expected error: {e}")
    except Exception as e:
        print(f"   Failed with unexpected error: {e}")
        
    # Example 4: Test getting a slice from the loaded sequence
    print("\n4. Getting a slice of the sequence from the FASTA file (positions 5 to 10)")
    try:
        # The sequence is ATGCATGCATGCATGC
        # 1-based positions 5 to 10 are: 'ATGCAT'
        seq_slice = get_genomic_slice(seq_record_fasta, 5, 10)
        print(f"   Success! Extracted slice: {seq_slice}")
    except Exception as e:
        print(f"   Failed: {e}")
        
    # Example 5: Test invalid coordinates for slicing
    print("\n5. Testing slice with invalid coordinates (start > end)")
    try:
        get_genomic_slice(seq_record_fasta, 10, 5)
    except ValueError as e:
        print(f"   Caught expected error: {e}")

    # Clean up the dummy file
    if os.path.exists(dummy_fasta_path):
        os.remove(dummy_fasta_path)
        print(f"\nCleaned up {dummy_fasta_path}")
