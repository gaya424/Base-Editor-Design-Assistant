import argparse
import os
import pandas as pd
import sys

# Add the project root to the system path to allow absolute imports for core modules
# This is crucial when running base_editor.py directly from the project root.
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(current_dir) # base_editor.py is in the root
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Import core modules
try:
    from core import sequence_parser
    from core import editor_specs
    from core import sgRNA_designer
    from core import scoring
    from core import visualization
except ImportError as e:
    print(f"Error importing core modules: {e}")
    print("Please ensure your project structure is correct and all core modules are present.")
    print("Expected structure: Base-Editor-Design-Assistant/core/*.py")
    sys.exit(1)

def main():
    """
    Main function to parse arguments, orchestrate the sgRNA design workflow,
    and generate output.
    """
    parser = argparse.ArgumentParser(
        description="A Python tool to help researchers design efficient and precise CRISPR base-editing experiments."
    )

    # Define command-line arguments
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the genomic DNA sequence file (FASTA, GenBank) or a plain text sequence."
    )
    parser.add_argument(
        "--edit",
        type=str,
        required=True,
        help="Desired base edit (e.g., C123T for C to T at 1-based position 123, A256G)."
    )
    parser.add_argument(
        "--editor",
        type=str,
        required=True,
        help="Base editor type (e.g., BE4max, ABE8e, SpCas9-NG). Must be defined in editor_specs.py."
    )
    parser.add_argument(
        "--pam",
        type=str,
        default="NGG",
        help="PAM constraint (e.g., NGG, NG). Default is NGG."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results",
        help="Directory to save output files (CSV, HTML plots). Default is 'results'."
    )

    args = parser.parse_args()

    # --- 1. Validate and Create Output Directory ---
    output_dir = args.output
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output will be saved to: {os.path.abspath(output_dir)}")
    except OSError as e:
        print(f"Error: Could not create output directory '{output_dir}'. Details: {e}")
        sys.exit(1)

    # --- 2. Load Genomic Sequence ---
    print(f"\nLoading genomic sequence from: {args.input}")
    try:
        genomic_sequence_record = sequence_parser.load_sequence(args.input)
        print(f"Successfully loaded sequence (ID: {genomic_sequence_record.id}, Length: {len(genomic_sequence_record.seq)} bp).")
    except (FileNotFoundError, ValueError, IOError, TypeError) as e:
        print(f"Error loading sequence: {e}")
        sys.exit(1)

    # --- 3. Get Editor Specifications ---
    print(f"Retrieving specifications for editor: {args.editor}")
    editor_spec = editor_specs.get_editor_spec(args.editor)
    if not editor_spec:
        print(f"Error: Editor '{args.editor}' not found in editor_specs.py. Please check the name.")
        sys.exit(1)
    print(f"Using editor: {editor_spec.get('description', args.editor)}")

    # --- 4. Design sgRNAs ---
    print(f"Designing sgRNAs for edit '{args.edit}' with PAM '{args.pam}'...")
    try:
        compatible_sgRNAs = sgRNA_designer.design_sgRNAs(
            sequence=genomic_sequence_record,
            desired_edit=args.edit,
            editor_name=args.editor,
            pam_constraint=args.pam
        )
        print(f"Found {len(compatible_sgRNAs)} compatible sgRNA candidates.")
    except ValueError as e:
        print(f"Error during sgRNA design: {e}")
        sys.exit(1)

    if not compatible_sgRNAs:
        print("No compatible sgRNAs found for the specified edit and editor. Exiting.")
        sys.exit(0) # Exit gracefully if no sgRNAs are found

    # --- 5. Score and Flag sgRNAs ---
    print("Calculating on-target efficiency, identifying bystander edits, and checking for PAM-blocking mutations...")
    for sgRNA in compatible_sgRNAs:
        # Calculate on-target efficiency (placeholder)
        sgRNA["on_target_efficiency"] = scoring.calculate_on_target_efficiency(
            sgRNA["sequence"], args.editor
        )
        
        # Identify bystander edits
        sgRNA["bystander_edits"] = scoring.identify_bystander_edits(
            sgRNA, editor_spec
        )
        
        # Check for PAM-blocking mutation
        sgRNA["pam_blocking_mutation"] = scoring.check_pam_blocking_mutation(
            sgRNA, editor_spec
        )
        
        # Placeholder for off-target warnings (CRISPRitz/Cas-OFFinder integration)
        # For now, it's None, but could be populated here if an external tool was called.
        sgRNA["off_target_warnings"] = None # Initialize as None, update if off-target logic is added

    # --- 6. Store Results in a pandas DataFrame and Sort/Rank ---
    print("Processing results and ranking sgRNAs...")
    results_df = pd.DataFrame(compatible_sgRNAs)

    # Sort by on-target efficiency (descending)
    results_df = results_df.sort_values(by="on_target_efficiency", ascending=False).reset_index(drop=True)
    results_df["rank"] = results_df.index + 1 # Add a 1-based rank column

    # Reorder columns for better readability in output
    desired_columns = [
        "rank", "sequence", "on_target_efficiency", "is_compatible",
        "original_target_base", "target_base_relative_pos",
        "bystander_edits", "pam_blocking_mutation", "off_target_warnings",
        "pam_sequence", "genomic_start", "genomic_end", "pam_start", "pam_end", "strand"
    ]
    # Ensure all desired columns exist, add missing ones if any
    existing_columns = results_df.columns.tolist()
    final_columns = [col for col in desired_columns if col in existing_columns]
    results_df = results_df[final_columns]

    # --- 7. Generate Output Files ---
    print("Generating output files...")

    # Save results to CSV
    csv_output_path = os.path.join(output_dir, "sgRNA_candidates.csv")
    results_df.to_csv(csv_output_path, index=False)
    print(f"  Saved sgRNA candidates to: {csv_output_path}")

    # Generate individual sgRNA visualizations for top candidates
    num_plots = min(5, len(results_df)) # Plot up to 5 top sgRNAs
    print(f"  Generating detailed visualizations for top {num_plots} sgRNAs...")
    for i in range(num_plots):
        sgRNA_to_plot = results_df.iloc[i].to_dict()
        plot_filename = os.path.join(output_dir, f"sgRNA_plot_rank_{sgRNA_to_plot['rank']}.html")
        visualization.plot_editing_window(sgRNA_to_plot, editor_spec, plot_filename, output_format="html")

    # Generate summary statistics plots
    print("  Generating summary statistics plots...")
    summary_plot_prefix = os.path.join(output_dir, "summary")
    visualization.plot_summary_statistics(results_df.to_dict('records'), summary_plot_prefix, output_format="html")

    # --- 8. Print Summary to Console ---
    print("\n--- Base Editor Design Assistant Complete ---")
    print(f"Analysis for edit: {args.edit} using editor: {args.editor}")
    print(f"Total compatible sgRNAs found: {len(compatible_sgRNAs)}")
    if not results_df.empty:
        print(f"Top 3 sgRNAs (by efficiency):")
        print(results_df[['rank', 'sequence', 'on_target_efficiency', 'bystander_edits']].head(3).to_string(index=False))
    print(f"\nAll results saved to: {os.path.abspath(output_dir)}")
    print("Please check the generated CSV and HTML files for detailed analysis.")

if __name__ == "__main__":
    main()
