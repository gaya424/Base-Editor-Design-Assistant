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

# Import necessary libraries for plotting
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.offline as pyo # For saving HTML files
import pandas as pd # For easier data manipulation in summary plots

# Import editor specifications (for descriptive names if needed)
from core.editor_specs import get_editor_spec

# Define a type alias for clarity for the sgRNA candidate dictionary structure.
SgRNACandidate = Dict[str, Any]

def plot_editing_window(
    sgRNA_info: SgRNACandidate,
    editor_spec: Dict[str, Any],
    output_path: str,
    output_format: str = "html"
) -> None:
    """
    Visualizes the sgRNA sequence, highlighting the editable window, target base,
    potential bystander edits, and the PAM sequence.

    Args:
        sgRNA_info (SgRNACandidate): A dictionary containing sgRNA candidate information.
                                     Expected keys: 'sequence', 'target_base_relative_pos',
                                     'pam_sequence', 'bystander_edits'.
        editor_spec (Dict[str, Any]): A dictionary containing the base editor's specifications.
                                      Expected keys: 'editable_window_start', 'editable_window_end',
                                      'description'.
        output_path (str): The full path including filename (e.g., "results/sgRNA_plot.html").
        output_format (str): The desired output format ('html' or 'pdf' - PDF requires additional setup).
    """
    sgRNA_sequence = sgRNA_info["sequence"].upper()
    target_base_relative_pos = sgRNA_info["target_base_relative_pos"]
    pam_sequence = sgRNA_info["pam_sequence"].upper()
    bystander_edits = sgRNA_info.get("bystander_edits", []) # Use .get() for safety
    
    editor_description = editor_spec.get("description", editor_spec.get("name", "Unknown Editor"))
    editable_window_start = editor_spec["editable_window_start"]
    editable_window_end = editor_spec["editable_window_end"]

    # Combine sgRNA sequence and PAM for full context visualization
    # The PAM is typically immediately 3' to the 20bp protospacer.
    full_sequence_display = sgRNA_sequence + pam_sequence
    
    # Create the figure
    fig = go.Figure()

    # Add the full sequence as text annotations
    for i, base in enumerate(full_sequence_display):
        color = 'black'
        if i == target_base_relative_pos:
            color = 'red' # Target base
        elif editable_window_start <= i <= editable_window_end:
            # Check if it's a bystander base (not the primary target but in window)
            is_bystander = False
            for b_edit in bystander_edits:
                if b_edit['relative_pos'] == i:
                    is_bystander = True
                    break
            if is_bystander:
                color = 'orange' # Bystander edit
            else:
                color = 'blue' # In editable window but not target/bystander
        elif i >= len(sgRNA_sequence):
            color = 'purple' # PAM sequence

        fig.add_annotation(
            x=i,
            y=0.5, # Y-position for the base text
            text=base,
            showarrow=False,
            font=dict(size=16, color=color, family="monospace"),
            xanchor='center',
            yanchor='middle'
        )

    # Add shapes for highlighting regions
    # Editable Window (relative to sgRNA)
    fig.add_shape(
        type="rect",
        x0=editable_window_start - 0.5, # -0.5 to center rectangle on base
        y0=0.2, y1=0.8,
        x1=editable_window_end + 0.5,   # +0.5 to center rectangle on base
        fillcolor="rgba(0, 255, 0, 0.2)", # Light green
        line_width=0,
        layer="below",
        name="Editable Window"
    )
    fig.add_annotation(
        x=(editable_window_start + editable_window_end) / 2,
        y=0.9,
        text="Editable Window",
        showarrow=False,
        font=dict(size=10, color="green")
    )

    # PAM Sequence
    pam_display_start = len(sgRNA_sequence)
    pam_display_end = len(full_sequence_display) - 1
    fig.add_shape(
        type="rect",
        x0=pam_display_start - 0.5,
        y0=0.2, y1=0.8,
        x1=pam_display_end + 0.5,
        fillcolor="rgba(255, 165, 0, 0.2)", # Light orange
        line_width=0,
        layer="below",
        name="PAM Sequence"
    )
    fig.add_annotation(
        x=(pam_display_start + pam_display_end) / 2,
        y=0.9,
        text=f"PAM ({pam_sequence})",
        showarrow=False,
        font=dict(size=10, color="orange")
    )

    # Target Base Mark (a distinct shape or line)
    fig.add_shape(
        type="line",
        x0=target_base_relative_pos, y0=0.1,
        x1=target_base_relative_pos, y1=0.9,
        line=dict(color="red", width=2, dash="dot"),
        name="Target Base"
    )
    fig.add_annotation(
        x=target_base_relative_pos,
        y=0.0,
        text="Target",
        showarrow=False,
        font=dict(size=10, color="red")
    )

    # Update layout
    fig.update_layout(
        title=f"sgRNA Visualization for {editor_description}<br><sup>Protospacer: {sgRNA_sequence}</sup>",
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-1, len(full_sequence_display)] # Adjust range to fit all bases
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[0, 1]
        ),
        height=300,
        width=800,
        plot_bgcolor='rgba(0,0,0,0)', # Transparent background
        paper_bgcolor='rgba(0,0,0,0)', # Transparent paper background
        margin=dict(l=20, r=20, t=80, b=20)
    )

    # Save the plot
    if output_format.lower() == "html":
        pyo.plot(fig, filename=output_path, auto_open=False, include_plotlyjs='cdn')
        print(f"  Saved sgRNA visualization to {output_path}")
    elif output_format.lower() == "pdf":
        # Plotly can export to PDF, but it requires kaleido or other renderers.
        # For simplicity, we'll just print a message for PDF as it's not directly supported
        # without external dependencies in a standard Python environment.
        print(f"  PDF output requested for {output_path}. Plotly PDF export requires 'kaleido' library (pip install kaleido). Skipping PDF export for now.")
        # Example if kaleido is installed: fig.write_image(output_path)
    else:
        print(f"  Unsupported output format: {output_format}. Please choose 'html' or 'pdf'.")


def plot_summary_statistics(
    sgRNA_list: List[SgRNACandidate],
    output_path_prefix: str, # e.g., "results/summary"
    output_format: str = "html"
) -> None:
    """
    Generates summary statistics plots for a list of sgRNA candidates.
    Includes GC content distribution, on-target efficiency distribution,
    and bystander edit counts.

    Args:
        sgRNA_list (List[SgRNACandidate]): A list of sgRNA candidate dictionaries.
        output_path_prefix (str): The prefix for output filenames (e.g., "results/summary").
                                  Plots will be saved as <prefix>_gc.html, <prefix>_efficiency.html, etc.
        output_format (str): The desired output format ('html' or 'pdf').
    """
    if not sgRNA_list:
        print("  No sgRNAs provided for summary statistics. Skipping plots.")
        return

    # Extract data for plotting
    gc_contents = []
    on_target_efficiencies = []
    bystander_counts = []
    off_target_warning_status = {"Warnings Present": 0, "No Warnings": 0}

    for sgRNA in sgRNA_list:
        seq = sgRNA["sequence"]
        if seq:
            gc_count = seq.count('G') + seq.count('C')
            gc_contents.append((gc_count / len(seq)) * 100) # Percentage
        
        if sgRNA["on_target_efficiency"] is not None:
            on_target_efficiencies.append(sgRNA["on_target_efficiency"])
        
        bystander_counts.append(len(sgRNA["bystander_edits"]))

        if sgRNA["off_target_warnings"]: # Assuming off_target_warnings is not None/empty if present
            off_target_warning_status["Warnings Present"] += 1
        else:
            off_target_warning_status["No Warnings"] += 1

    # --- Plot 1: GC Content Distribution (Histogram) ---
    if gc_contents:
        fig_gc = go.Figure(data=[go.Histogram(x=gc_contents, nbinsx=10, marker_color='lightblue')])
        fig_gc.update_layout(
            title_text='Distribution of sgRNA GC Content (%)',
            xaxis_title_text='GC Content (%)',
            yaxis_title_text='Number of sgRNAs',
            bargap=0.05
        )
        pyo.plot(fig_gc, filename=f"{output_path_prefix}_gc.html", auto_open=False, include_plotlyjs='cdn')
        print(f"  Saved GC content distribution to {output_path_prefix}_gc.html")

    # --- Plot 2: On-Target Efficiency Distribution (Histogram) ---
    if on_target_efficiencies:
        fig_eff = go.Figure(data=[go.Histogram(x=on_target_efficiencies, nbinsx=10, marker_color='lightcoral')])
        fig_eff.update_layout(
            title_text='Distribution of On-Target Efficiency Scores',
            xaxis_title_text='Efficiency Score (0.0-1.0)',
            yaxis_title_text='Number of sgRNAs',
            bargap=0.05
        )
        pyo.plot(fig_eff, filename=f"{output_path_prefix}_efficiency.html", auto_open=False, include_plotlyjs='cdn')
        print(f"  Saved on-target efficiency distribution to {output_path_prefix}_efficiency.html")

    # --- Plot 3: Off-Target Warning Summary (Pie Chart) ---
    # Only plot if there's at least one category with counts
    if any(off_target_warning_status.values()):
        labels = list(off_target_warning_status.keys())
        values = list(off_target_warning_status.values())
        fig_off_target = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig_off_target.update_layout(
            title_text='Off-Target Warning Summary',
            annotations=[dict(text='Warnings', x=0.5, y=0.5, font_size=12, showarrow=False)]
        )
        pyo.plot(fig_off_target, filename=f"{output_path_prefix}_off_target_summary.html", auto_open=False, include_plotlyjs='cdn')
        print(f"  Saved off-target warning summary to {output_path_prefix}_off_target_summary.html")

    # --- Plot 4: Bystander Edit Counts (Bar Chart) ---
    if bystander_counts:
        # Convert to a Series to easily count occurrences of each bystander count
        bystander_series = pd.Series(bystander_counts)
        bystander_counts_df = bystander_series.value_counts().sort_index().reset_index()
        bystander_counts_df.columns = ['Num_Bystanders', 'Count']

        fig_bystander = go.Figure(data=[go.Bar(
            x=bystander_counts_df['Num_Bystanders'],
            y=bystander_counts_df['Count'],
            marker_color='lightseagreen'
        )])
        fig_bystander.update_layout(
            title_text='Number of Bystander Edits per sgRNA',
            xaxis_title_text='Number of Bystander Edits',
            yaxis_title_text='Number of sgRNAs',
            xaxis=dict(tickmode='linear', tick0=0, dtick=1) # Ensure integer ticks
        )
        pyo.plot(fig_bystander, filename=f"{output_path_prefix}_bystander_counts.html", auto_open=False, include_plotlyjs='cdn')
        print(f"  Saved bystander edit counts to {output_path_prefix}_bystander_counts.html")

    if output_format.lower() == "pdf":
        print(f"  PDF output requested for summary plots. Plotly PDF export requires 'kaleido' library. Skipping PDF export for now.")
