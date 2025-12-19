import re
import matplotlib.patches as mpatches
from matplotlib.path import Path
from spectrum_utils.spectrum import MsmsSpectrum

def draw_peptide_annotation(
    spec: MsmsSpectrum,
    draw_peptide_kws: dict,
    ax,
):
    """
    Draw peptide sequence with y and b ion annotations.
    
    Parameters:
    -----------
    spec : MsmsSpectrum
        The spectrum object containing proforma and annotation data.
    draw_peptide_kws : dict
        Dictionary of keyword arguments for customizing the peptide drawing.
        Supported keys:
        - font_size (int): Font size for amino acids (default: 12)
        - show_ptms (bool): Whether to show PTMs in the sequence (default: False)
    ax : matplotlib.axes.Axes
        The axis to draw on
    """
    # Extract proforma string
    if not hasattr(spec, "proforma") or spec.proforma is None:
        return  # Nothing to draw if no proforma annotation
    
    proforma = spec.proforma
    if isinstance(proforma, list):
        if len(proforma) == 0:
            return
        proforma_str = str(proforma[0])
    else:
        proforma_str = str(proforma)
    
    # Get annotations
    annotations = spec.annotation

    
    # Extract drawing parameters
    font_size = draw_peptide_kws.get("font_size", 12)
    show_ptms = draw_peptide_kws.get("show_ptms", False)
    
    # Extract sequence from proforma string

    
    # Tokenize the sequence to handle modifications
    tokens = _tokenize_mod_sequence(proforma_str, show_ptms=show_ptms)
    
    # Force a draw first
    ax.figure.canvas.draw()
    
    # Get figure dimensions
    fig_width_inches = ax.figure.get_figwidth()
    fig_height_inches = ax.figure.get_figheight()
    
    # Measure token widths using visible text temporarily
    temp_text = ax.text(0.5, 0.5, ' ,', 
                        transform=ax.figure.transFigure,
                        fontfamily='monospace',
                        fontsize=font_size,
                        fontweight='bold',
                        alpha=0)
    ax.figure.canvas.draw()
    bbox = temp_text.get_window_extent(renderer=ax.figure.canvas.get_renderer())
    space_width = bbox.width

    token_widths_display = []
    for token in tokens:
        temp_text = ax.text(0.5, 0.5, token, 
                           transform=ax.figure.transFigure,
                           fontfamily='monospace',
                           fontsize=font_size,
                           fontweight='bold',
                           alpha=0)
        ax.figure.canvas.draw()
        bbox = temp_text.get_window_extent(renderer=ax.figure.canvas.get_renderer())
        token_widths_display.append(bbox.width + space_width)
    
    # Get figure width in pixels
    fig_width_px = ax.figure.get_window_extent().width
    
    # Convert each pixel width to axes fraction
    token_widths = [w / fig_width_px for w in token_widths_display]
    
    # Calculate total width and starting position
    total_width = sum(token_widths)
    peptide_x_start = 0.5 - (total_width / 2)
    
    # Scale char_height based on figure height
    char_height = (font_size / 72.0) / fig_height_inches
    peptide_y = 0.98 - (char_height * 1.5)
    
    # Get average char width for hook calculations
    avg_char_width = sum(token_widths) / len(tokens) if tokens else 0

    # Draw amino acids and track their positions
    token_positions = []
    current_x = peptide_x_start
    
    for i, token in enumerate(tokens):
        token_width = token_widths[i]
        x_center = current_x + (token_width / 2)
        ax.text(x_center, peptide_y, token, 
                transform=ax.transAxes,
                ha='center', va='center',
                fontfamily='monospace',
                fontsize=font_size,
                fontweight='bold')
        
        token_positions.append((current_x, x_center, current_x + token_width))
        current_x += token_width
    
    # Filter for b and y ions
    b_ions = []
    y_ions = []
    
    if annotations is not None:
        for peak_interp in annotations:
            if peak_interp is not None and hasattr(peak_interp, 'fragment_annotations'):
                for frag_ann in peak_interp.fragment_annotations:
                    if hasattr(frag_ann, 'ion_type') and frag_ann.ion_type is not None:
                        ion_type_str = frag_ann.ion_type
                        #test
                        # Check if it's a b or y ion
                        if len(ion_type_str) > 0 and ion_type_str[0] in 'by':
                            ion_letter = ion_type_str[0]
                            
                            # Extract numeric part (digits only, ignoring charge, neutral loss, etc.)
                            try:
                                ion_index = int(''.join(c for c in ion_type_str[1:] if c.isdigit()))
                                if ion_letter == 'b':
                                    b_ions.append(ion_index)
                                elif ion_letter == 'y':
                                    y_ions.append(ion_index)
                            except (ValueError, IndexError):
                                continue
    
    # Remove duplicates and sort
    b_ions = sorted(set(b_ions))
    y_ions = sorted(set(y_ions))
    
    # Scale linewidth based on font size
    linewidth = font_size * 0.15
    # Draw b-ion separators (top, blue)
    for b_num in b_ions:
        if 1 <= b_num < len(tokens):
            x_right = token_positions[b_num - 1][2]
            x_hook = x_right - (avg_char_width * 0.2)
            
            verts = [
                (x_hook, peptide_y + (char_height * 1.2)),
                (x_right, peptide_y + (char_height * 0.7)),
                (x_right, peptide_y + (char_height * 0.05))
            ]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO]
            path = Path(verts, codes)
            patch = mpatches.PathPatch(path, facecolor='none', 
                                    edgecolor='#1f77b4', 
                                    linewidth=linewidth,
                                    transform=ax.transAxes)
            ax.add_patch(patch)

    # Draw y-ion separators (bottom, red)
    for y_num in y_ions:
        if 1 <= y_num < len(tokens):
            i = len(tokens) - y_num
            x_left = token_positions[i - 1][2]
            x_hook = x_left + (avg_char_width * 0.2)
            
            verts = [
                (x_left, peptide_y - (char_height * 0.05)),
                (x_left, peptide_y - (char_height * 0.7)),
                (x_hook, peptide_y - (char_height * 1.2))
            ]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO]
            path = Path(verts, codes)
            patch = mpatches.PathPatch(path, facecolor='none', 
                                      edgecolor='#d62728', 
                                      linewidth=linewidth,
                                      transform=ax.transAxes)
            ax.add_patch(patch)


def _tokenize_mod_sequence(sequence: str, show_ptms: bool = False):
    """
    Tokenize a peptide sequence into individual amino acids and modifications.
    
    Parameters:
    -----------
    sequence : str
        The peptide sequence string (can include modifications in ProForma format)
    show_ptms : bool
        Whether to include PTM annotations in the tokens
        
    Returns:
    --------
    list
        List of tokens (amino acids and optionally their modifications)
    """
    sequence = sequence.split('/')[0]  # Remove charge state if present
    sequence = sequence.replace('-', '')  # Remove any '-' characters
    
    if not show_ptms:
        # Remove all [mod] and {mod} blocks and their contents
        sequence = re.sub(r'\[.*?\]', '', sequence)
        sequence = re.sub(r'\{.*?\}', '', sequence)
    
    tokens = []
    i = 0
    
    while i < len(sequence):
        if i == 0 and sequence[i] == '[':
            # N-terminal modification
            mod_end = sequence.find(']', i)
            if mod_end != -1 and mod_end + 1 < len(sequence):
                tokens.append(sequence[i:mod_end+2])
                i = mod_end + 2
            else:
                i += 1
        elif sequence[i] == '[':
            # Attach modification to previous token
            mod_end = sequence.find(']', i)
            if mod_end != -1:
                if tokens:
                    tokens[-1] += sequence[i:mod_end+1]
                else:
                    tokens.append(sequence[i:mod_end+1])
                i = mod_end + 1
            else:
                i += 1
        else:
            tokens.append(sequence[i])
            i += 1
    
    return tokens