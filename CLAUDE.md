# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**seismic-graph** is a Python library for analyzing and visualizing SEISMIC-RNA data. It processes mutation profiles from RNA structure probing experiments and generates interactive Plotly visualizations.

- **Main Package**: `seismic_graph/`
- **Python Version**: >=3.10
- **Documentation**: https://rouskinlab.github.io/seismic-graph
- **Webapp**: Separate Flask webapp at `/Users/casper/Local/HMS/Code/draw-app-5/` that provides a GUI for this library

## Common Commands

```bash
# Install dependencies
pip install -r requirements.txt

# Install package in development mode
python setup.py install
# or
make

# Run tests
pytest test -v
# or
make pytest

# Build documentation
make docs_update

# Clean documentation
make docs_clear

# Upgrade dependencies (careful - updates requirements.txt)
make upgrade-dependencies
```

## Code Architecture

### Core Components

#### 1. Study Class (`seismic_graph/study.py`)
- **Primary user-facing interface** for data analysis
- Wraps all plotting functions via `wrap_to_plotter()` method
- Manages DataFrame of mutation profiles with filtering by sample/reference/section/cluster
- Instantiate with: `Study(data=json_or_df, min_cov=0, filter_by='sample')`
- Key attributes:
  - `df`: pandas DataFrame with mutation profile data
  - `table`: LinFitTable for normalization across samples

#### 2. Plotter Module (`seismic_graph/plotter.py`)
- **Low-level plotting functions** that Study wraps
- Each function signature: `func(data: pd.DataFrame, **kwargs) -> dict`
- Returns: `{'fig': plotly.graph_objects.Figure, 'data': pd.DataFrame}`
- Main plot types:
  - `mutation_fraction()`: Bar chart of mutation rates by position and base
  - `experimental_variable_across_samples()`: Compare samples across experimental variable
  - `binding_affinity()`: Dose-response curves with Hill equation fitting
  - `compare_mutation_profiles()`: Correlation plots between samples
  - `base_coverage()`: Read coverage visualization

#### 3. Manipulator Module (`seismic_graph/manipulator.py`)
- **Data filtering and selection** before plotting
- Core function: `get_df(df, sample=None, reference=None, section=None, cluster=None, base_index=None, base_type=['A','C','G','T'], base_pairing=None, ...)`
- Filters by:
  - Mutation profile attributes (sample, reference, section, cluster)
  - Per-base attributes (base_index, base_type, base_pairing)
  - Minimum coverage (`min_cov`)
- **Important**: `base_index` is 1-indexed for users, converted to 0-indexed internally

#### 4. Utilities (`seismic_graph/util/`)
- `normalization.py`: LinFitTable class for linear regression-based normalization between samples
- `filtered_pearson.py`: Pearson correlation with filtering
- `dump.py`: JSON flattening and DataFrame preprocessing
- `misc.py`: Helper functions for assertions, plot saving, argument extraction

### Data Flow Pattern

```
SEISMIC JSON → flatten_json() → DataFrame → Study.df
                                              ↓
                                    manipulator.get_df() (filters data)
                                              ↓
                                    plotter.function() (generates plot)
                                              ↓
                            {'fig': Plotly Figure, 'data': DataFrame}
```

### Study.wrap_to_plotter() Architecture

The `wrap_to_plotter()` method is central to the Study class design:

```python
def wrap_to_plotter(self, func, loc, kwargs):
    # 1. Merge local variables and kwargs
    # 2. Call manipulator.get_df() with relevant filtering parameters
    # 3. Call plotter function with remaining parameters
    # 4. Return result
```

This pattern allows Study methods to:
- Accept combined filtering + plotting parameters
- Automatically route parameters to the right functions
- Maintain a clean user-facing API

### Plot Registration System

Plot functions in Study class use `@plot_info` decorator for webapp integration:

```python
@plot_info("binding_affinity", "Binding Affinity",
           requirements={"rows": {"min": 1}, "references": {"min": 1, "max": 1}, "columns": ["experimental_variable"]})
@save_plot
@doc_inherit(...)
def binding_affinity(self, experimental_variable, selected_binding_affinity, ...):
    """Binding Affinity Plot"""
    return self.wrap_to_plotter(plotter.binding_affinity, locals(), kwargs)
```

The webapp queries `Study.get_plots_list()` to build its UI dynamically.

### Plot Requirements System

Each plot can specify **data selection requirements** that must be met before the plot can be generated. Requirements are defined in the `@plot_info` decorator and automatically:
- Sent to the webapp in the `plot_options` API response
- Displayed in the webapp's "Choose Plot Type" panel with visual indicators (✓/✗)
- Used to enable/disable the Plot button based on current data selections

#### Requirements Data Structure

```python
requirements={
    "rows": {"min": 1, "max": 1},           # Number of DataFrame rows required
    "samples": {"min": 2, "max": 2},        # Number of unique samples required
    "references": {"min": 1, "max": 1},     # Number of unique references required
    "columns": ["experimental_variable"]     # Required columns in the data
}
```

**Constraint Format:**
- `{"min": 1, "max": 1}` - exactly one
- `{"min": 2, "max": 2}` - exactly two
- `{"min": 1}` or `{"min": 1, "max": None}` - one or more (unbounded)
- `{"min": 2, "max": 14}` - between 2 and 14 inclusive

#### Complete Requirements Reference

| Plot Function | Rows | Samples | References | Columns |
|--------------|------|---------|------------|---------|
| mutation_fraction | min:1, max:1 | - | - | - |
| mutation_fraction_identity | min:1, max:1 | - | - | - |
| mutation_fraction_delta | - | min:2, max:2 | min:1, max:1 | - |
| base_coverage | min:1, max:1 | - | - | - |
| mutation_per_read_per_reference | min:1, max:1 | - | - | - |
| one_pager | min:1, max:1 | - | - | - |
| correlation_by_refs_between_samples | min:1 | min:2, max:2 | - | - |
| pearson_correlation_histogram | min:1 | min:2, max:2 | - | - |
| experimental_variable_across_samples | min:1 | - | min:1, max:1 | experimental_variable |
| binding_affinity | min:1 | - | min:1, max:1 | experimental_variable |
| compare_mutation_profiles | min:2, max:14 | - | min:1, max:1 | - |
| mutations_per_read_per_sample | min:1 | - | - | - |
| num_aligned_reads_per_reference_frequency_distribution | min:1 | - | - | - |

**Key Design Principles:**
- **Rows**: Number of rows after filtering (each row = unique sample/reference/section/cluster combo)
- **Samples**: Number of unique sample values in the filtered data
- **References**: Number of unique reference values (reference = unique RNA sequence)
- **Columns**: Required columns like `experimental_variable` that must exist in the data
- **Empty requirements `{}`**: No constraints, plot always valid (if requirements parameter omitted)

#### Frontend Validation (Webapp)

The webapp frontend performs client-side validation using simple range checking:
```javascript
// Check if current value meets requirement
const isValid = (currentValue >= min) && (currentValue <= max);
```

This enables the webapp to:
- Disable the Plot button when requirements aren't met
- Show real-time visual feedback (green ✓ or black ✗) for each requirement
- Display current vs. expected values (e.g., "Exactly 2 (current: 1)")
- Prevent invalid plot requests from reaching the backend

## Adding a New Plot Function

1. **Create plot function in `plotter.py`**:
   ```python
   def my_new_plot(data: pd.DataFrame, param1, param2, normalize=False) -> dict:
       """Plot description"""
       # Process data
       fig = go.Figure(...)
       return {'fig': fig, 'data': data}
   ```

2. **Add wrapper method to Study class with requirements**:
   ```python
   @plot_info("my_new_plot", "My New Plot",
              requirements={"rows": {"min": 1}, "references": {"min": 1, "max": 1}})
   @save_plot
   @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
   def my_new_plot(self, param1, param2, **kwargs) -> dict:
       """My New Plot description"""
       return self.wrap_to_plotter(plotter.my_new_plot, locals(), kwargs)
   ```

3. **Determine requirements**:
   - Look for validation assertions in the plotter function (e.g., `assert len(df) == 1`)
   - Consider what data structure the plot logically needs
   - Use min/max format: `{"min": 1, "max": 1}` for exactly one, `{"min": 1}` for one or more
   - Add `columns` requirement if plot needs specific columns like `experimental_variable`

4. **Webapp automatically discovers** the new plot and its requirements via `get_plots_list()`

## Important Implementation Details

### Position Indexing
- **User-facing APIs**: 1-indexed positions (genomic convention)
- **Internal processing**: 0-indexed (Python convention)
- **Conversion**: `manipulator.get_df()` converts `base_index` from 1-indexed to 0-indexed

### Performance Optimization Pattern (binding_affinity example)
When users specify a subset of positions via `positions_to_plot`:
1. Filter positions early in DataFrame building: `build_dataframes_from_replicates(positions_to_plot=...)`
2. Only fit curves for those positions: `fit_hill_per_position(positions_to_fit=...)`
3. Avoids expensive computation on unwanted data

**Key principle**: Pass position filters down through the entire pipeline, not just at visualization.

### Experimental Variable Generalization
Plot functions like `binding_affinity()` work with any numeric experimental variable (concentration, temperature, pH, etc.):
- Extract unit from variable name: `'ASO_conc ("µM")'` → `{"name": "ASO_conc", "unit": "µM", "has_unit": True}`
- Use generic terminology: `exp_var`, `exp_value` (not `concentration`, `conc`)
- Make labels and hover text unit-aware: `f"{exp_value}{unit_str}"`

### Replicate Handling
When multiple samples have the same experimental variable value (replicates):
- Keep all replicates as separate columns: `"{value} ({sample_name})"`
- Fit Hill curves to all individual points (don't average)
- Track sample names in column labels for hover text

### DataFrame Schema
Key columns in Study.df:
- `sample`, `reference`, `section`, `cluster`: Identifiers (a row = unique combination)
- `sequence`: String of bases (e.g., "ACGTACGT")
- `structure`: Dot-bracket notation for secondary structure
- `sub_rate`: Numpy array of mutation rates per position
- `cov`: Numpy array of coverage per position
- `num_aligned`: Total aligned reads
- `min_cov`: Minimum coverage across positions (for filtering)

## Input Data Format

### JSON File Structure

seismic-graph processes JSON files formatted as `*__webapp.json` containing mutation profile data from SEISMIC-RNA experiments.

#### Hierarchical Data Organization

Data is nested using four levels: **Sample → Reference → Section → Cluster**

```
{
  "#sample": "sample_name",           # Top-level metadata
  "#User": "initials",
  "#Date": "YYYY-MM-DD",
  ... (other metadata keys starting with #)

  "reference_name_1": {               # Reference level (e.g., "C5", "C10")
    "#sequence": "ACGT...",
    "#num_aligned": 6026,

    "section_name": {                 # Section level (typically "full")
      "#section_start": 1,
      "#section_end": 67,
      "#positions": [1, 2, 3, ...],   # For backwards compatibility

      "cluster_name": {               # Cluster level (typically "average")
        "sub_rate": [...],            # Core data: mutation rates per position
        "cov": [...],                 # Core data: coverage per position
        "sub_hist": [1146, 1712, 805],
        "sub_N": [...],
        "sub_A": [...],
        "sub_C": [...],
        "sub_G": [...],
        "sub_T": [...],
        "del": [...],
        "ins": [...]
      }
    }
  },

  "reference_name_2": { ... }
}
```

#### Key Structure Rules

1. **Metadata Keys**: All keys starting with `#` are metadata
2. **Reference Names**: Any top-level key without `#` is a reference name
3. **Guaranteed Levels**: At minimum, expect section `"full"` and cluster `"average"`
4. **Multiple Entities**: One file = one sample; may contain multiple references, sections per reference, and clusters per section

#### Core Data Fields (used by seismic-graph)

- `sub_rate`: Array of mutation rates per position (primary analysis data)
- `cov`: Array of read coverage per position
- `#num_aligned`: Total number of aligned reads (reference level)
- `#sequence`: RNA sequence string (reference level)
- `#positions`: List from `range(#section_start, #section_end+1)` - for backwards compatibility

#### Position Indexing in JSON

- `#section_start` and `#section_end` use 1-indexed positions (genomic convention)
- `#positions` array contains 1-indexed positions
- Data arrays (`sub_rate`, `cov`, etc.) use 0-indexed Python arrays
- User-facing APIs convert from 1-indexed to 0-indexed internally

## Webapp Integration

The separate Flask webapp (`/Users/casper/Local/HMS/Code/draw-app-5/`) uses this library:
- Imports: `from seismic_graph import Study`
- Creates Study objects from user-uploaded data
- Calls Study plot methods with user-selected parameters
- Renders returned Plotly figures in browser
- Uses `build_plot_selections()` in `utils.py` to parse frontend requests

When modifying plots, consider how parameters map to webapp UI components (dropdowns, text fields, etc.).
