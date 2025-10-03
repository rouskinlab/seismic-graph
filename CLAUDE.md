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
@plot_info("binding_affinity", "Binding Affinity")
@save_plot
@doc_inherit(...)
def binding_affinity(self, experimental_variable, selected_binding_affinity, ...):
    """Binding Affinity Plot"""
    return self.wrap_to_plotter(plotter.binding_affinity, locals(), kwargs)
```

The webapp queries `Study.get_plots_list()` to build its UI dynamically.

## Adding a New Plot Function

1. **Create plot function in `plotter.py`**:
   ```python
   def my_new_plot(data: pd.DataFrame, param1, param2, normalize=False) -> dict:
       """Plot description"""
       # Process data
       fig = go.Figure(...)
       return {'fig': fig, 'data': data}
   ```

2. **Add wrapper method to Study class**:
   ```python
   @plot_info("my_new_plot", "My New Plot")
   @save_plot
   @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
   def my_new_plot(self, param1, param2, **kwargs) -> dict:
       """My New Plot description"""
       return self.wrap_to_plotter(plotter.my_new_plot, locals(), kwargs)
   ```

3. **Webapp automatically discovers** the new plot via `get_plots_list()`

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

## Webapp Integration

The separate Flask webapp (`/Users/casper/Local/HMS/Code/draw-app-5/`) uses this library:
- Imports: `from seismic_graph import Study`
- Creates Study objects from user-uploaded data
- Calls Study plot methods with user-selected parameters
- Renders returned Plotly figures in browser
- Uses `build_plot_selections()` in `utils.py` to parse frontend requests

When modifying plots, consider how parameters map to webapp UI components (dropdowns, text fields, etc.).
