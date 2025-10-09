# Visual Regression Testing for seismic-graph

Automated visual regression testing system to ensure plot outputs remain consistent across code changes.

## Overview

This testing framework automates the testing process. It generates all plot types, saves them as HTML files, and provides tools to compare outputs before and after code changes.

## Quick Start

**Prerequisites:** Activate the `seismic-graph` conda environment before running tests:
```bash
conda activate seismic-graph
```

### 1. Generate Baseline (Before Making Changes)

```bash
make visual-baseline
```

This creates "known good" baseline plots in `test/visual_regression/baseline/`.

### 2. Make Your Code Changes

Edit the seismic-graph codebase as needed.

### 3. Generate Current Test Output

```bash
make visual-test
```

This generates new plots with your changes in `test/visual_regression/current/`.

### 4. Compare Outputs

**Option A: Visual Comparison (Recommended)**
```bash
make visual-compare-open
```
Opens all baseline and current plots side-by-side in your browser for manual review.

**Option B: List Tests**
```bash
make visual-compare
```
Lists all available test cases.

**Option C: JSON Structure Diff**
```bash
make visual-compare-diff
```
Shows structural differences in the plot JSON (useful for catching data changes).

### 5. Approve Changes (If Everything Looks Good)

```bash
make visual-approve
```

This promotes the current outputs to be the new baseline. Commit the updated baseline files to git.

## Test Configuration

Edit [test_config.py](test_config.py) to:
- Add/remove plot test cases
- Change plot parameters
- Use different test datasets

## Directory Structure

```
test/
├── fixtures/                        # Test data files
│   ├── 20200117_IVT4_10C_trimmed__webapp.json
│   ├── 20200117_IVT4_30C_trimmed__webapp.json
│   └── 20200117_IVT4_42C_trimmed__webapp.json
└── visual_regression/
    ├── baseline/                    # "Known good" baseline plots (git tracked)
    │   ├── mutation_fraction_basic.html
    │   ├── base_coverage_all_positions.html
    │   └── ...
    ├── current/                     # Latest test run (git ignored)
    │   └── ...
    ├── test_config.py              # Test case configuration
    ├── test_runner.py              # Plot generation script
    ├── compare.py                  # Comparison utilities
    └── README.md                   # This file
```

## Command Reference

### Make Targets

```bash
make visual-help            # Show help message with all commands
make visual-baseline        # Generate baseline plots
make visual-test           # Generate current test plots
make visual-compare        # List all test cases
make visual-compare-open   # Open all tests in browser
make visual-compare-diff   # Show JSON structure differences
make visual-approve        # Approve current as new baseline
```

### Direct Script Usage

**Generate plots:**
```bash
python test/visual_regression/test_runner.py baseline
python test/visual_regression/test_runner.py current
```

**Compare outputs:**
```bash
python test/visual_regression/compare.py list                    # List tests
python test/visual_regression/compare.py open <test_name>        # Open specific test
python test/visual_regression/compare.py open-all                # Open all tests
python test/visual_regression/compare.py diff                    # Show differences
python test/visual_regression/compare.py approve                 # Approve current
```

## Typical Workflows

### Before Making Changes
```bash
# Generate baseline if you haven't already
make visual-baseline
```

### After Making Changes
```bash
# Test your changes
make visual-test

# Review visually
make visual-compare-open

# If everything looks good
make visual-approve
git add test/visual_regression/baseline/
git commit -m "Update visual regression baselines"
```

### Testing Specific Changes
```bash
# Generate current plots
make visual-test

# Open just one test for review
python test/visual_regression/compare.py open mutation_fraction_basic

# Check for structural changes
make visual-compare-diff
```

## Tips

1. **Commit baselines to git**: The `baseline/` directory should be tracked so you can see plot changes in code review.

2. **Ignore current/**: The `current/` directory is for temporary comparison and doesn't need to be committed.

3. **Review regularly**: After any significant changes to plotter.py, manipulator.py, or study.py, run the visual tests.

4. **Add new tests**: When adding a new plot function, add test cases to `test_config.py`.

5. **Use position filters**: For large datasets, use `positions_to_plot` or `base_index` parameters to speed up tests.

## Limitations

- **Manual review required**: This automates generation but still requires human eyes to verify plots look correct.
- **No pixel-perfect comparison**: HTML/Plotly outputs may have minor rendering differences. Focus on data accuracy.
- **Test data size**: Keep fixture data reasonably small for fast test runs.

## Future Enhancements

Potential improvements:
- Automated image diffing (pixel comparison)
- Automated JSON value assertions for key statistics
- Integration with pytest
- CI/CD integration to catch regressions automatically
- Performance benchmarking alongside visual testing

## Troubleshooting

**"No test data could be loaded"**
- Check that files exist in `test/fixtures/`
- Verify paths in `test_config.py`

**"Baseline directory not found"**
- Run `make visual-baseline` first

**Plot generation fails**
- Check that seismic-graph is properly installed: `python setup.py install`
- Verify test data is valid JSON in the expected format

**Browser doesn't open**
- The files are saved locally; you can manually open `test/visual_regression/baseline/*.html` and `current/*.html`
