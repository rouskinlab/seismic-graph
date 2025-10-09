#!/usr/bin/env python3
"""Visual regression test runner for seismic-graph.

Generates plots for all configured test cases and saves them as HTML files.
Run this script before and after code changes to compare visual outputs.

Usage:
    python test_runner.py baseline  # Generate baseline outputs
    python test_runner.py current   # Generate current outputs for comparison
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Add parent directory to path to import seismic_graph
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

import plotly.io as pio
from seismic_graph import Study
from test_config import (
    TEST_DATA_FILES,
    PLOT_TESTS,
)


def load_test_data(data_files):
    """Load test data from JSON files."""
    data = []
    for path in data_files:
        if not os.path.exists(path):
            print(f"Warning: Test data file not found: {path}")
            continue
        with open(path, 'r') as f:
            data.append(json.load(f))
    return data


def generate_plots(output_dir, study, test_configs):
    """Generate all plots and save as HTML files.

    Args:
        output_dir: Directory to save HTML files
        study: Study object with loaded data
        test_configs: List of (method_name, params, description) tuples
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    results = {
        'success': [],
        'failed': [],
    }

    for method_name, params, description in test_configs:
        print(f"Generating {description}...", end=' ')

        try:
            # Get the plot method from Study object
            plot_method = getattr(study, method_name)

            # Generate the plot
            result = plot_method(**params)

            # Save as HTML
            output_file = output_path / f"{description}.html"

            # Handle different return formats
            if 'fig' in result:
                # Standard format: {'fig': plotly.Figure, 'data': ...}
                pio.write_html(result['fig'], str(output_file))
            elif 'html' in result:
                # one_pager format: {'html': html_string, 'data': ...}
                with open(output_file, 'w') as f:
                    f.write(result['html'])
            else:
                raise ValueError(f"Unknown result format: {result.keys()}")

            print(f"✓ Saved to {output_file}")
            results['success'].append(description)

        except Exception as e:
            print(f"✗ Error: {str(e)}")
            results['failed'].append((description, str(e)))

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Generate visual regression test outputs for seismic-graph'
    )
    parser.add_argument(
        'output_type',
        choices=['baseline', 'current'],
        help='Type of output to generate (baseline = known good, current = test run)'
    )

    args = parser.parse_args()

    # Determine output directory
    script_dir = Path(__file__).parent
    output_dir = script_dir / args.output_type

    print(f"=== Visual Regression Test Runner ===")
    print(f"Output directory: {output_dir}")
    print()

    # Generate standard plots
    print("Loading test data...")
    data = load_test_data(TEST_DATA_FILES)

    if not data:
        print("Error: No test data could be loaded!")
        sys.exit(1)

    print(f"Loaded {len(data)} test data files")

    print("Creating Study object...")
    study = Study(data=data)
    print(f"Study contains {len(study.df)} rows")
    print()

    print(f"Generating {len(PLOT_TESTS)} standard plots...")
    results = generate_plots(output_dir, study, PLOT_TESTS)


    # Print summary
    print()
    print("=== Summary ===")
    print(f"✓ Successful: {len(results['success'])}")
    print(f"✗ Failed: {len(results['failed'])}")

    if results['failed']:
        print()
        print("Failed tests:")
        for description, error in results['failed']:
            print(f"  - {description}: {error}")
        sys.exit(1)

    print()
    print(f"All plots saved to: {output_dir}")
    print(f"Open {output_dir}/*.html to view the plots")


if __name__ == '__main__':
    main()
