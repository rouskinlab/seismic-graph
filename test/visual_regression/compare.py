#!/usr/bin/env python3
"""Compare baseline and current visual regression outputs.

This script helps compare plots generated before and after code changes.
It can open HTML files side-by-side for manual review, or perform automated
JSON structure comparisons.

Usage:
    python compare.py list         # List all test cases
    python compare.py open <test>  # Open baseline and current for a specific test
    python compare.py open-all     # Open all tests in browser
    python compare.py diff         # Show JSON structure differences
"""

import os
import sys
import json
import argparse
import webbrowser
from pathlib import Path
from difflib import unified_diff


def get_test_files(baseline_dir, current_dir):
    """Get list of test files that exist in both directories."""
    baseline_path = Path(baseline_dir)
    current_path = Path(current_dir)

    if not baseline_path.exists():
        print(f"Error: Baseline directory not found: {baseline_dir}")
        print("Run: python test_runner.py baseline")
        return []

    if not current_path.exists():
        print(f"Error: Current directory not found: {current_dir}")
        print("Run: python test_runner.py current")
        return []

    baseline_files = set(f.name for f in baseline_path.glob("*.html"))
    current_files = set(f.name for f in current_path.glob("*.html"))

    common_files = sorted(baseline_files & current_files)
    missing_in_current = baseline_files - current_files
    new_in_current = current_files - baseline_files

    return {
        'common': common_files,
        'missing_in_current': sorted(missing_in_current),
        'new_in_current': sorted(new_in_current),
    }


def list_tests(baseline_dir, current_dir):
    """List all available test cases."""
    files = get_test_files(baseline_dir, current_dir)

    if not files['common']:
        print("No common test files found.")
        return

    print("=== Available Test Cases ===")
    for i, filename in enumerate(files['common'], 1):
        test_name = filename.replace('.html', '')
        print(f"{i}. {test_name}")

    if files['missing_in_current']:
        print("\nMissing in current (only in baseline):")
        for filename in files['missing_in_current']:
            print(f"  - {filename}")

    if files['new_in_current']:
        print("\nNew in current (not in baseline):")
        for filename in files['new_in_current']:
            print(f"  - {filename}")


def open_test(baseline_dir, current_dir, test_name):
    """Open baseline and current versions of a test in the browser."""
    baseline_file = Path(baseline_dir) / f"{test_name}.html"
    current_file = Path(current_dir) / f"{test_name}.html"

    if not baseline_file.exists():
        print(f"Error: Baseline file not found: {baseline_file}")
        return False

    if not current_file.exists():
        print(f"Error: Current file not found: {current_file}")
        return False

    print(f"Opening {test_name}...")
    print(f"  Baseline: {baseline_file}")
    print(f"  Current:  {current_file}")

    webbrowser.open(f"file://{baseline_file.absolute()}")
    webbrowser.open(f"file://{current_file.absolute()}")

    return True


def open_all_tests(baseline_dir, current_dir):
    """Open all test cases in the browser."""
    files = get_test_files(baseline_dir, current_dir)

    if not files['common']:
        print("No test files to open.")
        return

    print(f"Opening {len(files['common'])} test cases...")
    print("(This will open many browser tabs!)")

    response = input("Continue? [y/N]: ")
    if response.lower() != 'y':
        print("Cancelled.")
        return

    for filename in files['common']:
        test_name = filename.replace('.html', '')
        open_test(baseline_dir, current_dir, test_name)


def extract_plot_json(html_file):
    """Extract Plotly JSON data from HTML file."""
    with open(html_file, 'r') as f:
        content = f.read()

    # Find the Plotly.newPlot call and extract JSON
    import re

    # Helper function to extract balanced brackets/braces
    def extract_balanced_structure(text, start_pos, start_char):
        """Extract a balanced JSON structure starting at start_pos."""
        end_char = ']' if start_char == '[' else '}'
        count = 0
        in_string = False
        escape = False

        for i in range(start_pos, len(text)):
            char = text[i]

            # Handle string literals (to ignore brackets inside strings)
            if char == '"' and not escape:
                in_string = not in_string
            elif char == '\\' and not escape:
                escape = True
                continue

            if not in_string:
                if char == start_char:
                    count += 1
                elif char == end_char:
                    count -= 1
                    if count == 0:
                        # Found matching bracket, try to parse JSON
                        json_str = text[start_pos:i+1]
                        try:
                            return json.loads(json_str), i+1
                        except json.JSONDecodeError as e:
                            return None, -1

            escape = False

        return None, -1

    # Look for Plotly.newPlot with the pattern:
    # Plotly.newPlot("div-id", [data], {layout}, {config})
    match = re.search(r'Plotly\.newPlot\s*\(\s*"[^"]*"\s*,\s*', content)
    if not match:
        return None

    # Position after the div ID and comma
    pos = match.end()

    # Extract data array
    # Skip whitespace to find the opening bracket
    while pos < len(content) and content[pos] in ' \t\n\r':
        pos += 1

    if pos >= len(content) or content[pos] != '[':
        return None

    data, pos = extract_balanced_structure(content, pos, '[')
    if data is None:
        return None

    # Skip whitespace and comma to find layout
    while pos < len(content) and content[pos] in ' \t\n\r,':
        pos += 1

    if pos >= len(content) or content[pos] != '{':
        # No layout found, just return data
        return {'data': data, 'layout': {}}

    layout, _ = extract_balanced_structure(content, pos, '{')
    if layout is None:
        return {'data': data, 'layout': {}}

    return {'data': data, 'layout': layout}


def diff_plots(baseline_dir, current_dir):
    """Show structural differences between baseline and current plots."""
    files = get_test_files(baseline_dir, current_dir)

    if not files['common']:
        print("No test files to compare.")
        return

    print("=== Comparing Plot Structures ===")
    differences_found = False

    for filename in files['common']:
        test_name = filename.replace('.html', '')
        baseline_file = Path(baseline_dir) / filename
        current_file = Path(current_dir) / filename

        baseline_json = extract_plot_json(baseline_file)
        current_json = extract_plot_json(current_file)

        if not baseline_json or not current_json:
            print(f"\n{test_name}: Could not extract JSON")
            continue

        # Compare JSON structures
        baseline_str = json.dumps(baseline_json, indent=2, sort_keys=True)
        current_str = json.dumps(current_json, indent=2, sort_keys=True)

        if baseline_str != current_str:
            differences_found = True
            print(f"\n{test_name}: DIFFERENCES FOUND")

            # Show diff
            diff = unified_diff(
                baseline_str.splitlines(),
                current_str.splitlines(),
                fromfile='baseline',
                tofile='current',
                lineterm=''
            )

            # Show first 50 lines of diff
            diff_lines = list(diff)
            for line in diff_lines[:50]:
                if line.startswith('+'):
                    print(f"  \033[32m{line}\033[0m")  # Green
                elif line.startswith('-'):
                    print(f"  \033[31m{line}\033[0m")  # Red
                else:
                    print(f"  {line}")

            if len(diff_lines) > 50:
                print(f"  ... ({len(diff_lines) - 50} more lines)")

        else:
            print(f"{test_name}: ✓ Identical")

    if not differences_found:
        print("\n✓ All plots have identical JSON structure!")


def approve_current(baseline_dir, current_dir):
    """Promote current outputs to baseline (after manual review)."""
    files = get_test_files(baseline_dir, current_dir)

    total_files = len(files['common']) + len(files['new_in_current'])

    if total_files == 0:
        print("No current files to approve.")
        return

    print("=== Approve Current as New Baseline ===")
    print(f"This will update {len(files['common'])} existing baseline files")
    print(f"and add {len(files['new_in_current'])} new baseline files.")

    if files['new_in_current']:
        print("\nNew files to be added:")
        for filename in files['new_in_current']:
            print(f"  + {filename}")

    print("\nMake sure you have manually reviewed the current outputs!")
    print()

    response = input("Continue? [y/N]: ")
    if response.lower() != 'y':
        print("Cancelled.")
        return

    baseline_path = Path(baseline_dir)
    current_path = Path(current_dir)

    import shutil

    # Update existing baseline files
    for filename in files['common']:
        current_file = current_path / filename
        baseline_file = baseline_path / filename
        shutil.copy2(current_file, baseline_file)
        print(f"✓ Updated {filename}")

    # Add new baseline files
    for filename in files['new_in_current']:
        current_file = current_path / filename
        baseline_file = baseline_path / filename
        shutil.copy2(current_file, baseline_file)
        print(f"✓ Added {filename}")

    print()
    print(f"Successfully updated {len(files['common'])} and added {len(files['new_in_current'])} baseline files!")
    print("You can now commit the updated baseline files to git.")


def main():
    script_dir = Path(__file__).parent
    baseline_dir = script_dir / 'baseline'
    current_dir = script_dir / 'current'

    parser = argparse.ArgumentParser(
        description='Compare visual regression test outputs'
    )
    parser.add_argument(
        'command',
        choices=['list', 'open', 'open-all', 'diff', 'approve'],
        help='Command to run'
    )
    parser.add_argument(
        'test_name',
        nargs='?',
        help='Test name (for "open" command)'
    )

    args = parser.parse_args()

    if args.command == 'list':
        list_tests(baseline_dir, current_dir)

    elif args.command == 'open':
        if not args.test_name:
            print("Error: Test name required for 'open' command")
            print("Usage: python compare.py open <test_name>")
            sys.exit(1)
        open_test(baseline_dir, current_dir, args.test_name)

    elif args.command == 'open-all':
        open_all_tests(baseline_dir, current_dir)

    elif args.command == 'diff':
        diff_plots(baseline_dir, current_dir)

    elif args.command == 'approve':
        approve_current(baseline_dir, current_dir)


if __name__ == '__main__':
    main()
