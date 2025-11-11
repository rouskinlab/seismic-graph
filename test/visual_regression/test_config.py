"""Configuration for visual regression testing.

Defines test data paths and parameters for each plot type to be tested.
"""

import os

# Path to test fixtures
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), '..', 'fixtures')

# Test data files (Temperature experiment data)
TEST_DATA_FILES = [
    os.path.join(FIXTURES_DIR, 'GFPIVDPOS109_05x_S34_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'GFPIVDPOS109_7x_S38_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'POS109_0x_S6_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'POS109_02x_S28_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'POS109_005x_S22_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'POS109_007x_S23_L001__webapp.json'),
    os.path.join(FIXTURES_DIR, 'POS109_035x_S25_L001__webapp.json'),
]

# Plot test configurations
# Each entry: (plot_method_name, plot_params_dict, test_description)

PLOT_TESTS = [
    # Basic mutation fraction plot (single sample)
    (
        'mutation_fraction',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'mutation_fraction_basic' 
    ),

    # Mutation fraction identity (single sample)
    (
        'mutation_fraction_identity',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'mutation_fraction_identity_basic'
    ),

    # Experimental variable across samples (ASO_conc ("µM")) - multi-sample
    (
        'experimental_variable_across_samples',
        {
            'experimental_variable': 'ASO_conc ("µM")',
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'exp_var_aso_conc'
    ),

    # Mutation Fraction Delta
    (
        'mutation_fraction_delta',
        {
            'sample' : ['GFPIVDPOS109_05x_S34_L001', 'GFPIVDPOS109_7x_S38_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'mutation_fraction_delta_2x'
    ),

    

    # Mutations per Read per Sample
    (
        'mutations_per_read_per_sample',
        {
            'sample' : ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'mutations_per_read_per_sample'
    ),


    # distribution of aligned reads per reference for each sample and all samples combined
    (
        'num_aligned_reads_per_reference_frequency_distribution',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'num_aligned_reads_per_reference_frequency_distribution'
    ),

    

        # Mutations per Read per Reference
    (
        'mutation_per_read_per_reference',
        {
            'sample' : ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'mutation_per_read_per_reference'
    ),

    # Base coverage with specific positions (single sample)
    (
        'base_coverage',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average'],
            'base_index': [12, 38, 33]
        },
        'base_coverage_selected_positions'
    ),

    # Base coverage all positions (single sample)
    (
        'base_coverage',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'base_coverage_all_positions'
    ),



    # Compare mutation profiles - specific positions (multi-sample)
    (
        'compare_mutation_profiles',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average'],
            'positions_to_compare': [38, 61]
        },
        'compare_profiles_positions_38_61'
    ),

    # Correlation by refs between samples
    (
        'correlation_by_refs_between_samples',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001', 'GFPIVDPOS109_7x_S38_L001'],
            'section': ['full'],
            'cluster': ['average']
        },
        'correlation_by_refs_between_samples'
    ),

    # Distribution of Sequence Lengths
    (
        'dist_of_seq_lengths',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001', 'GFPIVDPOS109_7x_S38_L001'],
            'section': ['full'],
            'cluster': ['average']
        },
        'dist_of_seq_lengths'
    ),


    # Compare mutation profiles - all positions (may be slow, multi-sample)
    (
        'compare_mutation_profiles',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'compare_profiles_all_positions'
    ),

    # One pager (single sample, single row)
    (
        'one_pager',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001'],
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'cluster': ['average']
        },
        'one_pager_reference-pBADupdate'
    ),

    # Binding affinity, no fitting
    (
        'binding_affinity',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'experimental_variable': 'ASO_conc ("µM")',
            'selected_binding_affinity': 'none',
        },
        'binding_affinity_no_fit'
    ),

    # Binding affinity with Hill fitting
    (
        'binding_affinity',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'experimental_variable': 'ASO_conc ("µM")',
            'selected_binding_affinity': 'hill',
        },
        'binding_affinity_hill_fit'
    ),

    # Binding affinity with specific positions, hill fitting
    (
        'binding_affinity',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'experimental_variable': 'ASO_conc ("µM")',
            'selected_binding_affinity': 'hill',
            'positions_to_plot': [175, 208, 800],
        },
        'binding_affinity_selected_positions'
    ),

    # Binding affinity with specific positions, no fitting
    (
        'binding_affinity',
        {
            'reference': ['reference-pBADupdate'],
            'section': ['full'],
            'experimental_variable': 'ASO_conc ("µM")',
            'selected_binding_affinity': 'none',
            'positions_to_plot': [175, 208, 800],
        },
        'binding_affinity_no_fit_selected_positions'
    ),

        # Pearson Correlation Histogram
    (
        'pearson_correlation_histogram',
        {
            'sample': ['GFPIVDPOS109_05x_S34_L001', 'GFPIVDPOS109_7x_S38_L001'],
            'section': ['full'],
        },
        'pearson_correlation_histogram'
    ),
]