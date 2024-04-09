
.. _plots_get_started:

Get started
=================

Install SEISMOGRAPH 
-------------------

Install SEISMOGRAPH with pip. **Use Python >= 3.10**

.. code::

    pip install seismic-graph


Import your data into a study object
------------------------------------

.. code::

    from seismic_graph import Study
    import json

    # Define the list of samples you want to load
    my_samples = ['my_sample_1', 'my_sample_2', 'my_sample_3']

    # Load the data as a list of dictionaries
    data = []
    for sample in my_samples:
        with open('{}.json'.format(sample), 'r') as f:
            data.append(json.load(f))

    # Create a study object
    my_study = Study(data)

    # Print the list of available samples
    print(my_study.get_samples())


Make a plot
-----------

.. code::

    study.mutation_fraction(
        sample='my_sample_1',
        reference = 'my_reference_1',
        section = 'my_section_1',
        cluster = 'pop_avg'
    )


.. raw:: html
    :file: plots_figs/mutation_fraction.html


.. note::

    We regularly update the list of available plots. Make sure that you have the last version of seismograph.


Special arguments
-----------------

- **pearson_filter_gap**: float, default=1
    Sometimes, a single outlier can lead to a very low Pearson correlation. 
    This parameter allows to filter out the outliers by setting a minimal gap between the Pearson correlation 
    of the data with the outlier and the Pearson correlation of the data without the outlier.
    If the difference is higher than this value, the outlier is removed from the Pearson correlation calculation.

- **normalize**: bool, default=False
    DMS-MaPseq has highly reproducible mutation rates, yet there is often a scale difference between samples.
    This parameter allows to normalize the mutation rates of the samples to the same scale.

    