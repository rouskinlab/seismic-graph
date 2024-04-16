
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


Outlier filtering
-----------------

- **pearson_filter_gap**: float, default=None
    This parameter allows filtering out outliers by setting a minimal gap between the Pearson correlation of the data 
    with the outlier and the Pearson correlation of the data without the outlier. If the difference is higher than 
    this value, the outlier is removed from the Pearson correlation calculation. 

    If the parameter is set to None (default), no outlier filtering is applied.

    Example:

    .. code::

        reference_data = [1.1, 1.8, 2.9, 4.1, 5.2, 6.0]
        data_with_outlier = [1., 2., 3., 4., 5., 100.]


    If ``outlier_filter_gap`` is set to 0.5, the outlier (100.) will be removed from the Pearson correlation calculation 
    because the difference between the correlation with and without the outlier is greater than 0.5.


Normalize
---------

- **normalize**: bool, default=False
    If you have multiple samples and this parameter is set to True, the data of the samples are divided by the 
    slope of the linear regression between the data of the first sample and the data of each other samples.
    
    Example:

    .. code::

        sample_1 = [1, 2, 3, 4, 5]
        sample_2 = [2, 4, 6, 8, 10]


    If normalize is set to True, the data of sample_2 will be divided by 2.

