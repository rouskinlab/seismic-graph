from . import manipulator, util, plotter
import pandas as pd
import numpy as np
from .util.dump import sort_dict, flatten_json, add_min_cov_field, remove_leading_pound, all_pos
import plotly.graph_objects as go
from custom_inherit import doc_inherit
from .util.docstring import style_child_takes_over_parent
import os
from .util.misc import save_plot, extract_args
import inspect
import json
import tqdm
from .util.normalization import LinFitTable
from sklearn.metrics import roc_auc_score

class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.

    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'])
    """

    attr_list = ["name", "samples"]

    @classmethod
    def from_seismic(cls, data, min_cov=0, filter_by="sample"):
        pass

    @classmethod
    def verbose_init(cls, data, min_cov=0, filter_by="sample"):
        pass
    
    @property
    def df(self):
        return self._df
    
    @df.setter
    def df(self, new_df):
        self._df = new_df
        self.table = LinFitTable(new_df)

    def __init__(self, data=None, min_cov=0, filter_by="sample") -> None:
        """Creates a Study object.

        Args:
            data (dict or list[dict] or pandas.DataFrame, optional): Data to use. Can be a dictionary or list of dictionaries containing SEISMIC-output jsons, or directly a pandas dataframe. Defaults to None.
            min_cov (int, optional): Minimum number of base coverage for a row to be filtered-in. Defaults to 0.
            filter_by (str, optional): Filter rows by sample or study. When filtered by study, if a row passes the filter, rows with the same 'reference', 'section' and 'cluster' fields for all other samples have a sufficient base coverage. Defaults to 'sample'.

        Example:
            >>> study = Study(data = {'sample':'mysample',{'reference1': {'section1': {'cluster1': {'sub_N': [100], 'cov': [1000]}}}}},
                              min_cov=1000,
                              filter_by='sample')
            >>> study = Study(data = (
                                        {'sample':'mysample',{'reference1': {'section1': {'cluster1': {'sub_N': [100], 'cov': [1000]}}}}},
                                        {'sample':'mysample2',{'reference1': {'section1': {'cluster1': {'sub_N': [99], 'cov': [1000]}}}}},
                                        ),
                              min_cov=1000,
                              filter_by='sample')"""
        if data is not None:
            df = pd.DataFrame()

            # If data is a list of json, concatenate them into a single dataframe
            if type(data) is not pd.DataFrame:
                # if data isn't iterable, make it a list
                if not hasattr(data, "__iter__") or isinstance(data, dict):
                    data = [data]

                for sample in data:
                    df = pd.concat(
                        [df, pd.DataFrame(flatten_json(sort_dict(sample)))], axis=0
                    )
            # df.to_csv('/Users/casper/Documents/output/initial_df.csv', index=False)
            df = all_pos(df)
            df = remove_leading_pound(df)
            if "min_cov" not in df.columns:
                df = add_min_cov_field(df)
            # df.to_csv('/Users/casper/Documents/output/processed_df.csv', index=False)
            # Use the dataframe (loaded or created from json)
            self.set_df(df, min_cov=min_cov, filter_by=filter_by)

        else:
            self.df = None

    def set_df(self, df, min_cov=0, filter_by="sample"):
        self.df = df.reset_index(drop=True)

        if min_cov > 0:
            self.df = self.df[self.df["min_cov"] >= min_cov]

        if filter_by == "study":
            self.filter_by_study()

        for attr in ["sample", "reference"]:
            self.df[attr] = self.df[attr].astype(str)

        for attr in ["section", "cluster"]:
            if attr not in self.df.columns:
                self.df[attr] = 0

        if hasattr(self.df, "deltaG"):
            self.df["deltaG"] = self.df["deltaG"].apply(
                lambda x: 0.0 if x == "void" else float(x)
            )

        # convert every cell that's a list into a numpy array
        for attr in self.df.columns:
            if not self.df[attr].empty:
                if isinstance(self.df[attr].iloc[0], list):
                    self.df[attr] = self.df[attr].apply(lambda x: np.array(x))


    def filter_by_study(self):
        df = self.df.groupby(["reference", "section", "cluster"]).filter(
            lambda x: len(self.df["sample"].unique()) == len(x["sample"].unique())
        )
        self.df = df

    def get_samples(self):
        return self.df["sample"].unique()

    def get_references(self, sample: str):
        return self.df[self.df["sample"] == sample]["reference"].unique()

    def get_sections(self, sample: str, reference: str):
        return self.df[
            (self.df["sample"] == sample) & (self.df["reference"] == reference)
        ]["section"].unique()

    def get_clusters(self, sample: str, reference: str, section: str):
        return self.df[
            (self.df["sample"] == sample)
            & (self.df["reference"] == reference)
            & (self.df["section"] == section)
        ]["cluster"].unique()

    def wrap_to_plotter(self, func, loc, kwargs):
        kwargs = {
            **{k: v for k, v in loc.items() if not k in ["self", "args", "kwargs"]},
            **kwargs,
        }

        """Wrapper for the plot functions."""
        return func(
            manipulator.get_df(
                self.df,
                **{
                    k: v
                    for k, v in kwargs.items()
                    if k in list(self.df.columns) + extract_args(manipulator.get_df)
                }
            ),
            **{k: v for k, v in kwargs.items() if k in extract_args(func)}
        )

    def plot_info(name, display_name):
        def decorator(func):
            func.plot_info = {
                "value": name,
                "label": display_name,
                "description": func.__doc__,
            }
            return func

        return decorator

    def get_plots_list(self):
        methods_info = []
        for name, func in inspect.getmembers(self, predicate=inspect.ismethod):
            if hasattr(func, "plot_info"):
                methods_info.append(func.plot_info)
        # return json.dumps(methods_info, indent=4)
        return methods_info

    def default_arguments_per_base(self):
        """Default arguments for the plot functions.

        Args:
            base_index (list, int, str, optional): Filter per-base attributes (sub_rate, sequence, etc) by base index, using 1-indexing. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Gives a  Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (sub_rate, sequence, etc) by base type. Defaults to ``['A','C','G','T']``.
            base_pairing (bool, optional): Filter per-base attributes (sub_rate, sequence, etc) by expected base pairing. True will keep only base pairs, False will keep only non-base pairs. Defaults to None.
            **kwargs: Additional arguments to pass to filter rows by. Ex: ``flank='flank_1'`` will keep only rows with ``flank==flank_1``.

        Returns:
            dict: {``'fig'``: a plotly figure, ``'data'``: a pandas dataframe}

        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    def default_arguments_single_row(self):
        """Default arguments for the mutiple rows plot functions.

        Args:
            sample (str, optional): Selects this sample. Defaults to None.
            reference (str, optional): Selects this reference. Defaults to None.
            section (str, optional): Selects this section. Defaults to ``full``.
            cluster (str, optional): Selects this cluster. Defaults to ``pop_avg``.

        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    def default_arguments_multi_rows(self):
        """Default arguments for the single row plot functions.

        Args:
            sample (list, str, optional): Filter rows by sample (a list of samples or just a sample). Defaults to None.
            reference (list, str, optional): Filter rows by reference (a list of references or just a reference). Defaults to None.
            section (list, str, optional): Filter rows by section (a list of sections or just a section). Defaults to None.
            cluster (list, str, optional): Filter rows by cluster (a list of clusters or just a cluster). Defaults to None.
            normalize (bool, optional): Fit one sample to the other to normalize the mutation fractions. Defaults to False.
        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def get_df(self, **kwargs):
        """Filter the dataframe by the given arguments."""
        return manipulator.get_df(self.df, **kwargs)

    ############################################################################################################
    # Plot functions                                                                                           #
    ############################################################################################################

    @plot_info("mutation_fraction", "Mutation Fraction")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_fraction(
        self, sample, reference, section=None, cluster=None, **kwargs
    ) -> dict:
        """Plot the mutation rates as histograms.

        Args:
            show_ci(bool, optional): Show confidence intervals. Defaults to True.

        """
        return self.wrap_to_plotter(plotter.mutation_fraction, locals(), kwargs)

    @plot_info("mutation_fraction_identity", "Mutation Fraction Identity")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_fraction_identity(
        self, sample, reference, section=None, cluster=None, **kwargs
    ) -> dict:
        """Plot the mutation rates as histograms.

        Args:
            show_ci(bool, optional): Show confidence intervals. Defaults to True.

        """
        return self.wrap_to_plotter(
            plotter.mutation_fraction_identity, locals(), kwargs
        )


    @plot_info(
        "experimental_variable_across_samples", "Experimental Variable Across Samples"
    )
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def experimental_variable_across_samples(
        self, experimental_variable, reference, section, **kwargs
    ) -> dict:
        """Plot a given experimental variable vs Mutation fraction across samples for a given reference and section.

        Args:
            experimental_variable (str): Name of the experimental variable to plot.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].

        """
        index_selected = True
        kwargs['table'] = self.table
        return self.wrap_to_plotter(
            plotter.experimental_variable_across_samples, locals(), kwargs
        )

    # @save_plot
    # @doc_inherit(save_plot, style=style_child_takes_over_parent)
    # @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    # def auc(self, **kwargs)->dict:
    #     """Plot the AUC for each mutation profile of the selected data.

    #     """
    #     unique_id = True
    #     return self.wrap_to_plotter(
    #         plotter.auc,
    #         locals(),
    #         kwargs
    #     )

    @plot_info(
        "num_aligned_reads_per_reference_frequency_distribution",
        "# Aligned Reads / Reference as Freq. Dist.",
    )
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def num_aligned_reads_per_reference_frequency_distribution(
        self, section=None, **kwargs
    ) -> dict:
        """Plot the number of aligned reads per reference as a frequency distribution. x axis is the number of aligned reads per reference, y axis is the count of reference that have this number of aligned reads."""
        return self.wrap_to_plotter(
            plotter.num_aligned_reads_per_reference_frequency_distribution,
            locals(),
            kwargs,
        )

    @plot_info("mutation_fraction_delta", "Mutation Fraction Delta")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    def mutation_fraction_delta(self, **kwargs) -> dict:
        """Plot the Mutation fraction difference between two mutation profiles.
        Returns:
            dict: {'fig': a plotly figure, 'data': a pandas dataframe}

        """
        kwargs['table'] = self.table
        df = manipulator.get_df(
            self.df,
            **{
                k: v
                for k, v in kwargs.items()
                if k
                in list(self.df.columns) + list(manipulator.get_df.__code__.co_varnames)
            }
        )
        assert len(df) > 1, "None or one row found for the mutation profiles."
        assert len(df) <= 2, "More than two row found for the mutation profiles."
        return plotter.mutation_fraction_delta(
            df.reset_index(drop=True),
            **{
                k: v
                for k, v in kwargs.items()
                if k in plotter.mutation_fraction_delta.__code__.co_varnames
            }
        )

    @plot_info("mutations_per_read_per_sample", "Mutations per Read per Sample")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def mutations_per_read_per_sample(self, sample, section=None, **kwargs) -> dict:
        """Plot the number of mutations per read per sample as an histogram."""
        return self.wrap_to_plotter(
            plotter.mutations_per_read_per_sample, locals(), kwargs
        )

    @plot_info("base_coverage", "Base Coverage")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def base_coverage(self, sample, reference, section=None, cluster=None, **kwargs):
        """Plot the base coverage of a single row of your dataframe."""
        return self.wrap_to_plotter(plotter.base_coverage, locals(), kwargs)

    @plot_info("mutation_per_read_per_reference", "Mutation per Read per Reference")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_per_read_per_reference(
        self, sample, reference, section=None, cluster=None, **kwargs
    ) -> dict:
        """Plot the number of mutations per read per reference as an histogram."""
        return self.wrap_to_plotter(
            plotter.mutation_per_read_per_reference, locals(), kwargs
        )

    @plot_info("compare_mutation_profiles", "Compare Mutation Profiles")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def compare_mutation_profiles(self, max_plots=100, max_axis=None, **kwargs):
        """Plot the mutation fraction of multiple mutation profiles.

        Args:
            max_plots: maximum number of plots to show.
            max_axis: maximum value of the x and y axis. If None, the maximum value of the data will be used if above 0.15, otherwise 0.15.
        """
        kwargs["unique_id"] = True
        kwargs['table'] = self.table

        return self.wrap_to_plotter(plotter.compare_mutation_profiles, locals(), kwargs)

    @plot_info(
        "correlation_by_refs_between_samples", "Correlation by Refs. Between Samples"
    )
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def correlation_by_refs_between_samples(self, **kwargs):
        """Plot the correlation between mutation profiles of multiple samples, for each reference.

        Args:
            sample (list, str, optional): Filter rows by sample (use exactly two samples). Defaults to None.
        """
        kwargs['table'] = self.table

        return self.wrap_to_plotter(
            plotter.correlation_by_refs_between_samples, locals(), kwargs
        )
        
    @plot_info('one_pager', 'One Pager')
    # @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def one_pager(self, **kwargs):
        """Plot a one pager of the sample/reference/section/cluster data.
        
        Args:
            xrange (list, float, optional): Range of the x axis for the mutation fraction plot. Defaults to [0, 0.15].
            plot_height_cov (int, optional): Height of the base coverage plot. Defaults to 250.
            plot_height_count (int, optional): Height of the mutation count plot. Defaults to 200.
            plot_height_bar (int, optional): Height of the mutation fraction bar plot. Defaults to 225.
            plot_width_first_col (int, optional): Width of the first column. Defaults to 600.
            margin (dict, optional): Plot margin. Defaults to dict(l=0, r=0, t=25, b=10).
            plot_width_bar (int, optional): Width of the mutation fraction bar plot. Defaults to 900.
        """
        return self.wrap_to_plotter(plotter.one_pager, locals(), kwargs)

    # @plot_info("dist_of_seq_lengths", "Distribution of Sequence Lengths")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def dist_of_seq_lengths(self, **kwargs) -> dict:
        """Plot a histogram of sequence lengths across all rows in the dataframe.
        """
        return self.wrap_to_plotter(
            plotter.dist_of_seq_lengths, locals(), kwargs
        )
    
    # @plot_info("percent_masked_histogram", "Percent Masked Histogram")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def percent_masked_histogram(self, **kwargs) -> dict:
        """Plot a histogram of the percentage of bases that are masked across all rows in the dataframe.

        """
        return self.wrap_to_plotter(
            plotter.percent_masked_histogram, locals(), kwargs
        )
    
    # @plot_info("f1_violin_by_category", "F1 Violin Plot by Category")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def f1_violin_by_category(self, **kwargs) -> dict:
        """
        Generate violin plots grouped by the specified category (default: 'family').
        If 'F1_no_dms' exists in the DataFrame, the violin plot is split to
        compare 'F1' (with DMS) vs 'F1_no_dms' (without DMS). Otherwise, a
        single violin is displayed for 'F1'.

        Args:
            category (str, optional): Column name used to group by. Defaults to 'family'.

        Returns:
            dict: {
                'fig': Plotly Figure,
                'data': DataFrame filtered of NaNs, possibly containing 'F1_no_dms' if present
            }
        """
        return self.wrap_to_plotter(
            plotter.f1_violin_by_category, locals(), kwargs
        )
    

    # @plot_info("f1_violin_by_category_symmetrical", "F1 Violin Plot by Category, Symmetrical")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def f1_violin_by_category_symmetrical(self, **kwargs) -> dict:
        """
        Generate symmetrical violin plots grouped by the specified category (default: 'family').

        Args:
            category (str, optional): Column name used to group by. Defaults to 'family'.

        Returns:
            dict: {
                'fig': Plotly Figure,
                'data': DataFrame filtered of NaNs
            }
        """
        return self.wrap_to_plotter(
            plotter.f1_violin_by_category_symmetrical, locals(), kwargs
        )


    def add_sections_from_library(self, library):
        """
        Add sections to the study df from a library

        Args:
            library: path to a csv file containing the library

        Returns:
            df (pandas.DataFrame): a dataframe with the sections added

        """
        lib = pd.read_csv(library)

        stack = []

        for (ss, se), rows in tqdm.tqdm(
            lib.groupby(["section_start", "section_end"]),
            total=len(lib.groupby(["section_start", "section_end"])),
        ):
            subdf = (
                self.get_df(
                    reference=rows["reference"].values,
                    section="full",
                    base_index=list(range(ss, 1 + se)),
                )
                .reset_index(drop=True)
                .copy()
            )
            if len(rows) == 0:
                continue

            for attr in rows.keys():
                for ref, g in subdf.groupby("reference"):
                    subdf.loc[g.index, attr] = rows.loc[
                        rows["reference"] == ref, attr
                    ].values[0]

            stack.append(subdf)

        df = pd.concat(stack, ignore_index=True)
        df = pd.concat(
            [df, self.get_df(section="full")], ignore_index=True
        ).reset_index(drop=True)
        df.drop_duplicates(
            subset=["sample", "reference", "section", "cluster"], inplace=True
        )
        df.dropna(
            subset=[
                "section_start",
                "section_end",
                "sample",
                "reference",
                "section",
                "cluster",
            ],
            inplace=True,
        )
        df.reset_index(drop=True, inplace=True)

        self.df = df

        return df

    def get_exp_env(self, sample):
        return self.get_df(sample=sample)["exp_env"].values[0]


    # @plot_info("compare_mutation_profiles_2", "Compare Mutation Profiles 2")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def compare_mutation_profiles_2(self, max_plots=100, max_axis=None, **kwargs):
        """Plot the mutation fraction of multiple mutation profiles.

        Args:
            max_plots: maximum number of plots to show.
            max_axis: maximum value of the x and y axis. If None, the maximum value of the data will be used if above 0.15, otherwise 0.15.
        """
        kwargs["unique_id"] = True
        kwargs['table'] = self.table

        return self.wrap_to_plotter(plotter.compare_mutation_profiles_2, locals(), kwargs)
    
    @plot_info("pearson_correlation_histogram", "Pearson Correlation Histogram")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def pearson_correlation_histogram(self, **kwargs):
        """Plot a histogram showing the distribution of Pearson R² correlations across references.

        Args:
            sample (list, str, optional): Filter rows by sample (use exactly two samples). Defaults to None.
        """
        return self.wrap_to_plotter(
            plotter.pearson_correlation_histogram, locals(), kwargs
        )
    
    def identify_high_mutation_sequences(self, **kwargs):
        """
        Identify sequences in the untreated data (NO DMS) that have high mutations.

        Args:
            signal_threshold_low (float): The lower threshold for mutation fraction.
            signal_threshold_high (float): The higher threshold for mutation fraction.

        Returns:
            pd.DataFrame: The dataframe with additional columns indicating acceptable sequences.
        """
        return self.wrap_to_plotter(
            plotter.identify_high_mutation_sequences, locals(), kwargs
        )
    
    @plot_info("max_sub_rate_histogram", "Histogram of Maximum Sub Rate per Reference")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def max_sub_rate_histogram(self, **kwargs):
        """
        Generate a histogram of maximum sub rates. For each row in the data, it finds the maximum mutation rate and generates a histogram.

        Returns:
            dict: Contains the Plotly figure "fig" and histogram data "histogram_data".
        """
        return self.wrap_to_plotter(
            plotter.max_sub_rate_histogram, locals(), kwargs
        )
    
    def sequence_length_vs_aligned_reads(self, **kwargs):
        """
        Generate a scatter plot with sequence length on the x-axis and number of aligned reads on the y-axis.
        Adds a line of best fit and indicates the R² value.

        Returns:
            dict: Contains the Plotly figure and data used for the plot.
        """
        return self.wrap_to_plotter(
            plotter.sequence_length_vs_aligned_reads, locals(), kwargs
        )
    
    def count_rows_with_all_nan_cov(self, **kwargs) -> dict:
        """Count the number of rows in the dataframe where 'cov' is all NaN.

        Returns:
            dict: {'count': integer count of rows where 'cov' is all NaN,
                   'rows': DataFrame with these rows}
        """
        return self.wrap_to_plotter(plotter.count_rows_with_all_nan_cov, locals(), kwargs)
    
    @plot_info("percent_masked_histogram", "Percent Masked Histogram")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def percent_masked_histogram(self, **kwargs) -> dict:
        """Plot a histogram showing the distribution of percentage of A's and C's bases masked.

        Args:
            potentially_reactive_bases (list of str, optional): List of bases considered potentially reactive. Defaults to ['A', 'C'].

        Returns:
            dict: {'fig': a Plotly figure, 'data': data}
        """
        return self.wrap_to_plotter(
            plotter.percent_masked_histogram, locals(), kwargs
        )
    

    def ensure_auroc_column(self):
        """
        If the 'auroc' column does not exist in self.df, compute it for every row
        in the entire DataFrame. Modify self.df in place.
        """
        if 'auroc' not in self.df.columns:
            
            def compute_auroc(structure, sub_rate):
                # Convert structure to paired/unpaired status
                structure_array = np.array([0 if c in '()' else 1 for c in structure])
                # Convert sub_rate to numpy array
                sub_rate = np.array(sub_rate, dtype=np.float64)
                # Check lengths
                if len(structure_array) != len(sub_rate):
                    return np.nan
                # Remove NaNs
                valid_positions = ~np.isnan(sub_rate)
                sub_rate = sub_rate[valid_positions]
                structure_array = structure_array[valid_positions]
                # Need at least two classes
                unique_classes = np.unique(structure_array)
                if len(unique_classes) < 2:
                    return np.nan
                # Compute AUROC
                try:
                    return roc_auc_score(structure_array, sub_rate)
                except ValueError:
                    return np.nan

            # Compute and add the 'auroc' column in-place
            self.df['auroc'] = self.df.apply(
                lambda row: compute_auroc(row['structure'], row['sub_rate']),
                axis=1
            )

    @plot_info("auroc_histogram", "AUROC Histogram")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def auroc_histogram(self, **kwargs):
        """
        Create a histogram of AUROC scores.

        - Ensures the 'auroc' column is present on the entire self.df
        - Then uses wrap_to_plotter, which filters the DataFrame, and calls plotter.auroc_histogram

        Returns:
            dict: {'fig': a Plotly figure, 'data': DataFrame with 'auroc' column added}
        """
        self.ensure_auroc_column()

        return self.wrap_to_plotter(
            plotter.auroc_histogram, locals(), kwargs
        )
    
    # @plot_info("refs_per_category", "Number of References per Category")
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def refs_per_category(self, category="sample", **kwargs) -> dict:
        """
        Generate a bar plot showing the number of references per category.

        Args:
            category (str): The column name to group by. Default is "sample".

        Returns:
            dict: {'fig': a Plotly figure, 'data': the DataFrame used for plotting}
        """
        return self.wrap_to_plotter(
            plotter.refs_per_category, locals(), kwargs
        )