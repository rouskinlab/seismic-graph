.. _SEISMICgraph: https://seismicrna.org/
.. |SEISMICgraph| replace:: **SEISMICgraph** 


===================================
SEISMICgraph: Web App Documentation
===================================

Overview
--------

|SEISMICgraph|_ is a webapp that enables researchers to utilize the seismic-graph Python package without learning the requisite coding skills.

SEISMICgraph accepts mutational profile data from SEISMIC-RNA, ShapeMapper2, or RNA Framework, and generates plots and analysis of that data.

Workflow for SEISMICgraph
-------------------------

The workflow for SEISMICgraph is as follows:

1. **Upload a dataset**, or use the sample dataset.
2. **Make Data Selections**:

   * Filter the data by making selections for which parts of the dataset to use.
   * As you filter the dataset, you will see that the number of **Rows Selected** updates with your adjustments. A row is an entry in the dataset, and it can be uniquely identified by its **Sample**, **Reference**, **Section**, and **Cluster**.
3. **Plot Selections**:

   * Choose what type of plot to generate from the pre-existing plot types.
   * Some plots require a specific number of rows to be selected. For example, **Mutation Fraction** requires a single entry to be selected. To select a single row, choose values for **Sample**, **Reference**, **Section**, and **Cluster** to uniquely identify the row of interest.
   * Some plots have options that are specific to them.
4. **Plot**:

   * The plot will generate. If the plot could not be generated, a message will be displayed indicating why it could not be generated.

Let's go through these steps using the sample dataset as an example.



SEISMICgraph: Step by Step walkthrough
--------------------------------------



Navigate to Website
^^^^^^^^^^^^^^^^^^^

SEISMICgraph is located at `seismicrna.org <https://seismicrna.org/>`_.

Upload Dataset
^^^^^^^^^^^^^^

In this example, we will use the **Sample Dataset**. In the light blue panel, choose **Use Sample Dataset**. The sample dataset will be loaded into the app. This usually takes about 15 seconds.

Once the dataset has been loaded, we will see a few changes to the interface: under **Data Selections** we now see **Rows Selected: 1401**, and **Section** and **Cluster** have both been populated with values:

- **Section:** 20-151
- **Cluster:** average

Filter Dataset using Data Selections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the filters **Sample**, **Reference**, **Section**, and **Cluster**, if there is only one value for that filter given the other filters, then that value will be automatically selected. For example, this dataset does not have clustering. The only **Cluster** is **average**, indicating the population average mutation rate. Since every selected row has the same value for **Cluster**, that value is displayed as the selection for **Cluster**, to make it clear that all selected rows are using population average.

For this example, we will generate a plot **Mutation Fraction**. **Mutation Fraction** displays the proportion of reads of a base that were mutations for every base in the reference sequence. It requires a single row to be selected, so we will use **Sample** and **Reference** to uniquely identify the row of interest.

By selecting the **Reference** of interest and clicking outside of the dropdown menu, the selection is sent, and when the loading indicator goes away, the **Rows Selected** updates to reflect the new selection. Once a **Sample** and **Reference** are selected, we can see that the dataset has been filtered down to one row. 

**Note:** Datasets with multiple **Sections** and/or **Clusters** will require selections for those filters as well to uniquely identify a single row.

Make Plot Selections
^^^^^^^^^^^^^^^^^^^^

- Choose **Mutation Fraction**.

Plot
^^^^

Click **Plot**. The plot will be generated in around 15 seconds. More computationally intensive plots like **Number of Aligned Reads per Reference as a Frequency Distribution** will take longer, especially on larger datasets.