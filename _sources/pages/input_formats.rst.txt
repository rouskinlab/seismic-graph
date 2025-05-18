=======================
Compatible Data Formats
=======================

Overview
--------

seismic-graph works natively with DMS-MaPseq data processed by SEISMIC-RNA. When SEISMIC-RNA is run with with the :code:`--export` flag, it generates a json file, named <your-sample-name>__webapp.json. This file contains all of the information that seismic-graph needs to generate a variety of graphs.

seismic-graph is also capable of using data from RNAFramework or ShapeMapper2. RNAFramework and ShapeMapper2 use different data formats that include some but not all of the information that is saved in the SEISMIC-RNA export json. seismic-graph is able to read data from RNAFramework and ShapeMapper2, however, due to the limitations of their export formats, seismic-graph is limited in the graphs it can generate from those sources.

seismic-graph works with one data format at a time, so all files uploaded together must be processed by the same tool. For example, if you want to work with multiple SEISMIC-RNA webapp.json files at the same time, you can do so, but you cannot  process output from SEISMIC-RNA and RNAFramework in the same session.


SEISMIC-RNA webapp json
-----------------------

In order to export data from SEISMIC-RNA for use with seismic-graph, just include the :code:`--export` flag in your SEISMIC-RNA command.

Documentation for SEISMIC-RNA is available at `rouskinlab.github.io/seismic-rna <https://rouskinlab.github.io/seismic-rna/#>`_.

RNAFramework
------------


ShapeMapper2
------------