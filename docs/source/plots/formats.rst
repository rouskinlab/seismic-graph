
Input formats
-------------

SEISMIC-graph is compatible with the SEISMIC-output format, the shape-mapper format, and the RNA framework format. These formats are defined below.

SEISMIC-OUTPUT format
+++++++++++++++++++++

The input format is a JSON file. The data is stored on a 4 levels hierarchical structure: sample, reference, section, and cluster. 
Each level contains metadata and instances of the next level. Metadata names start with a "#" sign. 

Levels
******

- **sample level**: contains metadata about the whole file (ex: temperature, salt concentration, sample name, etc) and a list of its references. 
- **reference level**: A reference is the name of a sequence (for example, in a fasta file). The reference level contains metadata about one specific whole sequence (ex: number of aligned reads, RNA family, etc), and a list of its sections. 
- **section level**: A section is a subset of a sequence. For example, a specific region of interest, or simply the whole sequence. The section level contains metadata about one specific region (ex: name, sequence, structure and free energy prediction using RNAstructure, etc) and a list of its clusters. 
- **cluster level**: A cluster is an aggregation of certain reads based on their mutation profile (ref DREEM, DRACO, and DANCEmap). The ``average`` cluster contains all reads. The reactivity signal, the base coverage and other signals are contained in the cluster. 

Cluster level data
******************

- ``positions``: An array giving the one indexed positions.
- ``cov``: An array giving the number of valid reads per position, non-including deleted bases. 
- ``del``: An array giving the number of deletions per position.
- ``info``: An array giving the number of valid reads per position, including deleted bases.
- ``ins``: An array giving the number of insertions per position.
- ``sub_{base}``: An array giving the number of reads that mutated to this base for each position. N includes all bases.
- ``sub_hist``: A histogram of the number of mutations per read.
- ``sub_rate``: An array giving the substitution rate for each position. This is computed by dividing sub_N by info. 

Example:
********

.. code-block:: json

  {
    "#sample": "this sample's name",
    "#temperature_K": 310, 
    "a reference" : { 
        "#num_aligned": 1998,
        "a section": { 
            "#sequence": "TTAAACCGGCCAACATCAA",
            "average": { 
                "name": "average",
                "positions": ["one indexed positions"],
                "cov": ["some values"],
                "del": ["some values"],
                "info": ["some values"],
                "ins": ["some values"],
                "sub_A": ["some values"],
                "sub_C": ["some values"],
                "sub_G": ["some values"],
                "sub_T": ["some values"],
                "sub_N": ["some values"],
                "sub_hist": ["some values"],
                "sub_rate": ["some values"]
              },
            "some other cluster": {   },
          },
        "some other section": {   },
      },
    "some other reference": {   },
  }



SHAPE-MAPPER format
+++++++++++++++++++

A definition of the shape-mapper format can be found `here <https://github.com/Weeks-UNC/shapemapper2/blob/master/docs/file_formats.md#name_rna_profiletxt>`__.


Format example
**************

.. code-block:: text
  
  # name: sample_name/reference_name.txt

  Nucleotide	Sequence	Modified_mutations	Modified_read_depth	Modified_effective_depth	Modified_rate	Modified_off_target_mapped_depth	Modified_low_mapq_mapped_depth	Modified_primer_pair_1_mapped_depth	Untreated_mutations	Untreated_read_depth	Untreated_effective_depth	Untreated_rate	Untreated_off_target_mapped_depth	Untreated_low_mapq_mapped_depth	Untreated_primer_pair_1_mapped_depth	Denatured_mutations	Denatured_read_depth	Denatured_effective_depth	Denatured_rate	Denatured_off_target_mapped_depth	Denatured_low_mapq_mapped_depth	Denatured_primer_pair_1_mapped_depth	Reactivity_profile	Std_err	HQ_profile	HQ_stderr	Norm_profile	Norm_stderr
  1	U	3	5835	4214	0.000712	127	0	5835	3	6037	4225	0.000710	76	0	6037	6	5943	4189	0.001432	166	0	5943	0.001294	0.405299	0.001294	0.405299	0.000579	0.181487
  2	C	13	5852	4264	0.003049	128	0	5852	16	6052	4415	0.003624	76	0	6052	29	5957	4266	0.006798	166	0	5957	-0.084618	0.182980	-0.084618	0.182980	-0.037891	0.081936
  3	G	86	5874	4756	0.018082	129	0	5874	95	6087	4933	0.019258	76	0	6087	102	5996	4502	0.022657	166	0	5996	-0.051889	0.122631	-0.051889	0.122631	-0.023235	0.054913
  4	G	39	5926	5138	0.007591	129	0	5926	53	6168	5417	0.009784	77	0	6168	199	6048	5041	0.039476	167	0	6048	-0.055565	0.046071	-0.055565	0.046071	-0.024881	0.020630
  5	G	27	5970	5378	0.005020	129	0	5970	50	6215	5784	0.008645	77	0	6215	55	6087	5297	0.010383	168	0	6087	-0.349032	0.157278	-0.349032	0.157278	-0.156292	0.070427


Conversion to SEISMIC format
****************************

Mapping:

The shape-mapper files are converted map 3 samples, ``modified``, ``untreated`` and ``denatured``. 

The columns are converted to a SEISMIC format as follows:

- ``Nucleotide``: ``positions``
- ``Sequence``: ``sequence``

Then for each sample:

- ``{sample}_read_depth``: ``cov`` 
- ``{sample}_effective_depth``: ``info`` 
- ``{sample}_mutations``: ``sub_N`` 
- ``{sample}_rate``: ``sub_rate`` 

Naming:

The sample are named ``{sample_name}_{sample_type}`` where ``sample_name`` is the folder name and ``sample_type`` is ``modified``, ``untreated`` or ``denatured``.
The reference is the file name. 
The section is ``full`` by default. 
The cluster is ``pop_avg`` by default. 


Example:

Using the shape-mapper file above, the SEISMIC format would be:

.. code-block:: json

  # sample: test_data
  {
  "#sample":"test_data",
  "example2":{
    "full":{
      "#sequence":"TCGGG",
      "#positions":[1,2,3,4,5],
      "#num_aligned":5970,
      "average":{
        "cov":[5835,5852,5874,5926,5970],
        "info":[4214,4264,4756,5138,5378],
        "sub_N":[3,13,86,39,27],
        "sub_rate":[0.000712,0.003049,0.018082,0.007591,0.00502]
        }
      }
    }

    # sample: test_data-untreated
    {
    "#sample":"test_data-untreated",
    "example2":{
      "full":{
        "#sequence":"TCGGG",
        "#positions":[1,2,3,4,5],
        "#num_aligned":6168,
        "average":{
          "cov":[6037,6052,6087,6168,6215],
          "info":[4225,4415,4933,5417,5784],
          "sub_N":[3,16,95,53,50],
          "sub_rate":[0.00071,0.003624,0.019258,0.009784,0.008645]
          }
        }
      }
    }

    # sample: test_data-denatured
    {
    "#sample":"test_data-denatured",
    "example2":{
      "full":{
        "#sequence":"TCGGG",
        "#positions":[1,2,3,4,5],
        "#num_aligned":6087,
        "average":{
          "cov":[5943,5957,5996,6048,6087],
          "info":[4189,4266,4502,5041,5297],
          "sub_N":[6,29,102,199,55],
          "sub_rate":[0.001432,0.006798,0.022657,0.039476,0.010383]
          }
        }
      }
    }


RNA framework format
++++++++++++++++++++

A definition of the RNA framework format can be found `here <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6144828/table/tbl1/?report=objectonly>`__.

Example:

**Folder structure**

.. code-block:: text

  sample_name/
    sample_name_R1_001_count_data.txt
    reference_1_name.xml
    reference_2_name.xml
    ...

**.txt format**

.. code-block:: text

  # name: sample_name_R1_001_count_data.txt

  HFE_amp_1
  AGATGATACAAAAAAACATGACTACATGATAAGTACAAGAGGAGACAGACGACAGTGTCCACAGCACCCGTTTCAGCACAGTTGGAGGAGAGGGGATAAGATTTATTGATGAAATTTGTGATTTGCATCGTGGTACAGAAAAGTTATGTGAATATAAAAGTGTAGAACAATGTCTTCCGATTTCGACAGGTTAGAAGATGGGGAAGAGCAGGCATTTTGGAGAAGGCGAGGGCGACG
  138,63,72,38,42,51,17,51,22,9,5,4,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  138,201,273,311,353,404,421,472,494,503,508,512,512,512,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,513,514,514,514,514,514,514,514,514,514,514,514,514,514,514,514,514,514,514,515,516,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,517,508,443,379,272,241,186,142,103,56,15,15,9,7,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

  HFE_amp_2
  GGCAAGCGAAAGATTTTGAAACTTTCCGAGAAGGGGGAACAGAGGGGTAAGGGGCTCCGGTTTAGACAGAGGAACGTGACAAAGAGACAGAAGTTGGGGCGAGCAGGCTTTCAGGAAGGATTCTTGATGAGGGGGAGGGGATAAACAGGGAGGAGAGAGAGGGGAATCGATAGCGGCGGGGCAGAAAGAAGAATAGAAGGGGGCCGCGAGGAGGGGAGAGTCGAAGGATGAGTGAAGGAGAAGGAAGGG
  440,66,62,33,11,9,7,7,7,14,3,4,2,3,15,18,11,1,6,5,4,0,2,2,1,0,1,2,2,0,0,0,4,4,0,0,0,0,3,0,3,2,0,1,1,1,2,1,0,1,0,0,1,0,3,0,0,0,2,5,0,0,3,0,1,1,0,0,0,0,2,3,1,2,0,0,18,4,2,0,1,0,2,1,0,0,2,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,3,0,0,0,1,0,2,0,1,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  440,506,568,601,612,621,628,635,642,656,659,663,665,668,683,701,712,713,719,724,728,728,730,732,733,733,734,736,738,738,738,738,742,746,746,746,746,746,749,749,752,754,754,755,756,757,759,760,760,761,761,761,762,762,765,765,765,765,767,772,772,772,775,775,776,777,777,777,777,777,779,782,783,785,785,785,803,807,809,809,810,810,812,813,813,813,815,815,815,816,816,816,816,816,817,817,818,818,818,819,819,819,820,820,820,821,824,824,824,824,825,825,827,827,828,829,830,830,830,830,830,830,831,831,831,831,831,831,831,831,830,830,830,830,830,830,829,829,829,829,829,829,829,829,829,828,827,819,600,362,340,253,224,219,207,205,192,186,176,172,168,166,164,157,147,129,118,117,112,103,103,103,101,98,98,98,97,95,93,93,88,87,87,84,84,84,84,82,79,79,77,76,74,74,74,72,71,71,70,70,70,70,69,69,66,66,60,60,59,59,59,59,56,56,54,54,54,54,54,54,48,45,44,44,44,43,41,25,22,21,21,20,18,18,17,16,15,15,15,15,15,15,15,15,13,13,11,0,0

**.xml format**

.. code-block:: xml
  
  # name: sample_name/reference_1_name.xml

  <?xml version="1.0" encoding="UTF-8"?>
    <data combined="FALSE" maxmutrate="0.2" norm="90% Winsorizing" offset="1000000000" reactive="AC" remap="0" scoring="Zubradt" tool="rf-norm" win="1000000000">
      <transcript id="HFE_amp_1" length="237">
        <sequence>
          AGATGATACAAAAAAACATGACTACATGATAAGTACAAGAGGAGACAGACGACAGTGTCC
          ACAGCACCCGTTTCAGCACAGTTGGAGGAGAGGGGATAAGATTTATTGATGAAATTTGTG
          ATTTGCATCGTGGTACAGAAAAGTTATGTGAATATAAAAGTGTAGAACAATGTCTTCCGA
          TTTCGACAGGTTAGAAGATGGGGAAGAGCAGGCATTTTGGAGAAGGCGAGGGCGACG
        </sequence>
        <reactivity>
          NaN,NaN,NaN,NaN,NaN,1.000,NaN,1.000,1.000,1.000,1.000,0.794,0.000,0.000,0.198,0.000,0.000,0.000,NaN,NaN,0.000,0.000,NaN,0.000,0.000,0.000,NaN,NaN,0.000,NaN,0.000,0.000,NaN,NaN,0.000,0.000,0.198,0.000,NaN,0.000,NaN,NaN,0.000,NaN,0.000,0.000,0.000,NaN,0.000,0.000,NaN,0.000,0.000,0.000,NaN,NaN,NaN,NaN,0.000,0.000,
          0.000,0.000,0.000,NaN,0.000,0.000,0.000,0.000,0.000,NaN,NaN,NaN,NaN,0.000,0.000,NaN,0.000,0.000,0.000,0.000,NaN,NaN,NaN,NaN,NaN,0.000,NaN,NaN,0.000,NaN,0.000,NaN,NaN,NaN,NaN,0.000,NaN,0.000,0.000,NaN,0.000,NaN,NaN,NaN,0.000,NaN,NaN,NaN,0.000,NaN,NaN,0.000,0.000,0.000,NaN,NaN,NaN,NaN,NaN,NaN,
          0.000,NaN,NaN,NaN,NaN,0.000,0.000,NaN,0.000,NaN,NaN,NaN,NaN,NaN,0.000,0.000,0.000,NaN,0.000,0.000,0.000,0.000,NaN,NaN,NaN,0.000,NaN,NaN,NaN,NaN,0.000,0.000,NaN,0.000,NaN,0.000,0.000,0.000,0.000,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,
          NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
        </reactivity>
      </transcript>
    </data>
