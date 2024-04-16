
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

We write the shape-mapper file above into a file ``my_example_sample/my_example_ref.txt``. The SEISMIC format is:

.. code-block:: text

  # sample: my_example_sample
  {
  "#sample":"my_example_sample",
  "my_example_ref":{
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

    # sample: my_example_sample-untreated
    {
    "#sample":"my_example_sample-untreated",
    "my_example_ref":{
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

    # sample: my_example_sample-denatured
    {
    "#sample":"my_example_sample-denatured",
    "my_example_ref":{
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


Format example
**************

.. code-block:: text

  # name: example_rnaf/sample_name_R1_001_count_data.txt
  HFE_amp_1
  AGATGATACAAAAAAACATGACTAC
  138,63,72,38,42,51,17,51,22,9,5,4,0,0,1,0,0,0,0,0,0,0,0,0,0
  138,201,273,311,353,404,421,472,494,503,508,512,512,512,513,513,513,513,513,513,513,513,513,513,513

  HFE_amp_2
  GGCAAGCGAAAGATTTTGAAACTTT
  440,66,62,33,11,9,7,7,7,14,3,4,2,3,15,18,11,1,6,5,4,0,2,2,1
  440,506,568,601,612,621,628,635,642,656,659,663,665,668,683,701,712,713,719,724,728,728,730,732,733


Conversion to SEISMIC format
****************************

Mapping for the .txt format:

- The first line is the reference name.
- The second line is the sequence.
- The third line is the number of mutations per position.
- The fourth line is the number of reads per position.



.. code-block::

  {
    "#sample": "example_rnaf",
    "HFE_amp_1": {
      "full": {
        "average": {
          "sequence": "AGATGATACAAAAAAACATGACTAC",
          "positions": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
          ],
          "cov": [138,201,273,311,353,404,421,472,494,503,508,512,512,512,513,513,513,513,513,513,513,513,513,513,513
          ],
          "sub_rate": [1,0.31343283582089554,0.26373626373626374,0.12218649517684887,0.11898016997167139,0.12623762376237624,0.040380047505938245,0.10805084745762712,0.044534412955465584,0.017892644135188866,0.00984251968503937,0.0078125,0,0,0.001949317738791423,0,0,0,0,0,0,0,0,0,0
          ]
        }
      }
    },
    "HFE_amp_2": {
      "full": {
        "average": {
          "sequence": "GGCAAGCGAAAGATTTTGAAACTTT",
          "positions": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
          ],
          "cov": [440,506,568,601,612,621,628,635,642,656,659,663,665,668,683,701,712,713,719,724,728,728,730,732,733
          ],
          "sub_rate": [1,0.13043478260869565,0.10915492957746478,0.05490848585690516,0.017973856209150325,0.014492753623188406,0.011146496815286623,0.011023622047244094,0.010903426791277258,0.021341463414634148,0.004552352048558422,0.006033182503770739,0.0030075187969924814,0.004491017964071856,0.021961932650073207,0.025677603423680456,0.01544943820224719,0.001402524544179523,0.008344923504867872,0.006906077348066298,0.005494505494505495,0,0.0027397260273972603,0.00273224043715847,0.001364256480218281
          ]
        }
      }
    }
  }