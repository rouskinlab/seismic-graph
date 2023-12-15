Input format
------------

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