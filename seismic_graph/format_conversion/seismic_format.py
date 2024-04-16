import os, sys
import pandas as pd
import xmltodict 
import numpy as np

UKN = np.nan
def find_files_and_names(path, name, exts):
    if name is None:
        if len(path.split("/")) > 1:
            name = path.split("/")[-2]
        else:
            name = 'N/A'
    
    files = [path] if os.path.isfile(path) else [f'{path}/{file}' for file in os.listdir(path) if file.split('.')[-1] in exts]
    return files, name

def probe_path_depth(path):
    """Return the depth of the path. 0 if the path is a file, 1 if the path is a directory containing files, 2 if the path is a directory containing directories containing files, etc.
    """
    
    if os.path.isfile(path):
        return 0

    if os.path.isdir(path):
        return 1 + max([probe_path_depth(f'{path}/{file}') for file in os.listdir(path)])


class SEISMICdata:
    def __init__(self, samples=None, **kwargs):
        self.samples = [] if samples is None else samples
        self.kwargs = kwargs
        
    def to_dict(self):
        return {
            **{'#'+k:v for k,v in self.kwargs.items()},
            **{sample.name: sample.to_dict() for sample in self.samples},
        }
        
    def __repr__(self) -> str:
        return '[' + '\n'.join([str(sample) for sample in self.samples]) + ']'
    
    def to_list(self):
        return [sample.to_dict() for sample in self.samples]
    
    
    @classmethod
    def from_shape(cls, path, name = None):
        """Create a SEISMICdata from a SHAPE file. 
        
        Args:
            path (str): Path to the SHAPE file or directory containing SHAPE files.
            name (str): Name of the sample. If None, the name will be the last part of the path.
        """
        sample_types = ['modified', 'untreated', 'denatured']
        if probe_path_depth(path) <= 1:
            return cls([
                Sample.from_shape(path, 
                                  name,
                                  sample_type=sample_type) 
                for sample_type in sample_types
                ])
        return cls([
            Sample.from_shape(
                sample_folder, 
                sample_folder.split('/')[-1],
                sample_type=sample_type) 
            for sample_folder in os.listdir(path)
            for sample_type in sample_types
            ])
    
    @classmethod
    def from_rnaf(cls, path, name = None):
        """Create a SEISMICdata from a RNAF file. 
        
        Args:
            path (str): Path to the directory containing RNAF files. Must contrain a .txt file and one or more .xml files.
            name (str): Name of the sample. If None, the name will be the last part of the path.
        """
        if probe_path_depth(path) <= 1:
            return cls([Sample.from_rnaf(path, name)])
        return cls([Sample.from_rnaf(sample_folder) for sample_folder in os.listdir(path)])
    

class Sample:
    def __init__(self, name, refs=None, **kwargs):
        self.refs = [] if refs is None else refs
        self.name = name
        self.kwargs = kwargs
        
    def to_dict(self):
        return {
            '#sample': self.name,
            **{'#'+k:v for k,v in self.kwargs.items()},
            **{ref.name: ref.to_dict() for ref in self.refs},
        }
        
    def __repr__(self) -> str:
        return f"SEISMICsample({self.name}, {self.kwargs}, \n\t" + '\n\t'.join([str(ref) for ref in self.refs]) + '\n)'
    
        
    @classmethod
    def from_shape(cls, path, prefix, sample_type):
        """Create a SEISMICsample from a SHAPE file. 
        
        Args:
            path (str): Path to the SHAPE file or directory containing SHAPE files.
            prefix (str): Name of the sample. If None, the name will be the last part of the path.
        """
        files, prefix = find_files_and_names(path, prefix, ['shape', 'SHAPE', 'txt'])
        name = prefix + '-' + sample_type if sample_type != 'modified' else prefix
        return cls(name, [Reference.from_shape(file, sample_type=sample_type) for file in files])
    
    @classmethod
    def from_rnaf(cls, path, name = None):
        """Create a SEISMICsample from a RNAF file. 
        
        Args:
            path (str): Path to the directory containing RNAF files. Must contrain a .txt file and one or more .xml files.
            name (str): Name of the sample. If None, the name will be the last part of the path.
        """
        files, _ = find_files_and_names(path, name, ['txt', 'xml'])
        name = name if name is not None else path.split('/')[-1].split('.')[0]
        refs = {}
        # parse the data
        for file in files:
            ext = file.split('.')[-1]
            if ext == 'xml':
                ref_name = file.split('/')[-1].split('.')[0]
                refs[ref_name] = Reference.from_rnaf_xml(file)
            elif ext == 'txt':
                with open(file, 'r') as f:
                    data = f.readlines()
                # add one reference per block of 4 lines
                for i, line in enumerate(data):
                    if i%5 == 0:
                        reference_name = line.strip()
                    if i%5 == 1:
                        sequence = line.strip()
                    if i%5 == 2:
                        mutations = [int(m) for m in line.split(',')]
                    if i%5 == 3:
                        cov = [int(c) for c in line.split(',')]
                        sub_rate = [float(m)/float(c) if c else 0 for m,c in zip(mutations, cov)]
                        refs[reference_name] = Reference(reference_name, 
                                            [Section('full', 
                                                    clusters=[Cluster('average', sequence=sequence, positions=list(range(1, len(sequence) + 1)), cov=cov, sub_rate=sub_rate)])
                                            ])
                    if i%5 == 4:
                        assert line.strip() == '', f"Line {i} is not empty."
            else:
                raise ValueError(f"File {file} does not have a valid extension.")
        return cls(name, list(refs.values()))
       

class Reference:
    def __init__(self, name, sections=None, **kwargs):
        self.name = name
        self.sections = [] if sections is None else sections 
        self.kwargs = kwargs    

    def to_dict(self):
        return {
            **{'#'+k:v.to_dict() for k,v in self.kwargs.items()},
            **{section.name: section.to_dict() for section in self.sections},
        }
        
    def __repr__(self) -> str:
        return f"Reference({self.name}, {self.kwargs}, \n\t\t" + '\n\t\t'.join([str(section) for section in self.sections]) + '\n)'
        
    @classmethod
    def from_shape(cls, file, sample_type):
        assert sample_type in ['modified', 'untreated', 'denatured'], "Sample type must be 'sample_type', 'untreated', or 'denatured'."
        sample_type = sample_type.capitalize()

        assert os.path.isfile(file), f"File {file} does not exist or is not a file."
        data = pd.read_csv(file, sep='\t')
        name = file.split('/')[-1].split('.')[0]
        
        # add lines where nucleotide is missing
        temp = pd.DataFrame({'Nucleotide': list(range(1, len(data)+1))})
        data = pd.merge(temp, data, on='Nucleotide', how='left') 
        
        # Pick columns that are relevant to the sample type 
        columns = {
            'Nucleotide':'positions', 
            'Sequence': 'sequence',
            f'{sample_type}_read_depth': 'cov',
            f'{sample_type}_effective_depth': 'info',
            f'{sample_type}_mutations': 'sub_N',
            f'{sample_type}_rate': 'sub_rate',
        }
        data = data[columns.keys()].rename(columns=columns).fillna(UKN)
        
        # Create a sequence string
        sequence = ''.join(data['sequence'].tolist()).upper().replace('U', 'T')
        positions = data['positions'].tolist()
        data = data.drop(columns=['sequence', 'positions'])
        
        # Store the data in a Section object
        section = Section('full', 
                          clusters=[Cluster('average', **{c:data[c].tolist() for c in data.columns})],
                          sequence=sequence,
                          positions=positions,
                          num_aligned=int(data['cov'].max()),  
                          )
        
        return cls(name, [section])
    
    
    @classmethod
    def from_rnaf_xml(cls, file):
        # load the data
        assert os.path.isfile(file), f"File {file} does not exist or is not a file."

        with open(file, 'r', encoding='utf-8') as f: 
            my_xml = f.read() 
        transcript = xmltodict.parse(my_xml)['data']['transcript']
        name = transcript['@id'].strip()
        sequence = transcript['sequence'].replace('\t', '').replace('\n', '').upper().replace('U', 'T')
        reactivities = np.fromstring(transcript['reactivity'].replace('\t', '').replace('\n', ''), sep=',')
        reactivities[np.isnan(reactivities)] = np.nan
        assert len(sequence) == len(reactivities), f"Length of sequence ({len(sequence)}) does not match length of reactivities ({len(reactivities)})."
        return cls(name, [Section('full',
                                  [Cluster('average', cov=[np.nan], sub_rate=reactivities.tolist())], 
                                  sequence=sequence,
                                  positions=list(range(1, len(sequence)+1)),
                                  num_aligned=np.nan,
                                  )])


class Section:
    def __init__(self, name, clusters=None, **kwargs):
        self.name = name
        self.kwargs = kwargs
        self.clusters = [] if clusters is None else clusters

    def to_dict(self):
        return {
            **{'#'+k:v for k,v in self.kwargs.items()},
            **{cluster.name: cluster.to_dict() for cluster in self.clusters},
        }
    
    def __repr__(self) -> str:
        return f"Section({self.name}, {self.kwargs}, \n\t\t\t" + '\n\t\t\t'.join([str(cluster) for cluster in self.clusters]) + '\n)'
        
        
class Cluster:
    def __init__(self, name, **kwargs):
        self.name = name
        self.kwargs = kwargs

    def to_dict(self):
        return self.kwargs
    
    @classmethod
    def from_rnaf(cls, name, rnaf):
        return cls(name, **rnaf)

    def __repr__(self) -> str:  
        return f"Cluster({self.name}, \n\t\t\t\t" + ', \n\t\t\t\t'.join([f'{k}={v}' for k,v in self.kwargs.items()]) + ')'