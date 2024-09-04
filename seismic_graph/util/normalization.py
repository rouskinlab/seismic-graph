import numpy as np

class LinFitTable:
    def __init__(self, df):
        self.df = df
        if type(df) == type(None) or len(df) == 0:
            self.samples = []
            self.matrix = np.zeros((0, 0))
        else:
            self.samples = list(df['sample'].unique())
            self.matrix = self._build_lin_reg_matrix()
        
    def __repr__(self) -> str:
        return self.matrix.__repr__()

    @staticmethod   
    def _extract_values_from_sample(df, sample):
        return np.concatenate(df[df['sample'] == sample]['sub_rate'].values)
    
    def _lin_reg_between_samples(self, df, sampleA, sampleB):
        """Fit B = slope*A and return slope"""
        df = df[df['sample'].isin([sampleA, sampleB])]
        df = df.groupby(['reference', 'section']).filter(lambda x: len(x) == 2)
        if not len(df):
            return np.nan
            # raise ValueError(f'No common references between {sampleA} and {sampleB}')
        valuesA = self._extract_values_from_sample(df, sampleA)
        valuesB = self._extract_values_from_sample(df, sampleB)
        # remove nans
        mask = ~np.isnan(valuesA) & ~np.isnan(valuesB)
        valuesA = valuesA[mask]
        valuesB = valuesB[mask]
        slope = np.dot(valuesA, valuesB) / np.dot(valuesA, valuesA) 
        if np.isnan(slope):
            raise ValueError(f'Linear regression between {sampleA} and {sampleB} failed')
        return slope

    def _build_lin_reg_matrix(self):
        n = len(self.samples)
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j:
                    matrix[i, j] = 1
                else:
                    matrix[i, j] = self._lin_reg_between_samples(self.df, self.samples[i], self.samples[j])
        return matrix
    
    def normalize_array(self, array, ref_sample, sample):
        """Normalize sample to ref_sample"""
        slope = self.matrix[self.samples.index(ref_sample), self.samples.index(sample)]
        if np.isnan(slope):
            raise ValueError(f'Couldn\'t read a matrix value for {ref_sample} and {sample}')
        return array / slope 

    def normalize_df(self, df, ref_sample):
        """Normalize all samples to reference"""
        df['sub_rate'] = df.apply(lambda row: self.normalize_array(row['sub_rate'], ref_sample, row['sample']), axis=1)
        return df