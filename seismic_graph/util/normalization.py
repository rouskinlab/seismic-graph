import numpy as np

class LinFitTable:
    def __init__(self, df):
        self.df = df
        if df is None or len(df) == 0:
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
        df_pair = df[df['sample'].isin([sampleA, sampleB])]
        # Filter for common 'reference' and 'section'
        df_common = df_pair.groupby(['reference', 'section']).filter(lambda x: set(x['sample']) == {sampleA, sampleB})
        if df_common.empty:
            # No common references between samples
            return np.nan
        valuesA = self._extract_values_from_sample(df_common, sampleA)
        valuesB = self._extract_values_from_sample(df_common, sampleB)
        # Remove NaNs
        mask = ~np.isnan(valuesA) & ~np.isnan(valuesB)
        valuesA = valuesA[mask]
        valuesB = valuesB[mask]
        if len(valuesA) == 0:
            # No valid data points to compute slope
            return np.nan
        numerator = np.dot(valuesA, valuesB)
        denominator = np.dot(valuesA, valuesA)
        if denominator == 0:
            # Avoid division by zero
            return np.nan
        return numerator / denominator

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
        if np.isnan(slope) or slope == 0:
            # Can't compute slope; proceed without normalizing
            print(f'Warning: Cannot normalize sample "{sample}" to reference sample "{ref_sample}" due to lack of common data.')
            return array
        return array / slope 

    def normalize_df(self, df, ref_sample):
        """Normalize all samples to reference"""
        df['sub_rate'] = df.apply(lambda row: self.normalize_array(row['sub_rate'], ref_sample, row['sample']), axis=1)
        return df