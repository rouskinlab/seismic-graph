import numpy as np
from scipy.stats import pearsonr
import pytest

class FilteredPearson():
    """Calculate the Pearson correlation coefficient between two variables, with the option to filter out significant variations.
    
    Args:
        x: array-like (1D or 2D)
            First variable. 
        y: array-like (1D or 2D)
            Second variable
        thresh: float
            If when removing a pair of points the Pearson correlation coefficient increases by more than this value, the pair is removed. No filtering if None.
            
    Returns:
        score: list
            Pearson correlation coefficient between x and y.
    """
    # HYPERPARAMETERS
    SIGNIFICANT_DIGITS = 3
    
    # if the difference between two datapoints is greater than this line, check what happens when you remove them
    @staticmethod
    def significant_variation_line(x):
        return .1*x + 0.005 
    
    # Now the class
    def __init__(self, x:np.ndarray, y:np.ndarray, thresh:float=None):
        """Calculate the Pearson correlation coefficient between two variables."""
        
        # process input arrays
        mask = np.logical_and(np.logical_and(~np.isnan(x), ~np.isnan(y)), np.logical_and(x != -1000., y != -1000.))
        x = self._array_preprocess(x, mask)
        y = self._array_preprocess(y, mask)
        
        # calculate raw score
        self.raw_score = []
        for i in range(x.shape[1]):
            self.raw_score.append(np.round(pearsonr(x[:,i], y[:,i])[0], self.SIGNIFICANT_DIGITS))
        
        # calculate filter score
        self.filtered_score = self._compute_filtered_pearson(x, y)
        
        # set score
        self.score = []
        for i in range(x.shape[1]):
            if thresh is not None and self.filtered_score[i] - self.raw_score[i] > thresh:
                self.score.append(self.filtered_score[i])
            else:
                self.score.append(self.raw_score[i])
        
    def __str__(self):
        return f'FilteredPearson: {self.score}. Raw score: {self.raw_score}. Filtered score: {self.filtered_score}'
    
    def __call__(self) -> list:
        return self.score
    
    def __getitem__(self, i:int) -> float:
        return self.score[i]
    
    @staticmethod
    def _array_preprocess(v, mask):
        if isinstance(v, list):
            v = np.array(v)
        if not isinstance(v, np.ndarray):
            raise ValueError(f'Expected numpy.ndarray, got {type(v)}')
        if v.ndim == 1:
            v = v.reshape(-1,1)
        v = v[mask, :]
        return v
    
    def _compute_filtered_pearson(self, x:np.ndarray, y:np.ndarray):    
        best_score = []
        for i in range(x.shape[1]):
            # initialize with raw score
            best_score.append(self.raw_score[i]) 
            
            # find significantly divergent points
            significant_variation_idx = np.where(np.abs(x[:,i] - y[:,i]) > self.significant_variation_line(x[:,i]))[0]
            
            # if not significant variations, this loop will not run
            for idx in significant_variation_idx: 
                score = pearsonr(np.delete(x, idx), np.delete(y, idx))[0]
                if score > best_score[-1]:
                    best_score[-1] = np.round(score, self.SIGNIFICANT_DIGITS)
        return best_score 
        
@pytest.mark.parametrize("x, y, thresh, expected", [
    (np.array([1,2,3,4,5]),     np.array([1,2,3,4,5]),       0.5,  [1.]),    # no significant variations
    (np.array([1,2,3,4,5,2]),   np.array([1,2,3,4,5,10]),    0.5,  [1.]),    # one significant variation
    (np.array([1,2,3,4,5,2,1]), np.array([1,2,3,4,5,10,10]), 0.5,  [-.151]), # two significant variations -> should not be removed
    (np.array([1,2,3,4,5,2]),   np.array([1,2,3,4,5,10]),    None, [.178])   # one significant variation, but no threshold -> should not be removed
])
def test_filtered_pearson(x, y, thresh, expected):
    assert FilteredPearson(x, y, thresh).score == expected
