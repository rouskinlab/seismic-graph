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
            If removing a pair of points increases the Pearson correlation coefficient by more than this value, the pair is removed. No filtering if None.
            
    Returns:
        score: list
            Pearson correlation coefficient between x and y.
    """
    # HYPERPARAMETERS
    SIGNIFICANT_DIGITS = 3
    
    @staticmethod
    def significant_variation_line(x):
        """Define the threshold for significant variation."""
        return 0.1 * x + 0.005 
    
    def __init__(self, x: np.ndarray, y: np.ndarray, thresh: float = None):
        """Initialize and compute the Pearson correlation coefficient."""
        # Process input arrays
        mask = np.logical_and(
            np.logical_and(~np.isnan(x), ~np.isnan(y)),
            np.logical_and(x != -1000.0, y != -1000.0)
        )
        x = self._array_preprocess(x, mask)
        y = self._array_preprocess(y, mask)
        
        # Calculate raw score
        self.raw_score = []
        for i in range(x.shape[1]):
            xi = x[:, i]
            yi = y[:, i]
            if len(xi) < 2:
                self.raw_score.append(np.nan)
                continue
            r, _ = pearsonr(xi, yi)
            self.raw_score.append(np.round(r, self.SIGNIFICANT_DIGITS))
        
        # Calculate filtered score
        self.filtered_score = self._compute_filtered_pearson(x, y)
        
        # Set final score based on threshold
        self.score = []
        for i in range(x.shape[1]):
            if np.isnan(self.raw_score[i]):
                self.score.append(np.nan)
            elif thresh is not None and self.filtered_score[i] - self.raw_score[i] > thresh:
                self.score.append(self.filtered_score[i])
            else:
                self.score.append(self.raw_score[i])
            
    def __str__(self):
        return f'FilteredPearson: {self.score}. Raw score: {self.raw_score}. Filtered score: {self.filtered_score}'
    
    def __call__(self) -> list:
        return self.score
    
    def __getitem__(self, i: int) -> float:
        return self.score[i]
    
    @staticmethod
    def _array_preprocess(v, mask):
        if isinstance(v, list):
            v = np.array(v)
        if not isinstance(v, np.ndarray):
            raise ValueError(f'Expected numpy.ndarray, got {type(v)}')
        if v.ndim == 1:
            v = v.reshape(-1, 1)
        v = v[mask, :]
        return v
    
    def _compute_filtered_pearson(self, x: np.ndarray, y: np.ndarray):
        best_score = []
        for i in range(x.shape[1]):
            xi = x[:, i]
            yi = y[:, i]
            # Initialize with raw score
            best_score.append(self.raw_score[i])
            
            if np.isnan(self.raw_score[i]):
                continue  # Skip if not enough data points
            
            # Find significantly divergent points
            significant_variation_idx = np.where(
                np.abs(xi - yi) > self.significant_variation_line(xi)
            )[0]
            
            # Try removing each significant variation to see if it improves the score
            for idx in significant_variation_idx:
                xi_filtered = np.delete(xi, idx)
                yi_filtered = np.delete(yi, idx)
                if len(xi_filtered) < 2:
                    continue  # Skip if not enough data points after deletion
                r, _ = pearsonr(xi_filtered, yi_filtered)
                if r > best_score[-1]:
                    best_score[-1] = np.round(r, self.SIGNIFICANT_DIGITS)
        return best_score
        
@pytest.mark.parametrize("x, y, thresh, expected", [
    (np.array([1,2,3,4,5]),     np.array([1,2,3,4,5]),       0.5,  [1.]),    # no significant variations
    (np.array([1,2,3,4,5,2]),   np.array([1,2,3,4,5,10]),    0.5,  [1.]),    # one significant variation
    (np.array([1,2,3,4,5,2,1]), np.array([1,2,3,4,5,10,10]), 0.5,  [-.151]), # two significant variations -> should not be removed
    (np.array([1,2,3,4,5,2]),   np.array([1,2,3,4,5,10]),    None, [.178])   # one significant variation, but no threshold -> should not be removed
])
def test_filtered_pearson(x, y, thresh, expected):
    assert FilteredPearson(x, y, thresh).score == expected
