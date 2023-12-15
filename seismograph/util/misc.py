import numpy as np
import datetime
from os.path import exists
import os
import matplotlib.pyplot as plt
from typing import Tuple, List
from plotly.validators.scatter.marker import SymbolValidator
from scipy.optimize import curve_fit
import scipy
import inspect

vals = SymbolValidator().values

def Setshape(x):
        vals = SymbolValidator().values
        return vals[3*x]


def make_path(path:str)->str:
    """Create directories until path exists on your computer. Turns the keyword 'date' into today's date.

    Args:
        path: series of directories that you want to create.
    
    Returns:
        Updated path with today's date instead of the keyword 'date'  
    """

    path = os.path.normpath(path)
    path=path.split(os.sep)
    try:
        path[path.index('date')] = str(datetime.datetime.now())[:10]
    except:
        'No date in path'
    full_path = ''
    for repo in path:
        full_path = full_path + f"{repo}/"
        if not exists(full_path):
            os.mkdir(full_path)
    return full_path



def gini(x:np.array)->float:
    """Returns Gini index

    Args:
        x (np.array): the array you want the Gini index from

    Returns:
        float: Gini index of the input array
    """
    # (Warning: This is a concise implementation, but it is O(n**2)
    # in time and memory, where n = len(x).  *Don't* pass in huge
    # samples!)

    # Mean absolute difference
    mad = np.abs(np.subtract.outer(x, x)).mean()
    # Relative mean absolute difference
    rmad = mad/np.mean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g

def savefig(file:str, close=True)->None:
    """Save a matplotlib figure and create the directory if it doesn't exists.

    Args:
        file: path+title.
        facecolor: color of the background 
    """

    path = make_path('/'.join(file.split('/')[:-1]))
    plt.savefig(path+file.split('/')[-1], bbox_inches='tight')
    if close:
        # Clear the current axes.
        plt.cla() 
        # Clear the current figure.
        plt.clf() 
        # Closes all the figure windows.
        plt.close('all')   


def define_figure(title:str, xlabel:str, ylabel:str, figsize:Tuple[float, float])->plt.figure:
    """Define title, labels and size of your figure.

    Args:
        title: matplotlib title
        xlabel: matplotlib xlabel
        ylabel: matplotlib ylabel
        figsize: matplotlib figsize
    """

    fig = plt.figure(figsize=figsize)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    return fig

import os


class Fit(object):
    def __init__(self) -> None:
        self.legend = ''
    
    def get_legend(self):
        return self.legend

    def predict(self, x, y, model, prefix='', suffix=''):
        fit = self.fit(x,y,model)
        m = eval(model)
        try:
            linreg  = scipy.stats.linregress(y,m(x,*fit))
            self.rvalue = round(linreg.rvalue,5)
        except:
            self.rvalue = 'error'

        self._generate_legend(fit, model, prefix, suffix)
        return np.sort(x), m(np.sort(x),*fit)

    def fit(self, x,y, model):
        fit = curve_fit(eval(model), x, y)[0]
        return fit

    def _generate_legend(self, fit, m, prefix, suffix):
        slice_m = lambda start, stop: ','.join(str(m).split(',')[start:stop])
        first_slice = slice_m(0,len(fit))+','+slice_m(len(fit), len(fit)+1).split(':')[0]
        second_slice = ','.join(m.split(',')[len(fit):])[2:]
        fun_args = [a.strip() for a in str(m).split(',')[1:len(fit)+1]]
        fun_args[-1] = fun_args[-1][0]
        for a,b in zip(fun_args,fit):
            second_slice = second_slice.replace(a.strip(),str(round(b,5)))
        self.legend = prefix+ second_slice + suffix +f'\n R2={self.rvalue}'

                


def save_plot(func,  *args, **kwargs):
    """Default arguments for saving a plotly figure.

    Args:
        to_html (str, optional): File name to save the figure as a html.
        to_png (str, optional): File name to save the figure as a png.
    """

    def wrapper(*args, **kwargs):
        out = func(*args, **kwargs)
        to_html, to_png = kwargs.pop('to_html', None), kwargs.pop('to_png', None)
        if to_html:
            out['fig'].write_html(to_html)
        if to_png:
            out['fig'].write_image(to_png)
        return out
    wrapper.__name__=func.__name__
    wrapper.__doc__=func.__doc__
    return wrapper   

def extract_args(func):
    """Extract the arguments of a function.

    Args:
        func (function): Function to extract the arguments from.

    Returns:
        list: List of the arguments of the function.
    """
    return inspect.getfullargspec(func).args


def assert_only_one_row(df):
    """Assert that the dataframe has only one row.

    Args:
        df (pd.DataFrame): Dataframe to check.

    Raises:
        ValueError: If the dataframe has more than one row.
    """
    if df.shape[0] > 1:
        raise ValueError("The dataframe has more than one row.")
    
    if df.shape[0] == 0:
        raise ValueError("The dataframe is empty.")
