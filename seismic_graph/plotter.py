import pandas as pd
import numpy as np

from .util import *
# from .util import assert_only_one_row, Fit
from .util.misc import assert_only_one_row, Fit
from .util.filtered_pearson import FilteredPearson

import plotly.graph_objects as go
from plotly.offline import plot, iplot

from itertools import cycle
from typing import Tuple, List
from sklearn import metrics
from sklearn.linear_model import LogisticRegression, LinearRegression
from plotly.subplots import make_subplots
import plotly.express as px
from dms_ci import dms_ci
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from .util.normalization import LinFitTable
from io import StringIO

cmap = dict(A="#F09869", C="#8875C7", G="#F7ED8F", T="#99C3EB",
                    N="#f0f0f0")


LIST_COLORS = ['red','green','blue','orange','purple','black','yellow','pink','brown','grey','cyan','magenta']

def mutation_fraction(df, show_ci:bool=False)->dict:
    assert len(df) == 1, "df must have only one row"
    mh = df.iloc[0].copy()
    
    traces, layouts = [], []
    mh['index_selected'] = [i+1 for i in range(len(mh['sequence']))] #TODO[i + 1 for i in mh.index_selected] # index starts at 1
    mh_unrolled = pd.DataFrame({'mut_rate':list(mh.sub_rate), 'base':list(mh.sequence), 'index_reset':list(range(1, 1+len(mh.index_selected))),'index_selected':mh.index_selected})

    err_min, err_max= dms_ci(mh['sub_rate'], len(mh['sub_rate'])*[mh['num_aligned']])
    mr = np.array(mh['sub_rate'])

    for bt in ['T', 'G', 'C', 'A']:
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate':list(df_loc.mut_rate),
                                        'base':list(df_loc.base), 
                                        'index': df_loc['index_selected']})
        traces.append(go.Bar(
            x= np.array(df_loc['index_reset']),
            y= np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text = hover_attr,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
            ))
        if show_ci:
            idx = [i for i, s in enumerate(mh['sequence']) if s == bt]
            traces[-1].update(
                        error_y=dict(
                            type='data',
                            array= (err_max - mr)[idx],
                            arrayminus = (mr - err_min)[idx],
                            visible=True,
                            symmetric=False,
                            thickness=1.5,
                            width=2,
                            color='black'
                        ))

    fig = go.Figure(data=traces)

    fig.update_layout(title=f"{mh['sample']} - {mh['reference']} - {mh['section']} - {mh['cluster']} - {mh['num_aligned']} reads",
                        xaxis=dict(title="Position"),
                        yaxis=dict(title="Mutation fraction"))
   
    fig.update_yaxes(
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True,
    )
    fig.update_xaxes(
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True,
    )
    
    # make the background white 
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'df':mh}


def mutation_fraction_identity(data, show_ci:bool=False)->dict:

    assert len(data) > 0, "The combination of sample, reference and section does not exist in the dataframe"
    assert len(data) == 1, "The combination of sample, reference and section is not unique in the dataframe"
    data = data.iloc[0].copy()
    
    df = pd.DataFrame(index = list(data['sequence']))
    fig = go.Figure()

    data['err_min'] = [dms_ci(p, data['num_aligned'])[0] for p in data['sub_rate']]
    data['err_max'] = [dms_ci(p, data['num_aligned'])[1] for p in data['sub_rate']]

    for base in ['A','C','G','T']:
        df[base] = np.array(data['sub_'+base])/np.array(data['info'])
        fig.add_trace( go.Bar(x=np.arange(1, 1+len(data['sequence'])), y=list(df[base]), marker_color=cmap[base], showlegend=True, name=base) )
    
    # add error bars to stacked_bar[-1]
    if show_ci:
        fig.data[-1].error_y=dict(
                type='data',
                array= [data['err_max'][i]-data['sub_rate'][i] for i in range(len(data['sequence']))],
                arrayminus = [data['sub_rate'][i]-data['err_min'][i] for i in range(len(data['sequence']))],
                visible=True,
                symmetric=False,
                width=2,
                color='black'
        )
        
    fig.update_layout(
        title='Mutation fraction identity - {} - {} - {} - {} - {} reads'.format(data['sample'], data['reference'], data['section'], data['cluster'], data['num_aligned']),
        xaxis_title='Position',
        yaxis_title='Mutation fraction',
        )

    # fig.update_xaxes(tickangle=0, 
    #         tickvals=np.arange(len(df.index)), ticktext=list(df.index), tickfont={'size':8})

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True,
    )
    

    fig.update_layout(barmode='stack')
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'df':df}

    
    
def experimental_variable_across_samples(data:pd.DataFrame, experimental_variable:str, table:LinFitTable, models:List[str]=None, normalize=False)->dict:

    fig = go.Figure()
    
    assert len(data) > 0, "No data to plot"
    assert experimental_variable in data.columns, "Experimental variable not found in data"
    assert len(data['sequence'].unique()) == 1, "More than one sequence found in data. CHeck that reference and section are unique"

    if normalize:
        data = table.normalize_df(data, data['sample'].iloc[0])

    df = pd.DataFrame(
        np.hstack([np.vstack(data['sub_rate'].values), np.vstack(data[experimental_variable].values)]),
        index=data[experimental_variable],
        columns=[c + str(idx+1) for c, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])] + [experimental_variable]
        ).sort_index()
    
    for col in df.columns:
        if col == experimental_variable:
            continue
        fig.add_trace(
            go.Scatter(
                x = df.index,
                y = df[col],
                mode='lines+markers',
                name=col,
                text= experimental_variable,
                hovertemplate = '<b>Experimental variable: %{x}<br>Mutation fraction: %{y}<extra></extra>',
            ),
        )
            
    fig.update_layout(
        title='Mutation rates across experimental variable - {}'.format(experimental_variable),
        xaxis_title=experimental_variable,
        yaxis_title='Mutation fraction',
    )

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True,
    )

    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'data':df}



def auc(df:pd.DataFrame,  savefile=None, auto_open=False, use_iplot=True, title=None)->dict:
    def make_roc_curve(X, y, y_pred, fig, title):
        fpr, tpr, thresholds = metrics.roc_curve(y, y_pred)
        roc_auc = metrics.auc(fpr, tpr)
        fig.add_trace(go.Scatter(x=fpr, y=tpr,
                            mode='lines',
                            name=title,
                            line=dict(width=2)))
        return fig


    fig = go.Figure()
    for row in df.iterrows():
        X = row[1]['sub_rate'].reshape(-1, 1)
        y = np.array([1 if c == '.' else 0 for c in row[1]['structure']]).reshape(-1, 1)
        y_pred = LogisticRegression().fit(X, y.ravel()).predict_proba(X)[:,1]
        fig = make_roc_curve(X, y, y_pred, fig, row[1]['unique_id'])


    fig.add_trace(
        go.Scatter(x=[0, 1], y=[0, 1], line=dict(color='black', width=2, dash='dash'), showlegend=False)
    )
    fig.update_layout(
        title='ROC Curve',
        xaxis_title='False Positive Rate',
        yaxis_title='True Positive Rate')

    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_xaxes(constrain='domain', autorange=True)
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'df':df}


def mutation_fraction_delta(df, table:LinFitTable, normalize=False, savefile=None, auto_open=False, use_iplot=True, title=None)->dict:
    assert len(df) == 2, "df must have 2 row"
    mp_attr = ['sample', 'reference', 'section', 'cluster']
    
    if normalize and len(df['sample'].unique()) > 1:
        df = table.normalize_df(df, df['sample'].iloc[0])
    
    df['unique_id'] = df.apply(lambda row: ' - '.join([str(row[attr]) for attr in mp_attr]), axis=1)

    mh = pd.Series(
        {
            'sub_rate': df['sub_rate'].values[0] - df['sub_rate'].values[1],
            'sequence': ''.join([c1 if c1 == c2 else '-' for c1,c2 in zip(df['sequence'].values[0],df['sequence'].values[1])]),
            'title': "{} - {} reads vs {} - {} reads".format(df['unique_id'].values[0], df['num_aligned'].values[0], df['unique_id'].values[1], df['num_aligned'].values[1])
        }
    )
    
    traces, layouts = [], []
    mh_unrolled = pd.DataFrame({'mut_rate':list(mh.sub_rate), 'base':list(mh.sequence), 'index_reset':list(range(1, 1+len(mh.sequence)))})

    for bt in set(mh['sequence']):
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate':list(df_loc.mut_rate),
                                        'base':list(df_loc.base)})
        traces.append(go.Bar(
            x= np.array(df_loc['index_reset']),
            y= np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text = hover_attr,
            hovertemplate = ''.join(["<b>"+ha+": %{text["+str(i)+"]}<br>" for i, ha in enumerate(hover_attr)]),
            ))
    

    fig = go.Figure(data=traces, 
                    layout=go.Layout(
                        title=go.layout.Title(text=mh['title']),
                        xaxis=dict(title="Position"),
                        yaxis=dict(title="Mutation fraction", range=[0, 0.1])))

    fig.update_yaxes(
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True
    )
    fig.update_xaxes(
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True
    )
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'df':mh}

               
def _mutations_per_read_subplot(data):

    hist = np.sum(np.stack(data.values), axis=0)
    if (hist[0]==0 and len(hist)<2):
        return go.Bar( x=[], y=[], showlegend=False, marker_color='indianred')
    else:

        bin_edges = np.arange(0, np.max(np.argwhere(hist != 0 )), 1)
        return go.Bar( x=bin_edges, y=hist, showlegend=False, marker_color='indianred')

def mutations_per_read_per_sample(data):

    unique_samples = data['sample'].unique()
    fig = make_subplots(rows=len(unique_samples), cols=1, vertical_spacing=0.4/len(unique_samples),
                        subplot_titles=['Number of mutations per read - {}'.format(sample) for sample in unique_samples])
    for i_s, sample in enumerate(unique_samples):
        
        fig.add_trace(_mutations_per_read_subplot(data[data['sample']==sample]['sub_hist'].reset_index(drop=True)),
                      row=i_s+1, col=1 )
        # fig.update_yaxes(title='Count')
        # fig.update_xaxes(dtick=10)
        fig.update_yaxes(
            title='Count',
            gridcolor='lightgray',
            linewidth=1,
            linecolor='black',
            mirror=True,
        )
        fig.update_xaxes(
            dtick=10,
            linewidth=1,
            linecolor='black',
            mirror=True,
            autorange=True,
        )

    fig.update_layout(autosize=True, height=len(unique_samples)*500, 
                      title='Number of mutation per read across samples', 
                      plot_bgcolor='white',paper_bgcolor='white')

    return {
        'fig':fig,
        'data':data
        }
    
def num_aligned_reads_per_reference_frequency_distribution(data):
    """Plot a histogram showing the distribution of aligned reads per reference for each sample and all samples combined.

    Args:
        data (pd.DataFrame): DataFrame containing 'num_aligned' and 'sample' columns.

    Returns:
        dict: {'fig': a Plotly figure, 'data': data}
    """
    # Check that required columns are present
    if 'num_aligned' not in data.columns or 'sample' not in data.columns:
        raise ValueError("Data must contain 'num_aligned' and 'sample' columns.")
    
    # Get unique samples
    unique_samples = data['sample'].unique()
    
    # Prepare the sample list, adding 'All Samples Combined'
    sample_list = list(unique_samples)
    if len(sample_list) > 1:
        sample_list.append('All Samples Combined')
    
    # Collect num_aligned data for each sample
    num_aligned_dict = {}
    for sample in unique_samples:
        sample_data = data[data['sample'] == sample]['num_aligned']
        num_aligned_dict[sample] = sample_data
    
    # If there are multiple samples, add combined data
    if len(unique_samples) > 1:
        combined_data = data['num_aligned']
        num_aligned_dict['All Samples Combined'] = combined_data
    
    # Initialize figure
    fig = go.Figure()
    traceTrack = []
    
    # Add histogram traces for each sample
    for sample in sample_list:
        num_aligned = num_aligned_dict[sample]
        hist_trace = go.Histogram(
            x=num_aligned,
            showlegend=False,
            marker_color='indianred',
            hovertemplate="Number of aligned reads: %{x}<br>Count: %{y}<extra></extra>",
            xbins=dict(size=1000),
            visible=(sample == 'All Samples Combined')  # Make "All Samples Combined" visible by default
        )
        fig.add_trace(hist_trace)
        traceTrack.append(sample)
    
    # Create dropdown menu
    buttons = []
    for idx, sample in enumerate(sample_list):
        visibility = [traceSample == sample for traceSample in traceTrack]
        button = dict(
            label=sample,
            method='update',
            args=[
                {'visible': visibility},
                {
                    'title': f"Distribution of Number of Aligned Reads - {sample}",
                    'xaxis': {'title': "Number of aligned reads"},
                    'yaxis': {'title': "Count"}
                }
            ]
        )
        buttons.append(button)
    
    # Update layout with dropdown menu, setting "All Samples Combined" as the default option
    fig.update_layout(
        updatemenus=[
            dict(
                active=sample_list.index('All Samples Combined'),  # Set default to "All Samples Combined"
                buttons=buttons,
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
        title=f"Distribution of Number of Aligned Reads - All Samples Combined",  # Set initial title to "All Samples Combined"
        xaxis=dict(title="Number of aligned reads"),
        yaxis=dict(title="Count"),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axes
    fig.update_yaxes(
        title='Count',
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True,
    )
    
    return {
        'fig': fig,
        'data': data
    }
    
    
def mutation_per_read_per_reference(data):
    assert len(data) == 1, "data must have 1 row"
    data = data.iloc[0]
    sample, reference = data['sample'], data['reference']
    data = data['sub_hist']
    
    # normalize by the number of reads
    fig = px.bar(x=np.arange(0,len(data)), y=data)

    fig.update_layout(barmode='stack')
    fig.update_layout(title='Number of mutations per read - {} - {} (N={})'.format(sample, reference, np.sum(data)))
    # fig.update_yaxes(title='Count')
    # fig.update_xaxes(title='Number of mutations per read')
    fig.update_yaxes(
        title='Count',
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        title='Number of mutations per read',
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True,
    )
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {
        'fig':fig,
        'data':data
        }



def base_coverage(data):
    
    assert_only_one_row(data)
    data = data.iloc[0]
    
    fig = go.Figure(
        go.Bar(
            x=np.arange(1, 1+len(data['cov'])),
            y=data['cov'],
            showlegend=False,
            marker_color='indianred',
            ),
        layout=go.Layout(
            title=go.layout.Title(text='Base coverage - {} - {}'.format(data['sample'], data['reference'])),
            xaxis=dict(title="Position"),
            yaxis=dict(title="Count")
            )
        )
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True,
    )
    
    return {'fig':fig, 'data':data}



def compare_mutation_profiles(data, table:LinFitTable, normalize=False, max_plots = 100, max_axis=None, pearson_filter_gap = None):
    if normalize:
        data = table.normalize_df(data, data['sample'].iloc[0])

    # total number of plots
    totNumPlots =  int((data.shape[0] *  (1+data.shape[0]))/2)
    if totNumPlots > max_plots:
        print('Too many plots: {} rows combined together make {} plots when using these arguments. Filtered dataset is: \n\n{}.'.format(data.shape[0], totNumPlots, data[['sample','reference','section','cluster','sub_rate_x']]))
        return {'fig': go.Figure(), 'data': data}
    
    if data.shape[0] == 0:
        print('No data found for this combination of arguments')
        return {'fig': go.Figure(), 'data': data}    

    if data.shape[0] == 1:
        print('Only one row found for this combination of arguments.')
        return {'fig': go.Figure(), 'data': data}        

    # Assert sequence is the same for all rows
    assert data['sequence'].nunique() == 1, 'Sequence is not the same for all rows. Select a subset of the data that has the same sequence.'

    # sort by unique_id
    data = data.sort_values('unique_id').reset_index(drop=True)

    # plot these
    fig = go.Figure()
    
    makePlotLabel = lambda x, y: '{} vs {}'.format(x, y)
    
    traceTrack = []
    annotationTrack = []

    if max_axis is not None:
        maxValue = max_axis
    else:
        maxValue = 0
        for idx1, row1 in data.iloc[:-1].iterrows():
            for _, row2 in data.iloc[idx1+1:].iterrows():
                x, y = row1['sub_rate'], row2['sub_rate']
                x, y = np.array(x), np.array(y)
                x, y = x[~np.isnan(x)], y[~np.isnan(y)]
                maxValue = max(x.max(), y.max(), 0.14, maxValue)
    maxValue += 0.01

    for idx1, row1 in data.iloc[:-1].iterrows():
        for _, row2 in data.iloc[idx1+1:].iterrows():
            x, y = row1['sub_rate'], row2['sub_rate']
            x, y = np.array(x), np.array(y)
            mask = np.logical_and(np.logical_and(~np.isnan(x), ~np.isnan(y)), np.logical_and(x != -1000., y != -1000.))
            x, y = np.round(x[mask], 4), np.round(y[mask], 4)
            xlabel, ylabel = row1['unique_id'], row2['unique_id']
            plotLabel = makePlotLabel(xlabel, ylabel)
            
            text = []
            for seq_idx in np.where(mask)[0]:
                text.append(f"position: {seq_idx}<br>base: {data['sequence'].iloc[0][seq_idx]}") 
            
            # plot x vs y then the linear regression line, then the 1:1 line, then the R2 and RMSE values
            fig.add_trace(go.Scatter(x=x, y=y, mode='markers', name='mutation fraction', visible=False, text=text, hovertemplate='%{text} <br>mut frac x: %{x} <br>mut frac y: %{y}'),)
            traceTrack.append(plotLabel)
            
            # linear regression line
            model = LinearRegression()            
            model.fit(x.reshape(-1,1), y)
            slope, intercept, r_value = model.coef_[0], model.intercept_, model.score(x.reshape(-1,1), y)
            fig.add_trace(go.Scatter(x=np.linspace(0, maxValue), y=slope*np.linspace(0, maxValue)+intercept, mode='lines', name='linear regression:<br>y = {}x + {}'.format(round(slope,4), round(intercept,4)), visible=False))
            traceTrack.append(plotLabel)
            
            # 1:1 line
            fig.add_trace(go.Scatter(x=[0, maxValue], y=[0, maxValue], mode='lines', name='Line of Identity', visible=False))
            traceTrack.append(plotLabel)
            
            # Determination, RMSE, Pearson and linear regression line
            coef_of_determination = np.round(r_value, 4)
            pearson = FilteredPearson(x, y, pearson_filter_gap)[0]
            annot = 'Pearson (R) = {}  <br> Coef. of determ. (R\u00B2) = {} <br> RMSE = {} <br> Lin Reg: y = {}x + {} '.format(pearson, coef_of_determination, round(np.sqrt(np.mean((y - (slope*x+intercept))**2)),4), round(slope,4), round(intercept,4))
            fig.add_annotation(visible = False, y=0.13,text=annot, showarrow=False)
            annotationTrack.append(annot)
            
    # add button to select a specific plot using the traceTrack and a dropdown menu
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=list([
                    # make the annotation visible and the plot visible
                    dict(label=plot,
                         method="update",
                         args=[
                            {
                                "visible": [
                                    True if traceTrack[j] == plot else False for j in range(len(traceTrack))
                                ]                                
                            },
                            {   
                                "annotations": [dict(visible=True if np.unique(traceTrack)[j] == plot else False, y=0.13, text=annotationTrack[idx], showarrow=False) for j in range(len(annotationTrack))],
                                "xaxis.title": "{} sub rate".format(plot.split(' vs ')[0]),
                                "yaxis.title": "{} sub rate".format(plot.split(' vs ')[1]),
                                "title": "{} mutation fraction correlation".format(plot)
                            }
                        ]
                    )
                    for idx, plot in enumerate(np.unique(traceTrack))
                ]),
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
    )

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            range=[0, maxValue],
            scaleanchor="y",
            scaleratio=1,
            linecolor='lightgray',
            linewidth=1,
        ),
        yaxis=dict(
            linecolor='lightgray',
            linewidth=1,
        ),
        xaxis_title="{} sub rate".format(xlabel),
        yaxis_title="{} sub rate".format(ylabel),
        # font=dict(
        #     size=15
        # ),
        margin=dict(
            l=0,
            r=0,
            b=50,
            t=50,
            pad=4
        ),
        title = plotLabel + ' mutation fraction correlation',
    )
    

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
        showgrid=True,
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='lightgray',
    )
    fig.update_xaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
        showgrid=True,
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='lightgray',
        # autorange=True,
    )
    # make the first plot visible
    fig.data[0].visible = True
    fig.data[1].visible = True
    fig.data[2].visible = True

    # make the first annotation visible
    fig.layout.annotations[0].visible = True
    
    return {'fig': fig, 'data': data}



def correlation_by_refs_between_samples(df:pd.DataFrame, table:LinFitTable, normalize=False, pearson_filter_gap = None)->dict:
        
    assert len(df['sample'].unique()) == 2, "only two samples are allowed"
    s1, s2 = df['sample'].unique()
    if normalize:
        df = table.normalize_df(df, s1)

    # only keep the references that are shared by both samples
    refs = df.groupby('sample')['reference'].unique()
    refs = list(set(refs.iloc[0]).intersection(set(refs.iloc[1])))
    df = df[df['reference'].isin(refs)]
    
    # calculate correlation
    scores = {}
    for ref in refs:
        df_ref = df[df['reference'] == ref]
        v1, v2 = df_ref[df_ref['sample'] == s1]['sub_rate'].values[0], df_ref[df_ref['sample'] == s2]['sub_rate'].values[0]
        if len(v1) <= 1 or len(v2) <= 1:
            continue

        r = FilteredPearson(v1, v2, pearson_filter_gap)[0]
        r_squared = r**2
        scores[ref] = r_squared
    
    # sort by correlation
    scores_df = pd.DataFrame.from_dict(scores, orient='index', columns=['r_squared']).reset_index().rename(columns={'index':'reference'})
    scores_df = scores_df.sort_values(by='r_squared', ascending=True).reset_index(drop=True)

    scores_df['rank'] = scores_df.index
    
    # Convert the DataFrame to CSV format in a string
    csv_output = StringIO()
    scores_df.to_csv(csv_output, index=False)
    csv_output.seek(0)
    csv_data = csv_output.getvalue()

    print(f"CSV Data:\n")
    print(csv_data)

    # plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=scores_df.index, y=scores_df['r_squared'], text = scores_df['reference'], mode='markers', hovertemplate='reference: %{text}<br>R²: %{y:.2f}'))

    # set layout
    fig.update_layout(
        title = f"R² between samples {s1} and {s2}",
        xaxis_title = "reference rank",
        yaxis_title = "R²",
    )

    return {'fig':fig, 'scores':scores_df, 'scores_csv':csv_data}


from .one_pager_utils import one_pager_html_template, export_fig, make_table

def one_pager(df, xrange=[0, 0.15], plot_height_cov=250, plot_height_count=200, plot_height_bar=225, 
                  plot_width_first_col=500, margin=dict(l=0, r=0, t=25, b=10), plot_width_bar=850, **kwargs):
    
    html_figs = {}
    assert len(df) == 1, "df must have only one row"
    row = df.iloc[0].copy()
    # make html tables
    html_figs['table'] = make_table(row)
    fig = base_coverage(df)['fig']\
        .update_layout(height=plot_height_cov, width=plot_width_first_col, title= 'Coverage')
    fig.update_layout(margin=margin)
    
    html_figs['coverage'] = export_fig(fig, format='html', name='coverage')
    
    # make coverage plot
    fig = mutation_per_read_per_reference(df)['fig']\
        .update_layout(height=plot_height_count, width=plot_width_first_col)
    fig.update_xaxes(range=[-0.5, max(10, np.argwhere(row['sub_hist'] > 0).max())])
    fig.update_layout(title = 'Mutations per read')
    fig.update_layout(margin=margin)
    html_figs['mutations per read'] = export_fig(fig, format='html', name='mutations per read')
    
    # make mutation fraction plot
    fig = mutation_fraction(df, show_ci=True)['fig']\
        .update_layout(height=plot_height_bar, width=plot_width_bar, title = 'Mutation fraction with DMS')
    fig.update_yaxes(range=xrange)
    fig.update_layout(margin=margin)
    html_figs['mutation fraction dms'] = export_fig(fig, format='html', name='mutation fraction dms')

    # make mutation fraction identity plot
    fig = mutation_fraction_identity(df, show_ci=0)['fig']\
        .update_layout(height=plot_height_bar, width=plot_width_bar, title = 'Mutation fraction identity with DMS')
    fig.update_yaxes(range=xrange)
    fig.update_layout(margin=margin)

    # use plotly_cdn to avoid downloading plotly.js
    html_figs['mutation fraction dms identity'] = export_fig(fig, format='html', name='mutation fraction dms identity')
    
    return {'html': one_pager_html_template(html_figs, row), 'data': row}


def dist_of_seq_lengths(data:pd.DataFrame):
    """Plot a histogram of sequence lengths across all rows in the dataframe.

    Args:
        data (pd.DataFrame): DataFrame containing 'sequence' column.

    Returns:
        dict: {'fig': a plotly figure}
    """

    if 'sequence' not in data.columns:
        raise ValueError("Data must contain a 'sequence' column.")
    
    lengths = data['sequence'].str.len()

    fig = go.Figure(
        go.Histogram(
            x=lengths,
            marker_color='indianred',
            hovertemplate="Sequence length: %{x}<br>Count: %{y}<extra></extra>",
            nbinsx=50
        )
    )

    fig.update_layout(
        title="Distribution of Sequence Lengths",
        xaxis=dict(title="Sequence Length"),
        yaxis=dict(title="Count"),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )

    return {
        'fig': fig,
        'data': data
    }

def percent_masked_histogram(data):
    """Plot a histogram showing the distribution of percentage of A's and C's bases masked.

    Args:
        data (pd.DataFrame): DataFrame containing 'cov' and 'sample' columns.

    Returns:
        dict: {'fig': a Plotly figure, 'data': data}
    """
    # Check required columns
    if 'cov' not in data.columns or 'sample' not in data.columns:
        raise ValueError("Data must contain 'cov' and 'sample' columns.")
    
    # Compute masked_percentage for each row
    def compute_masked_percentage(row):
        cov = row['cov']
        sequence = row['sequence']
        if cov is None or len(cov) == 0:
            return np.nan
        if len(sequence) != len(cov):
            ref = row['reference']
            raise ValueError(f"Length of coverage array does not match length of sequence for reference {ref}")
        cov_array = np.array(cov)
        ac_positions = np.isin(list(sequence), ['A', 'C'])
        ac_cov_array = cov_array[ac_positions]
        num_masked = np.isnan(ac_cov_array).sum()
        total_ac_bases = len(ac_cov_array)
        
        return 100.0 * num_masked / total_ac_bases if total_ac_bases > 0 else np.nan

    data['masked_percentage'] = data.apply(compute_masked_percentage, axis=1)
    
    # Get unique samples
    unique_samples = data['sample'].unique()
    
    # Prepare sample list
    sample_list = list(unique_samples)
    if len(sample_list) > 1:
        sample_list.append('All Samples Combined')
    
    # Collect masked_percentage data for each sample
    masked_percentage_dict = {}
    for sample in unique_samples:
        sample_data = data[data['sample'] == sample]
        masked_percentage = sample_data['masked_percentage'].dropna()
        masked_percentage_dict[sample] = masked_percentage

    if len(unique_samples) > 1:
        # Combined data
        masked_percentage_combined = data['masked_percentage'].dropna()
        masked_percentage_dict['All Samples Combined'] = masked_percentage_combined

    # Determine the global x-axis range
    all_masked_percentages = pd.concat(masked_percentage_dict.values())
    x_range = [0, all_masked_percentages.max()]

    # Initialize figure
    fig = go.Figure()
    traceTrack = []
    
    # Add histogram traces for each sample
    for sample in sample_list:
        masked_percentage = masked_percentage_dict[sample]
        hist_trace = go.Histogram(
            x=masked_percentage,
            marker_color='indianred',
            hovertemplate="Percentage of A's & C's Bases Masked: %{x}<br>Count: %{y}<extra></extra>",
            nbinsx=50,  # Adjust bins as needed
            visible=(sample == 'All Samples Combined')  # Make "All Samples Combined" visible by default
        )
        fig.add_trace(hist_trace)
        traceTrack.append(sample)
    
    # Create dropdown menu
    buttons = []
    for idx, sample in enumerate(sample_list):
        visibility = [traceSample == sample for traceSample in traceTrack]
        button = dict(
            label=sample,
            method='update',
            args=[
                {'visible': visibility},
                {
                    'title': f"Masked Percentage Histogram - {sample}",
                    'xaxis': {'title': "Percentage of A & C Masked", 'range': x_range},
                    'yaxis': {'title': "Count"}
                }
            ]
        )
        buttons.append(button)
    
    # Update layout with dropdown menu, setting "All Samples Combined" as the default option
    fig.update_layout(
        updatemenus=[
            dict(
                active=sample_list.index('All Samples Combined'),  # Set default to "All Samples Combined"
                buttons=buttons,
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
        title=f"Masked Percentage Histogram - All Samples Combined",  # Set initial title to "All Samples Combined"
        xaxis=dict(title="Percentage of A & C Masked", range=x_range),
        yaxis=dict(title="Count"),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axes
    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=False,
    )
    
    return {
        'fig': fig,
        'data': data
    }

def f1_violin_by_family(data):
    """Generate a violin plot for each family showing the distribution of F1 scores.

    Args:
        data (pd.DataFrame): DataFrame containing 'family' and 'F1' columns.

    Returns:
        dict: {'fig': a Plotly figure, 'data': data}
    """
    # Check required columns
    if 'family' not in data.columns or 'F1' not in data.columns:
        raise ValueError("Data must contain 'family' and 'F1' columns.")
    
    # Drop rows with missing values in 'family' or 'F1'
    data = data.dropna(subset=['family', 'F1'])

    # Get unique families
    families = data['family'].unique()

    # Initialize figure
    fig = go.Figure()

    # Add a violin trace for each family
    for family in families:
        family_data = data[data['family'] == family]['F1']
        fig.add_trace(
            go.Violin(
                y=family_data,
                name=family,
                box_visible=True,
                meanline_visible=True,
                line_color='black',
                fillcolor='#1f77b4',
                opacity=0.7,
                marker=dict(color='#1f77b4')
            )
        )

    # Update layout
    fig.update_layout(
        title="Distribution of F1 Scores by Family",
        yaxis=dict(title="F1 Score"),
        xaxis=dict(title="Family"),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        tickangle=45,
        tickmode='array',
        tickvals=[family for family in families],
        ticktext=[family for family in families]
    )

    return {
        'fig': fig,
        'data': data
    }