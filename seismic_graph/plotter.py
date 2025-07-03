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
import copy
from scipy.stats import pearsonr

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
    fig.update_traces(marker_line_width=0)

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
    fig.update_traces(marker_line_width=0)
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
    fig.update_traces(marker_line_width=0)
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
    fig.update_traces(marker_line_width=0)

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
    
    # Prepare the sample list, adding 'All Samples Combined' if multiple samples
    sample_list = list(unique_samples)
    has_combined = False
    if len(sample_list) > 1:
        sample_list.append('All Samples Combined')
        has_combined = True
    
    # Collect num_aligned data and counts for each sample
    num_aligned_dict = {}
    sample_counts = {}
    for sample in unique_samples:
        sample_data = data[data['sample'] == sample]['num_aligned']
        num_aligned_dict[sample] = sample_data
        sample_counts[sample] = len(sample_data)
    
    # If there are multiple samples, add combined data
    if has_combined:
        combined_data = data['num_aligned']
        num_aligned_dict['All Samples Combined'] = combined_data
        sample_counts['All Samples Combined'] = len(combined_data)
    
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
            visible=False  # Initially set all traces to not visible
        )
        fig.add_trace(hist_trace)
        traceTrack.append(sample)
    
    # Determine the default view
    if has_combined:
        default_sample = 'All Samples Combined'
    else:
        default_sample = sample_list[0]
    
    # Set the default trace to be visible
    for i, sample in enumerate(sample_list):
        if sample == default_sample:
            fig.data[i].visible = True
            break
    
    # Create dropdown menu
    buttons = []
    for idx, sample in enumerate(sample_list):
        visibility = [traceSample == sample for traceSample in traceTrack]
        # Get the count for this sample
        count = sample_counts.get(sample, 0)
        button = dict(
            label=f"{sample} ({count} rows)",
            method='update',
            args=[
                {'visible': visibility},
                {
                    'title': f"Distribution of Number of Aligned Reads - {sample} ({count} rows)",
                    'xaxis': {'title': "Number of aligned reads"},
                    'yaxis': {'title': "Count"}
                }
            ]
        )
        buttons.append(button)
    
    # Update layout with dropdown menu
    fig.update_layout(
        updatemenus=[
            dict(
                active=sample_list.index(default_sample),  # Set default active sample
                buttons=buttons,
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
        title=f"Distribution of Number of Aligned Reads - {default_sample} ({sample_counts[default_sample]} rows)",  # Set initial title
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
    fig.update_traces(marker_line_width=0)
    
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
    fig.update_traces(marker_line_width=0)
    return {'fig':fig, 'data':data}



def compare_mutation_profiles(
        data,
        table:LinFitTable,
        normalize=False,
        max_plots=100,
        max_axis=None,
        pearson_filter_gap=None,
        positions_to_compare=None
    ):
    if normalize:
        data = table.normalize_df(data, data['sample'].iloc[0])

    # total number of plots
    totNumPlots = int((data.shape[0] * (data.shape[0] - 1)) / 2)
    if totNumPlots > max_plots:
        print('Too many plots: {} rows combined together make {} plots when using these arguments. Filtered dataset is: \n\n{}.'.format(
            data.shape[0], totNumPlots, data[['sample', 'reference', 'section', 'cluster', 'sub_rate_x']]))
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

    # Initialize figure
    fig = go.Figure()

    makePlotLabel = lambda x, y: '{} vs {}'.format(x, y)

    # Create position filter mask if positions_to_compare is provided
    position_mask = None
    if positions_to_compare is not None and len(positions_to_compare) > 0:
        # Convert 1-indexed positions to 0-indexed for array indexing
        positions_0_indexed = [pos - 1 for pos in positions_to_compare if pos > 0]
        sequence_length = len(data['sequence'].iloc[0])
        # Create boolean mask for valid positions
        position_mask = np.zeros(sequence_length, dtype=bool)
        for pos in positions_0_indexed:
            if pos < sequence_length:
                position_mask[pos] = True

    if max_axis is not None:
        maxValue = max_axis
    else:
        maxValue = 0
        for idx1, row1 in data.iloc[:-1].iterrows():
            for _, row2 in data.iloc[idx1 + 1:].iterrows():
                x, y = row1['sub_rate'], row2['sub_rate']
                x, y = np.array(x), np.array(y)
                mask = np.logical_and(~np.isnan(x), ~np.isnan(y))
                
                # Apply position filter if specified
                if position_mask is not None:
                    mask = np.logical_and(mask, position_mask)
                    
                x, y = x[mask], y[mask]
                if len(x) > 0 and len(y) > 0:  # Check if any data remains after filtering
                    maxValue = max(x.max(), y.max(), 0.14, maxValue)
    maxValue += 0.01

    # Generate identical ticks for both axes
    tick_start = 0.0
    tick_end = np.round(maxValue, 2)
    tick_step = np.round((tick_end - tick_start) / 5, 2)
    # guard against zero or extremely small step
    if tick_step <= 0:
        # if tick_end > 0, spread it evenly in 5 steps; otherwise fall back to 0.01
        tick_step = (tick_end - tick_start) / 5 if tick_end > tick_start else 0.01
    tick_vals = np.arange(tick_start, tick_end + tick_step, tick_step)

    # Index to keep track of annotations
    traceTrack = []
    # stats_dict = {}
    annotation_dict = {}
    annotation_index = 0

    for idx1, row1 in data.iloc[:-1].iterrows():
        for idx2, row2 in data.iloc[idx1 + 1:].iterrows():
            x, y = row1['sub_rate'], row2['sub_rate']
            x, y = np.array(x), np.array(y)
            mask = np.logical_and(np.logical_and(~np.isnan(x), ~np.isnan(y)),
                                  np.logical_and(x != -1000., y != -1000.))
            
            # Apply position filter if specified
            if position_mask is not None:
                mask = np.logical_and(mask, position_mask)
                
            x, y = np.round(x[mask], 4), np.round(y[mask], 4)
            xlabel, ylabel = row1['unique_id'], row2['unique_id']
            plotLabel = makePlotLabel(xlabel, ylabel)

            # Skip if no data points remain after filtering
            if len(x) == 0 or len(y) == 0:
                continue

            text = []
            filtered_indices = np.where(mask)[0]
            for seq_idx in filtered_indices:
                # Convert back to 1-indexed for display
                display_position = seq_idx + 1
                text.append(f"position: {display_position}<br>base: {data['sequence'].iloc[0][seq_idx]}")

            # Plot x vs y
            fig.add_trace(go.Scatter(
                x=x, y=y, mode='markers', name='mutation fraction',
                visible=False, text=text,
                hovertemplate='%{text} <br>mut frac x: %{x} <br>mut frac y: %{y}'))

            # Keep track of plot labels for traces
            traceTrack.append(plotLabel)

            # Linear regression line
            model = LinearRegression()
            model.fit(x.reshape(-1, 1), y)
            slope = model.coef_[0]
            intercept = model.intercept_
            r_squared = model.score(x.reshape(-1, 1), y)
            pearson_r = np.corrcoef(x, y)[0, 1]

            # Create regression line
            regression_line_x = np.linspace(0, maxValue)
            regression_line_y = slope * regression_line_x + intercept
            fig.add_trace(go.Scatter(
                x=regression_line_x, y=regression_line_y,
                mode='lines', name=f'Linear regression: y = {round(slope, 4)}x + {round(intercept, 4)}',
                visible=False, line=dict(color='red')))
            traceTrack.append(plotLabel)

            # 1:1 line
            fig.add_trace(go.Scatter(
                x=[0, maxValue], y=[0, maxValue],
                mode='lines', name='Line of Identity',
                visible=False, line=dict(color='black', dash='dash')))
            traceTrack.append(plotLabel)

            # Collect stats for output
            # stats_dict[plotLabel] = {
            #     'Slope': round(slope, 3),
            #     'R2': round(r_squared, 3),
            #     'R': round(pearson_r, 3)
            # }

            # Annotation for statistics
            annot_text = (f'Slope = {round(slope, 3)}<br>'
                          f'R² = {round(r_squared, 3)}<br>'
                          f'Pearson R = {round(pearson_r, 3)}')

            # Store the annotation text for this plotLabel
            annotation_dict[plotLabel] = annot_text

            # Add an annotation to the figure (will control visibility later)
            fig.add_annotation(visible=False, x=0.05 * maxValue, y=0.95 * maxValue,
                               xanchor='left', yanchor='top',
                               text=annot_text, showarrow=False)

            annotation_index += 1

    # Add dropdown to select specific plots
    unique_plots = np.unique(traceTrack)
    buttons = []
    for idx, plot in enumerate(unique_plots):
        # Determine which traces to show
        visibility = []
        num_traces_per_plot = 3  # scatter, regression line, identity line
        for trace_label in traceTrack:
            if trace_label == plot:
                visibility.append(True)
            else:
                visibility.append(False)

        # Set annotations visibility
        annotations_visibility = []
        for i, ann in enumerate(fig.layout.annotations):
            if i == idx:
                annotations_visibility.append(True)
            else:
                annotations_visibility.append(False)

        # Update annotation visibility
        updated_annotations = []
        for i, ann in enumerate(fig.layout.annotations):
            ann_copy = copy.deepcopy(ann)
            ann_copy['visible'] = annotations_visibility[i]
            updated_annotations.append(ann_copy)

        buttons.append(dict(
            label=plot,
            method='update',
            args=[
                {'visible': visibility},
                {
                    'annotations': updated_annotations,
                    'xaxis.title': f"{plot.split(' vs ')[0]} sub rate",
                    'yaxis.title': f"{plot.split(' vs ')[1]} sub rate",
                    'title': f"{plot} mutation fraction correlation"
                }
            ]
        ))

    # Handle case where no valid plots were created
    if len(unique_plots) == 0:
        print('No valid plots could be created with the given position filter.')
        return {'fig': go.Figure(), 'data': data}

    # Determine the default plot (first plot)
    default_plot = unique_plots[0]
    default_visibility = [traceLabel == default_plot for traceLabel in traceTrack]

    # Set default annotation visibility
    default_annotations = []
    for i, ann in enumerate(fig.layout.annotations):
        ann_copy = copy.deepcopy(ann)
        if i == 0:
            ann_copy['visible'] = True
        else:
            ann_copy['visible'] = False
        default_annotations.append(ann_copy)

    # Update figure layout
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=buttons,
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=1000,
        height=1000,
        xaxis=dict(
            range=[0, maxValue],
            dtick=tick_step,
            tickvals=tick_vals,
            ticktext=[str(round(val, 2)) for val in tick_vals],
            showgrid=False,
            ticks='outside',
            showline=True,
            linecolor='black',
            mirror=True,
            zeroline=False,
            # scaleratio=1,
            constrain="domain",
        ),
        yaxis=dict(
            range=[0, maxValue],
            dtick=tick_step,
            tickvals=tick_vals,
            ticktext=[str(round(val, 2)) for val in tick_vals],
            showgrid=False,
            ticks='outside',
            showline=True,
            linecolor='black',
            mirror=True,
            zeroline=False,
            scaleanchor='x',
            scaleratio=1,
            constrain="domain",
        ),
        # margin=dict(
        #     l=0,
        #     r=0,
        #     b=0,
        #     t=0,
        #     pad=4
        # ),
        annotations=default_annotations,
        xaxis_title=f"{default_plot.split(' vs ')[0]} sub rate",
        yaxis_title=f"{default_plot.split(' vs ')[1]} sub rate",
        title=f"{default_plot} mutation fraction correlation"
    )

    # Set the default plot to be visible
    for i, trace in enumerate(fig.data):
        fig.data[i].visible = default_visibility[i]

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
        data (pd.DataFrame): DataFrame containing 'cov', 'sequence', 'reference', and 'sample' columns.

    Returns:
        dict: {'fig': a Plotly figure, 'data': data}
    """
    # Check required columns
    required_columns = ['cov', 'sample', 'sequence', 'reference']
    if not all(col in data.columns for col in required_columns):
        raise ValueError(f"Data must contain columns: {', '.join(required_columns)}")
    
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
    has_combined = False
    if len(sample_list) > 1:
        sample_list.append('All Samples Combined')
        has_combined = True
    
    # Collect masked_percentage data and counts for each sample
    masked_percentage_dict = {}
    sample_counts = {}
    for sample in unique_samples:
        sample_data = data[data['sample'] == sample]
        masked_percentage = sample_data['masked_percentage'].dropna()
        masked_percentage_dict[sample] = masked_percentage
        sample_counts[sample] = len(masked_percentage)
    
    if has_combined:
        # Combined data
        masked_percentage_combined = data['masked_percentage'].dropna()
        masked_percentage_dict['All Samples Combined'] = masked_percentage_combined
        sample_counts['All Samples Combined'] = len(masked_percentage_combined)

    # Restrict bins to be between 0 and 100%
    x_range = [0, 102]

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
            xbins=dict(start=0, end=102, size=2),
            visible=False  # Initially set to False
        )
        fig.add_trace(hist_trace)
        traceTrack.append(sample)
    
    # Determine the default view
    if has_combined:
        default_sample = 'All Samples Combined'
    else:
        default_sample = sample_list[0]
    
    # Set the default trace to be visible
    for i, sample in enumerate(sample_list):
        if sample == default_sample:
            fig.data[i].visible = True
            break
    
    # Create dropdown menu
    buttons = []
    for idx, sample in enumerate(sample_list):
        visibility = [traceSample == sample for traceSample in traceTrack]
        count = sample_counts.get(sample, 0)
        button = dict(
            label=f"{sample} ({count} rows)",
            method='update',
            args=[
                {'visible': visibility},
                {
                    'title': f"Masked Percentage Histogram - {sample} ({count} rows)",
                    'xaxis': {'title': "Percentage of A & C Masked", 'range': x_range},
                    'yaxis': {'title': "Count"}
                }
            ]
        )
        buttons.append(button)
    
    # Update layout with dropdown menu
    fig.update_layout(
        updatemenus=[
            dict(
                active=sample_list.index(default_sample),  # Set default active sample
                buttons=buttons,
                x=0.5,
                y=-0.15,
                xanchor='center',
                yanchor='top',
                direction='up'
            )
        ],
        title=f"Masked Percentage Histogram - {default_sample} ({sample_counts[default_sample]} rows)",  # Set initial title
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

def pearson_correlation_histogram(df: pd.DataFrame) -> dict:
    """Generate a histogram of Pearson R² correlations between two samples for each reference.

    Args:
        df (pd.DataFrame): DataFrame containing the mutation rates.

    Returns:
        dict: Dictionary containing the figure, the computed R² values, and the CSV data.
    """
    # Ensure exactly two samples are selected
    assert len(df['sample'].unique()) == 2, "Exactly two samples are required for this plot."
    s1, s2 = df['sample'].unique()

    refs = df['reference'].unique()
    results = []
    num_plotted = 0

    for ref in refs:
        ref_df = df[df['reference'] == ref]
        num_rows = len(ref_df)
        entry = {
            'Reference': ref,
            'Sample 1': None,
            'Sample 2': None,
            'R2': None,
            'No Replicate Found': False,
            'Insufficient Coverage to Calculate R2': None
        }

        if s1 in ref_df['sample'].values:
            entry['Sample 1'] = s1
        else:
            entry['No Replicate Found'] = True
        if s2 in ref_df['sample'].values:
            entry['Sample 2'] = s2
        else:
            entry['No Replicate Found'] = True

        if num_rows > 2:
            raise ValueError(f"More than two rows found for reference '{ref}'. Ensure only two samples are selected.")
        elif num_rows == 2:
            # Exactly two rows
            v1 = np.array(ref_df.iloc[0]['sub_rate'])
            v2 = np.array(ref_df.iloc[1]['sub_rate'])

            # Mask invalid values
            mask = np.logical_and(
                np.logical_and(~np.isnan(v1), ~np.isnan(v2)),
                np.logical_and(v1 != -1000., v2 != -1000.)
            )
            x = v1[mask]
            y = v2[mask]

            if len(x) < 2:
                # Insufficient data
                entry['Insufficient Coverage to Calculate R2'] = True

            else:
                # Compute Pearson correlation
                entry['Insufficient Coverage to Calculate R2'] = False
                r, _ = pearsonr(x, y)
                r_squared = r ** 2
                entry['R2'] = r_squared
                num_plotted += 1

        results.append(entry)

    # Create DataFrame from results
    results_df = pd.DataFrame(results).sort_values(by='R2', ascending=False, na_position='last')

    # Calculate statistics
    plotted_r2 = results_df['R2'].dropna()
    num_no_replicate = len(results_df[results_df['No Replicate Found'] == True])
    num_insufficient_data = len(results_df[results_df['Insufficient Coverage to Calculate R2'] == True])
    num_references = len(refs)
    num_not_plotted = num_references - num_plotted

    if not plotted_r2.empty:
        highest_r2 = plotted_r2.max()
        median_r2 = plotted_r2.median()
        lowest_r2 = plotted_r2.min()
        # Format R² values to four decimal places
        highest_r2_str = f"{highest_r2:.4f}"
        median_r2_str = f"{median_r2:.4f}"
        lowest_r2_str = f"{lowest_r2:.4f}"
    else:
        highest_r2_str = median_r2_str = lowest_r2_str = "N/A"

    # Prepare annotations for the plot
    stats_text = (
        f"References Plotted: {num_plotted}<br>"
        f"References Not Plotted: {num_not_plotted}<br>"
        f"&nbsp;&nbsp;- Due to No Replicate Found: {num_no_replicate}<br>"
        f"&nbsp;&nbsp;- Due to Insufficient Coverage to Calculate R²: {num_insufficient_data}<br>"
        f"Highest R²: {highest_r2_str}<br>"
        f"Median R²: {median_r2_str}<br>"
        f"Lowest R²: {lowest_r2_str}"
    )

    # Plot histogram with bin size of 0.01 from 0 to 1.01
    fig = go.Figure()
    if num_plotted > 0:
        fig.add_trace(
            go.Histogram(
                x=plotted_r2,
                xbins=dict(start=0, end=1.0001, size=0.010001),  # Define bin size and range
                marker_color='indianred',
                hovertemplate="R²: %{x:.2f}<br>Count: %{y}<extra></extra>"
            )
        )
    else:
        fig.add_trace(go.Histogram())  # Add an empty histogram

    fig.update_layout(
        title=f"Distribution of Pearson R² between samples '{s1}' and '{s2}'",
        xaxis_title="R²",
        yaxis_title="Count",
        plot_bgcolor='white',
        paper_bgcolor='white',
    )
    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
    )
    fig.update_xaxes(
        range=[0, 1.0001],  # Ensure x-axis matches bin range
        dtick=0.1,  # Set ticks to appear every 0.1 interval
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=False,  # Keep the range fixed to the defined bins
    )

    # Add annotation to display statistics
    fig.add_annotation(
        x=0.15,  # Position to the right of the plot
        y=0.5,
        xref='paper',
        yref='paper',
        text=stats_text,
        showarrow=False,
        align='left',
        bordercolor='black',
        borderwidth=1,
        bgcolor='white',
        font=dict(size=12),
    )

    # Convert DataFrame to CSV string
    csv_data = results_df.to_csv(index=False)

    return {'fig': fig, 'data': results_df, 'scores_csv': csv_data}

def compare_mutation_matt(data):
    # if data.shape[0] == 0:
    #     print('No data found for this combination of arguments')
    #     return {'fig': go.Figure(), 'data': data}

    # if data.shape[0] == 1:
    #     print('Only one row found for this combination of arguments.')
    #     return {'fig': go.Figure(), 'data': data}
    
    # # Assert sequence is the same for all rows
    # assert data['sequence'].nunique() == 1, 'Sequence is not the same for all rows. Select a subset of the data that has the same sequence.'

    return {'fig': go.Figure(), 'data': data}
