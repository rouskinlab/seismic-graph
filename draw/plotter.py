import pandas as pd
import numpy as np

from .util import *
# from .util import assert_only_one_row, Fit
from .util.misc import assert_only_one_row, Fit

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


cmap = dict(A="#F09869", C="#8875C7", G="#F7ED8F", T="#99C3EB",
                    N="#f0f0f0")


LIST_COLORS = ['red','green','blue','orange','purple','black','yellow','pink','brown','grey','cyan','magenta']

def mutation_fraction(df, show_ci:bool=True)->dict:
    assert len(df) == 1, "df must have only one row"
    mh = df.iloc[0].copy()
    
    traces, layouts = [], []
    mh['index_selected'] = [i+1 for i in range(len(mh['sequence']))] #TODO[i + 1 for i in mh.index_selected] # index starts at 1
    mh_unrolled = pd.DataFrame({'mut_rate':list(mh.sub_rate), 'base':list(mh.sequence), 'index_reset':list(range(1, 1+len(mh.index_selected))),'index_selected':mh.index_selected})

    for bt in set(mh['sequence']):
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
            traces[-1].update(
                        error_y=dict(
                        type='data',
                        symmetric=False,
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
                thickness=1.5,
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

    
def experimental_variable_across_samples(data:pd.DataFrame, experimental_variable:str, models:List[str]=[])->dict:

    fig = go.Figure()
    
    assert len(data) > 0, "No data to plot"
    assert experimental_variable in data.columns, "Experimental variable not found in data"
    assert len(data['sequence'].unique()) == 1, "More than one sequence found in data. CHeck that reference and section are unique"

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
    fig.update_xaxes(constrain='domain')
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {'fig':fig, 'df':df}


def mutation_fraction_delta(df, savefile=None, auto_open=False, use_iplot=True, title=None)->dict:
    assert len(df) == 2, "df must have 2 row"
    mp_attr = ['sample', 'reference', 'section', 'cluster']
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

    fig.update_layout(autosize=True, height=len(unique_samples)*500, title='Number of mutation per read across samples')

    
    fig.update_layout(plot_bgcolor='white',paper_bgcolor='white')

    return {
        'fig':fig,
        'data':data
        }
    

def num_aligned_reads_per_reference_frequency_distribution(data):
    assert len(samples:=data['sample'].unique()) == 1, "data must have 1 sample"
    data = data['num_aligned'].values
    fig = go.Figure(
            go.Histogram(
                x=data, 
                showlegend=False, 
                marker_color='indianred',
                hovertemplate="Number of aligned reads: %{x}<br>Count: %{y}<extra></extra>"
                ),
            layout=go.Layout(
                title=go.layout.Title(text='{} - Reads per reference count'.format(samples[0])),
                xaxis=dict(title="Number of aligned reads"),
                yaxis=dict(title="Count"),
                plot_bgcolor='white',
                paper_bgcolor='white'
                )          
            )
    
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
    return {
        'fig' : fig,
        'data' : data
    }
    

# def num_aligned_reads_per_reference_frequency_distribution(data):
#     assert len(samples:=data['sample'].unique()) == 1, "data must have 1 sample"
#     data = data['num_aligned'].values
#     return {
#             'fig':go.Figure(
#                 go.Histogram(
#                     x=data, 
#                     showlegend=False, 
#                     marker_color='indianred',
#                     hovertemplate="Number of aligned reads: %{x}<br>Count: %{y}<extra></extra>"
#                     ),
#                 layout=go.Layout(
#                     title=go.layout.Title(text='{} - Reads per reference count'.format(samples[0])),
#                     xaxis=dict(title="Number of aligned reads"),
#                     yaxis=dict(title="Count"),
#                     plot_bgcolor='white',
#                     paper_bgcolor='white'
#                     )          
#                 ),
#             'data':data
#             }
    
    
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



def compare_mutation_profiles(data, max_plots = 100, max_axis=None):

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
            x, y = x[~np.isnan(x)], y[~np.isnan(y)]
            xlabel, ylabel = row1['unique_id'], row2['unique_id']
            plotLabel = makePlotLabel(xlabel, ylabel)
            
            # plot x vs y then the linear regression line, then the 1:1 line, then the R2 and RMSE values
            fig.add_trace(go.Scatter(x=x, y=y, mode='markers', name='mutation fraction', text=plotLabel, visible=False),)
            traceTrack.append(plotLabel)
            
            # linear regression line
            model = LinearRegression()            
            model.fit(x.reshape(-1,1), y)
            slope, intercept, r_value = model.coef_[0], model.intercept_, model.score(x.reshape(-1,1), y)
            fig.add_trace(go.Scatter(x=x, y=slope*x+intercept, mode='lines', name='linear regression: y = {}x + {}'.format(round(slope,4), round(intercept,4)), visible=False))
            traceTrack.append(plotLabel)
            
            # 1:1 line
            fig.add_trace(go.Scatter(x=[0, maxValue], y=[0, maxValue], mode='lines', name='Line of Identity', visible=False))
            traceTrack.append(plotLabel)
            
            # R2, RMSE and linear regression line
            annot = 'R2 = {} <br> RMSE = {} <br> Lin Reg: y = {}x + {}'.format(round(r_value**2,4), round(np.sqrt(np.mean((y - (slope*x+intercept))**2)),2), round(slope,4), round(intercept,4))
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
                x=1,
                y=0.4,
                xanchor='left',
                yanchor='top',
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



def correlation_by_refs_between_samples(df:pd.DataFrame)->dict:
    assert len(df['sample'].unique()) == 2, "only two samples are allowed"
    s1, s2 = df['sample'].unique()

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
        corr, _ = pearsonr(v1, v2)
        scores[ref] = corr
    
    # sort by correlation
    scores = pd.DataFrame.from_dict(scores, orient='index', columns=['correlation']).reset_index().rename(columns={'index':'reference'})
    scores = scores.sort_values(by='correlation', ascending=True).reset_index(drop=True)

    # plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=scores.index, y=scores['correlation'], text = scores['reference'], mode='markers', hovertemplate='reference: %{text}<br>correlation: %{y:.2f}'))

    # set layout
    fig.update_layout(
        title = f"correlation between samples {s1} and {s2}",
        xaxis_title = "reference rank",
        yaxis_title = "correlation",
    )

    return {'fig':fig, 'scores':scores}