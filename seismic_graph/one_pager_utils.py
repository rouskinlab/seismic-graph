import pandas as pd
import json, os
import plotly.graph_objects as go


margin = dict(l=0, r=0, t=25, b=10)

to_html_args = dict(full_html=False, include_plotlyjs='cdn')


def make_table(row):
    table_items = {
        'sample': "Sample",
        'reference': "Reference",
        'user': "User",
        'date': "Date",
        'exp_env': "Experimental environment",
        'temperature_K': "Temperature [K]",
        'temperature_C': "Temperature [C]",
        'DMS_conc_mM': "DMS concentration [mM]",
        'DMS_concentration_pct_volume': "DMS Concentration",
        'Metabolite': "Metabolite",
        'Modification': "Modification",
    }
    row_lower = {k.lower(): v for k, v in row.items()}

    # Format temperature and DMS values
    if 'temperature_k' in row_lower:
        row_lower['temperature_k'] = f"{row_lower['temperature_k']}&deg; K"

    if 'temperature_c' in row_lower:
        row_lower['temperature_c'] = f"{row_lower['temperature_c']}&deg; C"

    if 'dms_concentration_pct_volume' in row_lower:
        row_lower['dms_concentration_pct_volume'] = f"{row_lower['dms_concentration_pct_volume']}% v/v"

    # Construct the table only for rows that have values
    table_html = '<table>' +\
        '\n'.join([f'<tr><td>{table_items[k]}</td><td>{row_lower[k.lower()]}</td></tr>' 
                   for k in table_items.keys() if k.lower() in row_lower and row_lower[k.lower()]])\
                   + '</table>'

    return table_html

def one_pager_html_template(html_figs, row):
    return f'''
<html>
<head>
<title>{row['sample']} - {row['reference']}</title>
<style>
.row {{
    display: flex;
    flex-direction: row;
    max-width: 1200px;
    max-height: 450px;
}}

table {{
    width: 38%;
    border-collapse: collapse;
    text-align: center;
    border-spacing: 0 3px;

}}
table, th, td {{
    border: 1px solid black;
}}
th, td {{
    padding: 3px;
    text-align: left;
    white-space: nowrap;
}}
.fig-container {{
    margin: 0;
    padding: 0;
}}

</style>
</head>
<body>
<h1>{row['sample']} - {row['reference']}</h1>
<tr>

<div class="row">
<div class="col">
<h2>Experimental info table</h2>
<tr>
<div class="fig-container">{html_figs['table']}</div>
<tr>
</div>

<div class="col">

{html_figs['coverage']}
{html_figs['mutations per read']}
</tr>
</div>
</div>
<tr>
{html_figs['mutation fraction dms']}
<tr>
{html_figs['mutation fraction dms identity']}
<tr>

</body>

</html>
'''

import base64

def export_fig(fig: go.Figure, format: str = 'html', name: str = 'fig'):
    if format == 'html':
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    elif format == 'png':
        img = fig.to_image(format='png', engine='kaleido')
        img = base64.b64encode(img)
        # embed to html
        return f'''<div> 
    /* {name} */
    <img src="data:image/png;base64,{img.decode()}" />
    </div>'''
    else:
        raise ValueError(f'Unknown format: {format}')

