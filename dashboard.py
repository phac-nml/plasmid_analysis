# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash, sys
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import plotly.express as px
from datetime import date
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
from sklearn.linear_model import LinearRegression

data_sources = {
    'plasmid':'/Users/jrobertson/PycharmProjects/salmonella_plasmid_dynamics/results/plasmid.txt',
    'serovar':'/Users/jrobertson/PycharmProjects/salmonella_plasmid_dynamics/results/serovar.txt',
    'genes':'/Users/jrobertson/PycharmProjects/salmonella_plasmid_dynamics/results/genes.txt'
}


def read_sample_data(file):
    df = pd.read_csv(file,sep="\t",header=0)
    return df

def figure_bar_plot_serotype_plasmid_prevalence(df):
    threshold = 1
    series = df['total_samples']
    perc = (series / series.sum() * 100)

    #build dataframe of larger samples
    df['percent_total_samples'] = perc
    subset = df[df['percent_total_samples'] >= threshold]
    subset = subset[['serovar','count_plasmid_positive_samples','count_plasmid_negative_samples']]

    #get data for all of the remaining rows
    count_samples = len(df) - len(subset)
    sum_plasmid_positive_samples = df['count_plasmid_positive_samples'].sum() - subset['count_plasmid_positive_samples'].sum()
    sum_plasmid_negative_samples = df['count_plasmid_negative_samples'].sum() - subset['count_plasmid_negative_samples'].sum()
    row_name = "Others ({})".format(count_samples)
    subset = subset.rename(columns={'count_plasmid_positive_samples': 'positive', 'count_plasmid_negative_samples': 'negative'})
    new_row = pd.DataFrame([[row_name,sum_plasmid_positive_samples,sum_plasmid_negative_samples]],
                           columns=['serovar','positive','negative'])
    subset = subset.append(new_row)
    print(subset)
    fig = px.bar(subset, x="serovar", y=['positive','negative'],
                 title="Plasmid prevalance by serovar",
                 color_discrete_sequence=px.colors.qualitative.D3,
                 labels={"count_plasmid_positive_samples": "Plasmid positive samples",
                         "count_plasmid_negative_samples": "Plasmid negative samples",
                         "serovar":"Serovar",
                         "value": "Sample count",
                         "variable": "legend"},

                 )
    fig.update_layout(
        font=dict(
            size=18,
         )
    )
    return fig


def figure_histogram_plasmid_prevalence(df):
    fig = px.histogram(df, x="total_samples",log_y=True,)
    return fig




def figure_gene_sunburst(df):
    num_cols = 3
    classes = df.groupby(['drug_class'], sort=True)['total'].sum().sort_values(ascending=False).keys()
    count_classes = len(classes)
    num_rows = int(count_classes / num_cols) + 1

    specs = list()
    for i in range(0,num_rows):
        row = [{"type": "sunburst"}] * num_cols
        specs.append(row)

    #initialize figure
    fig = make_subplots(rows=num_rows, cols=num_cols, specs=specs,horizontal_spacing = 0.001,vertical_spacing=0.0)


    #create a subplot for each drug class
    row_counter = 1
    col_counter = 1
    for c in classes:
        #print("{}\t{}".format(row_counter,col_counter))
        subset = df[df['drug_class'] == c]
        proto = px.sunburst(subset,
                          path=['drug_class', 'gene_id'], values='total',
                          color='proportion',
                          color_continuous_scale='RdBu', color_continuous_midpoint=0.5)

        fig.add_trace(
            go.Sunburst(
                labels=proto['data'][0]['labels'].tolist(),
                parents=proto['data'][0]['parents'].tolist(),
                values=proto['data'][0]['values'],
                marker=proto['data'][0]['marker'],
                branchvalues="total",
                textinfo='label+value',


            ),row=row_counter,col=col_counter
        )

        if col_counter == num_cols:
            col_counter =0
            row_counter+=1
        col_counter += 1
    fig.update_traces(textinfo="label+value")
    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    ),height=1400, width=1000, title_text="Figure 1: Resistance gene associations",
        font=dict(
            family="Arial",
            size=24,
        ),
        coloraxis=dict(colorscale='Portland')

    )
    return fig




external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

gene_df = read_sample_data(data_sources['genes'])
serovar_df = read_sample_data(data_sources['serovar'])
plasmid_df = read_sample_data(data_sources['plasmid'])
#figure_bar_plot_serotype_plasmid_prevalence(serovar_df).show()

figure_histogram_plasmid_prevalence(plasmid_df).show()
sys.exit()

fig = figure_gene_sunburst(gene_df)
fig.show()
sys.exit()



app.layout = html.Div(children=[
    html.H1(children='Salmonella plasmid dynamics interactive visualizations'),
    html.Div(children=[
        html.Div([
            html.H3('A'),
            dcc.Graph(id='f1a', figure=fig)
        ], className="one-third column"),

        html.Div([
            html.H3('B'),
            dcc.Graph(id='f1b', figure=fig)
        ], className="one-third column"),

        html.Div([
            html.H3('C'),
            dcc.Graph(id='f1c', figure=fig)
        ], className="one-third column"),

        html.Div([
            html.H3('D'),
            dcc.Graph(id='f1d', figure=fig)
        ], className="one-third column"),
    ])

], className="row")


if __name__ == '__main__':
    app.run_server(debug=True)