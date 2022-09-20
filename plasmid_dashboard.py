# -*- coding: utf-8 -*-

# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
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



def read_sample_data(file):
    df = pd.read_csv(file,sep="\t",header=0)
    return df


def testing(df):
    path = ['drug_class',
            'gene_id', ]
    fig = make_subplots(rows=1,cols=2)
    fig.add_trace( go.sunburst(
        data_frame=df,
        path=path,
        values='total',
        maxdepth=3,
        color='proportion',
        color_continuous_scale='RdBu'
    ))

    fig.add_trace( go.sunburst(
        data_frame=df,
        path=path,
        values='total',
        maxdepth=3,
        color='proportion'
    ))

    fig.update_layout(
        font=dict(
            size=18,
        )
    )
    fig.show()
    return fig

def create_sunburst(df):

    path = ['drug_class',
                'gene_id',]


    fig = px.sunburst(
        data_frame=df ,
        path=path ,
        values='total',
        maxdepth=3,
        color='proportion',
        color_continuous_scale='RdBu'
    )

    fig.update_layout(
        margin=dict(t=1, l=1, r=1, b=1)
    )
    fig.update_traces(textinfo="label+value")

    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=18,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    fig.show(config=config)
    return fig

def scater_plot(df):

    fig = px.scatter(df[df['plasmid']>10], x="plasmid entropy", y="serovar entropy",
                     size="total",
                     color='human proportion',
                     hover_name="gene_id", log_x=False,size_max=60,labels={
                         "serovar entropy": "Serovar entropy",
                         "plasmid entropy": "Plasmid entropy",
                        "human proportion": "Human proportion",

                     },)
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=18,
        )
    )
    fig.show(config=config)
    return

def scater_plot_chao1(df):

    fig = px.scatter(df[df['plasmid_chao1']>=0], x="total_samples", y="plasmid_chao1",
                     color='proportion_resistant',
                     hover_name="serovar", log_x=False,trendline="ols",                     )
    fig.update_layout(
        font=dict(
            size=18,
        )
    )
    fig.show()
    return


def scater_plot_mobility(df):
    print(stats.ttest_ind(df[df['overall_mobility'] == 'mobilizable']['serovar_entropy'],df[df['overall_mobility'] == 'non-mobilizable']['serovar_entropy']))
    print(stats.ttest_ind(df[df['overall_mobility'] == 'conjugative']['serovar_entropy'], df[df['overall_mobility'] == 'non-mobilizable']['serovar_entropy']))
    print(stats.ttest_ind(df[df['overall_mobility'] == 'mobilizable']['serovar_entropy'], df[df['overall_mobility'] == 'conjugative']['serovar_entropy']))
    df = df[df['serovar_entropy'] >=0]
    fig = px.scatter(df, x="serovar_entropy", y="total_samples",
                     color='overall_mobility',
                     labels={
                         "serovar_entropy": "Serovar entropy",
                         "total_samples": "log(10) Total cluster members",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=True,size_max=60,marginal_x="violin")
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=18,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ))
    fig.show(config=config)
    return

def scater_plot_plasmid_resistance(df):
    fig = px.scatter(df, x="proportion_resistant", y="serovar_entropy",
                     size='total_samples',color='overall_mobility',
                     labels={
                         "proportion_resistant": "Proportion resistant",
                         "serovar_entropy": "Serovar entropy",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=False,size_max=60)
    config = {
        'toImageButtonOptions': {
            'format': 'png',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 600,
            'width': 1200,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.update_layout(
        font=dict(
            size=18,
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1)
    )
    fig.show(config=config)
    return


def create_violin_plot(df):

    fig = px.violin(df, y="serovar_entropy", x="overall_mobility",
                    hover_data=df.columns)
    fig.show()
    return



plasmid_mob_serovar_df = read_sample_data('/Users/jrobertson/Desktop/plasmid_serovar_entropy.txt')
scater_plot_mobility(plasmid_mob_serovar_df)
#serovar_data_df = read_sample_data('/Users/jrobertson/Desktop/serovar_plasmid_info.txt')
#scater_plot_chao1(serovar_data_df)

scater_plot_plasmid_resistance(plasmid_mob_serovar_df)
create_violin_plot(plasmid_mob_serovar_df)
df = read_sample_data('/Users/jrobertson/PycharmProjects/salmonella_plasmid_dynamics/results/genes.txt')
gene_df = read_sample_data('/Users/jrobertson/Desktop/salmonella_res_genes.txt')
create_sunburst(df)
scater_plot(gene_df)
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)




#if __name__ == '__main__':
    #app.run_server(debug=True)