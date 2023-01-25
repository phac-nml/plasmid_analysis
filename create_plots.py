import pandas as pd
import os, logging
import plotly.express as px
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)

def parse_args():
    "Parse the input arguments, use '-h' for help"

    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Visualize the data produced by the analysis tool",
        formatter_class=CustomFormatter)
    parser.add_argument('--genes', type=str, required=True, help='Resistance gene summary information')
    parser.add_argument('--serovar', type=str, required=True, help='Serovar summary information')
    parser.add_argument('--plasmid', type=str, required=True, help='Plasmid summary information')
    parser.add_argument('--vignette',type=str, required=True, help='Specific gene vignette data')
    parser.add_argument('--outdir', type=str, required=True, help='Output results to this directory')
    return parser.parse_args()


def read_sample_data(file):
    return pd.read_csv(file,header=0, encoding = "UTF-8",sep="\t",low_memory=False)

def resistance_gene_sunburst(df):
    path = ['resistance',
            'gene_id', ]

    fig = px.sunburst(
        data_frame=df,
        path=path,
        values='total',
        maxdepth=3,
        color='proportion',
        color_continuous_scale='RdBu',

    )

    fig.update_layout(
        margin=dict(t=1, l=1, r=1, b=1)
    )
    fig.update_traces(textinfo="label+value")

    fig.update_layout(
        coloraxis_colorbar_title='Proportion plasmid',
        coloraxis_colorbar_x=1

    )
    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",
    )

    return fig

def resistance_gene_scatter(df):
    df = df[df['total'] > 10]
    df = df[df['serovar entropy'] >= 0]
    fig = px.scatter(df, x="plasmid entropy", y="serovar entropy",
                     size="total",
                     color='human proportion',
                     hover_name="gene_id", log_x=False, size_max=60, labels={
            "serovar entropy": "Serovar entropy",
            "plasmid entropy": "Plasmid entropy",
            "human proportion": "Human proportion",

        }, )

    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",
    )

    return fig

def plasmid_mobility_scater(df):
    df = df[df['total_samples'] > 10]
    df = df[df['serovar_entropy'] >= 0]
    fig = px.scatter(df, x="serovar_entropy", y="total_samples",
                     color='overall_mobility',
                     labels={
                         "serovar_entropy": "Serovar entropy",
                         "total_samples": "log(10) Total cluster members",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=True, size_max=60, marginal_x="box", marginal_y="box")

    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ))
    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",
    )
    return fig

def plasmid_resistance_scatter(df):
    df = df[df['total_samples'] > 10]
    df = df[df['serovar_entropy'] >= 0]
    fig = px.scatter(df, x="serovar_entropy", y="proportion_resistant",
                     size='total_samples', color='overall_mobility',
                     labels={
                         "proportion_resistant": "Proportion resistant",
                         "serovar_entropy": "Serovar entropy",
                         'overall_mobility': "Overall mobility"
                     },
                     hover_name="plasmid_id", log_y=False, size_max=60)

    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1)
    )
    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",
    )
    return fig

def serovar_plasmid_frac_bar(df,outfile):
    df.to_csv(outfile,header=True,sep='\t')
    fig = px.bar(df, x="serovar",
                 y=["positive", "negative"], log_y=True,
                 color_discrete_map={
                     "positive": "#2978A0", "negative": "#F17300"
                 },
                 )
    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        legend_title="Plasmid presence",
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",

    )
    fig.update_yaxes(title="log(10) Sample count")
    fig.update_xaxes(title="Serovar")
    return fig

def collapse_serovars(df,num_serovars=10):
    df = df.sort_values('total_samples',ascending = False)
    df.reset_index(inplace=True)
    counts = {}
    num_labels = 0
    positive = 0
    negative = 0
    for index,row in df.iterrows():
        serovar = row['serovar']
        count_plasmid_negative_samples = row['count_plasmid_negative_samples']
        count_plasmid_positive_samples = row['count_plasmid_positive_samples']
        if index < num_serovars:
            counts[serovar] = {'serovar':serovar,'positive':count_plasmid_positive_samples,'negative':count_plasmid_negative_samples}
        else:
            num_labels+=1
            positive+=count_plasmid_positive_samples
            negative+=count_plasmid_negative_samples
    counts["{}_others".format(num_labels)] = {'serovar':"{}_others".format(num_labels),'positive':positive,'negative':negative}
    return pd.DataFrame.from_dict(counts,orient='index')

def vignette_sunburst(df):
    path = ['serovar',
            'molecule_type','primary_id','secondary_id' ]

    fig = px.sunburst(
        data_frame=df,
        path=path,
        values='count',
        maxdepth=4,

    )

    fig.update_layout(
        margin=dict(t=1, l=1, r=1, b=1)
    )


    fig.update_layout(
        font=dict(
            size=18,
        ),
    )
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        template="simple_white",

    )
    fig.update_traces(leaf=dict(opacity=1),    marker_line_width=20)
    return fig



def main():
    args = parse_args()

    outdir = args.outdir
    # initialize analysis directory
    if not os.path.isdir(args.outdir):
        logging.info("Creating output directory {}".format(args.outdir))
        os.makedirs(args.outdir, 0o755)

    #Figure files
    figure_1 = os.path.join(outdir,"Figure.1.html")
    figure_2 = os.path.join(outdir, "Figure.2.html")
    figure_3 = os.path.join(outdir, "Figure.3.html")
    figure_4 = os.path.join(outdir, "Figure.4.html")
    figure_5 = os.path.join(outdir, "Figure.5.html")
    figure_6 = os.path.join(outdir, "Figure.6.html")
    serovar_perc_file = os.path.join(outdir, "plasmid.frac.txt")

    #Read data
    gene_df = read_sample_data(args.genes)
    serovar_df = read_sample_data(args.serovar)
    plasmid_df = read_sample_data(args.plasmid)
    spec_gene_df = read_sample_data(args.vignette)

    #Figure 1 - Plasmid prevelance in different serotypes
    fig1 = serovar_plasmid_frac_bar(collapse_serovars(serovar_df),serovar_perc_file)

    fig1.write_image(os.path.join(outdir, "Figure.1.png"), scale=6, width=1600, height=900 )
    fh = open(figure_1,'w')
    fh.write(fig1.to_html())
    fh.close()

    #Figure 2 - Plasmid abundance and serovar entropy based on mobility
    fig2 = plasmid_mobility_scater(plasmid_df)
    #fig2.update_layout(xaxis_range=[-0.1, 4.5])
    fig2.write_image(os.path.join(outdir, "Figure.2.png"), scale=6, width=1600, height=900 )
    fh = open(figure_2,'w')
    fh.write(fig2.to_html())
    fh.close()


    #Figure 3 - Sunburst chart of gene abundance and plasmid fraction
    fig3 = resistance_gene_sunburst(gene_df)
    fig3.write_image(os.path.join(outdir, "Figure.3.png"), scale=6, width=1600, height=1600 )
    fh = open(figure_3,'w')
    fh.write(fig3.to_html())
    fh.close()

    #Figure 4 - Serovar and plasmid entropy of resistance genes
    fig4 = resistance_gene_scatter(gene_df)
    fig4.update_layout(xaxis_range=[-0.1, 3.5])
    fig4.write_image(os.path.join(outdir, "Figure.4.png"), scale=6, width=1600, height=900 )

    fh = open(figure_4,'w')
    fh.write(fig4.to_html())
    fh.close()

    #Figure 5 - blaCMY-2 Sunburst chart
    fig5 = vignette_sunburst(spec_gene_df)
    fig5.write_image(os.path.join(outdir, "Figure.5.png"), scale=6, width=1600, height=1600 )
    fig5.write_image(os.path.join(outdir, "Figure.5.pdf"), scale=6, width=1600, height=1600)
    fh = open(figure_5,'w')
    fh.write(fig5.to_html())
    fh.close()

    #Figure 6 - Plasmid serovar entropy resistance fraction
    fig6 = plasmid_resistance_scatter(plasmid_df)
    fig6.write_image(os.path.join(outdir, "Figure.6.png"), scale=6, width=1600, height=900 )
    fh = open(figure_6,'w')
    fh.write(fig6.to_html())
    fh.close()

    return


if __name__== '__main__':
    main()

