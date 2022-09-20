#!/usr/bin/env python
import logging
import os
import sys
from argparse import (ArgumentParser, FileType)
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score,adjusted_rand_score
from scipy.stats import entropy
from statistics import mean
from skbio.diversity.alpha import simpson
from collections import Counter
import scipy.stats as stats
import skbio.diversity.alpha as alpha
import scipy


def calc_MI(category_1,category_2):
    return adjusted_mutual_info_score(category_1,category_2,average_method='arithmetic')

def calc_fisherExact(contingency_table):
    return stats.fisher_exact(contingency_table)

def calc_shanon_entropy(value_list):

    total = sum(value_list)
    values = []
    for v in value_list:
        values.append(v/total)
    return entropy(values)

def calc_homogeneity(value_list):
    return 1 - simpson(value_list)

def calc_adjusted_rand(category_1,category_2):
    return adjusted_rand_score(category_1,category_2)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Filter overlapping queries')
    parser.add_argument('--abricate', type=str, required=True, help='Abricate report')
    parser.add_argument('--contigs', type=str, required=True, help='MOB-suite contig report')
    parser.add_argument('--mobtyper', type=str, required=True, help='MOB-suite mobtyper report')
    parser.add_argument('--metadata', type=str, required=True, help='Sample metadata')
    parser.add_argument('--outdir', type=str, required=True, help='Output results to this directory')
    return parser.parse_args()


def read_data_into_df(file):
    return pd.read_csv(file,header=0, encoding = "UTF-8",sep="\t")

def process_abricate(file):
    data = read_data_into_df(file)
    records = {}
    for index, row in data.iterrows():

        sample_id = row['#FILE']
        contig_id = str(row['SEQUENCE'])
        gene_id = row['GENE']
        resistance = row['RESISTANCE']

        if not sample_id in records:
            records[sample_id] = {}
        if not contig_id in records[sample_id]:
            records[sample_id][contig_id] = {}
        records[sample_id][contig_id][gene_id] = resistance
    return records

def process_contigs(file):
    data = read_data_into_df(file)
    records = {}
    for index, row in data.iterrows():
        sample_id = row['sample_id']
        contig_id = str(row['contig_id'].split("_")[0])
        molecule_type = row['molecule_type']
        size = row['size']
        rep_type = row['rep_type(s)']
        mob_type = row['relaxase_type(s)']
        if mob_type != '-'  and row['mpf_type'] != '-':
            predicted_mobility = 'conjugative'
        elif mob_type != '-':
            predicted_mobility = 'mobilizable'
        else:
            predicted_mobility = '-'

        primary_cluster_id = row['primary_cluster_id']
        secondary_cluster_id = row['secondary_cluster_id']

        if not sample_id in records:
            records[sample_id] = {'chromosome':{},'plasmid':{}}
        if not molecule_type in records[sample_id]:
            records[sample_id][molecule_type] = {}
        if molecule_type == 'chromosome':
            records[sample_id][molecule_type][contig_id] = {'size':size,'genes':[]}
        else:
            records[sample_id][molecule_type][contig_id] = {'size': size,
                                                            'rep_type':rep_type,
                                                            'mob_type':mob_type,
                                                            'predicted_mobility':predicted_mobility,
                                                            'primary_cluster_id':primary_cluster_id,
                                                            'secondary_cluster_id':secondary_cluster_id,
                                                            'genes':[]}

    return records

def process_mobtyper(file):
    data = read_data_into_df(file)
    records = {}
    for index, row in data.iterrows():
        sample_id = row['sample_id']
        plasmid_id = row['plasmid_id']
        if not sample_id in records:
            records[sample_id] = {}

        records[sample_id][plasmid_id] = {
            'size': int(row['size']),
            'rep_type(s)':row['rep_type(s)'],
            'relaxase_type(s)':row['relaxase_type(s)'],
            'predicted_mobility':row['predicted_mobility'],
            'primary_cluster_id':row['primary_cluster_id'],
            'secondary_cluster_id':row['secondary_cluster_id'],
            'predicted_host_range_overall_rank':row['predicted_host_range_overall_rank'],
            'predicted_host_range_overall_name':row['predicted_host_range_overall_name'],
            'serovar':'',
            'resistance_genes':{},
        }

    return records


def process_metadata(file):
    data = read_data_into_df(file)
    records = {}
    for index, row in data.iterrows():
        sample_id = row['sample_id']
        serovar = row['serovar']
        continent = row['continent']
        country = row['country']
        year = row['collection_year']
        primary_sample_category = row['primary_sample_category']
        secondary_sample_category = row['secondary_sample_category']
        associated_taxa = str(row['associated_taxa']).split(';')
        records[sample_id] = {
            'serovar':serovar,
            'continent':continent,
            'country':country,
            'year':year,
            'primary_sample_category':primary_sample_category,
            'secondary_sample_category':secondary_sample_category,
            'associated_taxa':associated_taxa
        }
    return records

def get_meta_counts(records,key):
    counts = {}
    for sample_id in records:
        data = records[sample_id][key]
        if not data in counts:
            counts[data] = 0
        counts[data]+=1
    return counts

def get_gene_counts(records):
    counts = {}
    for sample_id in records:
        for contig_id in records[sample_id]:
            for gene_id in records[sample_id][contig_id]:
                if not gene_id in counts:
                    counts[gene_id] = 0
                counts[gene_id]+=1
    return counts

def associate_genes_with_plasmids(abricate,contigs):
    counts = {}
    for sample_id in abricate:
        if not sample_id in contigs:
            #print("skip...{}...not in contigs".format(sample_id))
            continue

        for contig_id in abricate[sample_id]:
            if not 'chromosome' in contigs[sample_id]:
                continue

            if contig_id in contigs[sample_id]['chromosome']:
                molecule_type = 'chromosome'
            else:
                molecule_type = 'plasmid'

            if not molecule_type  in contigs[sample_id]:
                continue

            for gene_id in abricate[sample_id][contig_id]:
                if not gene_id in counts:
                    counts[gene_id] = {
                        'chromosome':0,
                        'plasmid':{}
                    }
                if molecule_type == 'chromosome':
                    counts[gene_id]['chromosome']+=1
                else:
                    if gene_id == 'tet(G)':
                        print("====>{}".format(sample_id))
                        print(contigs[sample_id][molecule_type][contig_id])

                    clust_id = contigs[sample_id][molecule_type][contig_id]['primary_cluster_id']
                    if not clust_id in counts[gene_id]['plasmid']:
                        counts[gene_id]['plasmid'][clust_id] = 0
                    counts[gene_id]['plasmid'][clust_id] += 1
    return counts

def associate_serovars_with_plasmids(metadata,contigs):
    counts = {}
    for sample_id in contigs:
        if not sample_id in metadata:
            continue
        serovar = metadata[sample_id]['serovar']
        if not serovar in counts:
            counts[serovar] = {'total_samples':0,'plasmid_samples':0,'num_plasmids_per_sample':[],'plasmids':{}}
        counts[serovar]['total_samples']+=1
        if len(contigs[sample_id]['plasmid']) == 0:
            continue
        else:
            counts[serovar]['plasmid_samples']+=1
        found_plasmids = {}
        for contig_id in contigs[sample_id]['plasmid']:
            primary_cluster_id = contigs[sample_id]['plasmid'][contig_id]['primary_cluster_id']
            predicted_mobility = contigs[sample_id]['plasmid'][contig_id]['predicted_mobility']
            if not primary_cluster_id in found_plasmids and primary_cluster_id != '-':
                found_plasmids[primary_cluster_id] = predicted_mobility

        counts[serovar]['num_plasmids_per_sample'].append(len(found_plasmids))
        if primary_cluster_id not in counts[serovar]['plasmids']:
            counts[serovar]['plasmids'][primary_cluster_id] =0
        counts[serovar]['plasmids'][primary_cluster_id]+=1
    return counts

def build_plasmid_serovar_network(serovar_plasmid_associations):
    nodes = {}
    serovar_list = list(serovar_plasmid_associations.keys())
    for i in range(0,len(serovar_list)):
        serovar_1 = serovar_list[i]
        if not serovar_1 in nodes:
            nodes[serovar_1] = {}
        plasmids_1 = list(serovar_plasmid_associations[serovar_1]['plasmids'].keys())
        for k in range(i+1,len(serovar_list)):
            serovar_2 = serovar_list[k]
            if not serovar_2 in nodes[serovar_1]:
                nodes[serovar_1][serovar_2] = 0
            plasmids_2 = list(serovar_plasmid_associations[serovar_2]['plasmids'].keys())
            shared_plasmids = len(set(plasmids_1).intersection( set(plasmids_2)))
            if shared_plasmids > 0:
                nodes[serovar_1][serovar_2] += shared_plasmids

    edges = {}

    for s1 in nodes:
        for s2 in nodes[s1]:
            value = nodes[s1][s2]
            if value > 0:
                if not s1 in edges:
                    edges[s1] = {}
                edges[s1][s2] = value
    return edges

def build_plasmid_gene_network(gene_counts,gene_plasmid_associations):

    edges = {}
    gene_list = list(gene_counts.keys())
    #print(gene_plasmid_associations)
    for i in range(0, len(gene_list)):
        gene1 = gene_list[i]

        if not gene1 in gene_plasmid_associations:
            # print("skip {}".format(gene1))
            continue

        plasmids_1 = list(gene_plasmid_associations[gene1]['plasmid'].keys())
        if len(plasmids_1) <= 1:
            continue
       # print(plasmids_1)
        for k in range(0, len(plasmids_1) - 1):
            p1 = plasmids_1[k]
            p2 = plasmids_1[k + 1]

            if not p1 in edges:
                edges[p1] = {}

            if not p2 in edges[p1]:
                edges[p1][p2] = 0

            edges[p1][p2] += 1
    #print(edges)
    return edges


def build_gene_plasmid_network(gene_counts,gene_plasmid_associations):
    edges = {}
    gene_list = list(gene_counts.keys())
    for i in range(0,len( gene_list)):
        gene1 = gene_list[i]

        if not gene1 in gene_plasmid_associations:
           # print("skip {}".format(gene1))
            continue
        plasmids_1 = gene_plasmid_associations[gene1]['plasmid']
        if len(plasmids_1) == 0:
           # print("skip {} no plasmids".format(gene1))
            #print(gene_plasmid_associations[gene1])

            continue
        for k in range(i+1,len(gene_list)):
            gene2 = gene_list[k]
            if not gene2 in gene_plasmid_associations:
                continue
            plasmids_2 = gene_plasmid_associations[gene2]['plasmid']
            if len(plasmids_2) == 0:
                continue

            if not gene1 in edges:
                edges[gene1] = {}

            if not gene2 in edges[gene1]:
                edges[gene1][gene2] = 0

            edges[gene1][gene2] = len(set(plasmids_1).intersection(set(plasmids_2)))

    return edges



def write_plasmid_network_files(outdir,serovar_counts,edges,nodes_filename,edges_filename):
    edge_file = open(os.path.join(outdir,edges_filename),'w')
    nodes_file = open(os.path.join(outdir, nodes_filename), 'w')
    nodes_file.write("{}\t{}\n".format("id","size"))
    for serovar in serovar_counts:
        nodes_file.write("{}\t{}\n".format(serovar,serovar_counts[serovar]))
    nodes_file.close()

    edge_file.write("{}\t{}\t{}\n".format("Source", "Target","Weight"))
    for s1 in edges:
        for s2 in edges[s1]:
            edge_file.write("{}\t{}\t{}\n".format(s1,s2,edges[s1][s2]))

    edge_file.close()


def get_plasmid_counts(contigs):
    counts = {}
    for sample_id in contigs:
        plasmids = contigs[sample_id]['plasmid']
        if len(plasmids) == 0:
            continue
        clusters = {}
        for contig_id in plasmids:
            cluster_id = plasmids[contig_id]['primary_cluster_id']
            clusters[cluster_id] = ''
        clusters = list(clusters.keys())

        for cluster_id in clusters:
            if not cluster_id in counts:
                counts[cluster_id] = 0
            counts[cluster_id]+=1

    return counts


def get_gene_molecule_associations(gene_plasmid_associations):
    resistances = {}
    for gene_id in gene_plasmid_associations:
        if not gene_id in resistances:
            resistances[gene_id] = {'chromosome':0,'plasmid':0}
        if not 'chromosome' in gene_plasmid_associations[gene_id] or \
            not 'plasmid' in gene_plasmid_associations[gene_id]:
            continue
        resistances[gene_id]['chromosome'] = gene_plasmid_associations[gene_id]['chromosome']
        resistances[gene_id]['plasmid'] = sum(list(gene_plasmid_associations[gene_id]['plasmid'].values()))
    return resistances

def associate_genes_serovars(abricate,metadata):
    genes = {}
    for sample_id in abricate:
        if not sample_id in metadata:
            continue
        serovar = metadata[sample_id]['serovar']
        sample_source = metadata[sample_id]['primary_sample_category']
        if sample_source != 'human':
            sample_source = 'non-human'

        for contig_id in abricate[sample_id]:
            res_genes = abricate[sample_id][contig_id]
            for gene_id in res_genes:
                if not gene_id in genes:
                    genes[gene_id] = {'serovar':{},'sample_source':{'human':0,'non-human':0}}
                if not serovar in genes[gene_id]['serovar']:
                    genes[gene_id]['serovar'][serovar] = 0
                genes[gene_id]['serovar'][serovar]+=1

                genes[gene_id]['sample_source'][sample_source]+=1
    return genes

def add_serovar_to_mobtyper(mobtyper,metadata):
    for sample_id in mobtyper:
        if sample_id not in metadata:
            continue
        serovar = metadata[sample_id]['serovar']
        for plasmid_id in mobtyper[sample_id]:
            mobtyper[sample_id][plasmid_id]['serovar'] = serovar
    return mobtyper

def add_resistance_genes(abricate,mobtyper,contigs):
    for sample_id in abricate:
        if not sample_id in mobtyper or not sample_id in contigs:
            continue
        for contig_id in abricate[sample_id]:
            if len(contigs[sample_id]['plasmid']) == 0:
                break

            if not contig_id in contigs[sample_id]['plasmid']:
                continue

            plasmid_id = contigs[sample_id]['plasmid'][contig_id]['primary_cluster_id']
            for gene_id in abricate[sample_id][contig_id]:
                if not gene_id in mobtyper[sample_id][plasmid_id]['resistance_genes']:
                    mobtyper[sample_id][plasmid_id]['resistance_genes'][gene_id] = 0
                mobtyper[sample_id][plasmid_id]['resistance_genes'][gene_id]+=1
    return mobtyper

def add_metadata_to_mobtyper(mobtyper,metadata):
    for sample_id in mobtyper:

        if not sample_id in metadata:
            continue

        data = metadata[sample_id]
        for plasmid_id in mobtyper[sample_id]:
            mobtyper[sample_id][plasmid_id]['metadata'] = {
                'continent': 'unknown',
                'country': 'unknown',
                'serovar': 'unknown',
                'primary_sample_category': 'unknown',
                'secondary_sample_category': 'unknown',
                'associated_taxa': ['unknown'],
                'year': 'unknown'
            }
            mobtyper[sample_id][plasmid_id]['metadata']['continent'] = data['continent']
            mobtyper[sample_id][plasmid_id]['metadata']['country'] = data['country']
            mobtyper[sample_id][plasmid_id]['metadata']['serovar'] = data['serovar']
            mobtyper[sample_id][plasmid_id]['metadata']['primary_sample_category'] = data['primary_sample_category']
            mobtyper[sample_id][plasmid_id]['metadata']['secondary_sample_category'] = data['secondary_sample_category']
            mobtyper[sample_id][plasmid_id]['metadata']['associated_taxa'] = data['associated_taxa']
            mobtyper[sample_id][plasmid_id]['metadata']['year'] = data['year']
    return mobtyper

def mobtyper_plasmid_summarize(mobtyper):
    summary = {}
    for sample_id in mobtyper:
        plasmids = mobtyper[sample_id]
        for plasmid_id in plasmids:
            data = plasmids[plasmid_id]
            if not plasmid_id in summary:
                summary[plasmid_id] = {
                'replicons':{},
                'relaxases':{},
                'overall_mobility':'',
                'mobility': {'conjugative': 0, 'mobilizable': 0, 'non-mobilizable': 0},
                'overall_serovar':'',
                'serovar':{},
                'continent':{},
                'country':{},
                'primary_sample_category':{},
                'secondary_sample_category':{},
                'associated_taxa':{},
                'earliest_year':0,
                'year':{},
                'samples':[],
                'total_samples':0,
                'num_resistant': 0,
                'proportion_resistant':0,
                'resistance_genes': {},
                'serovar_entropy':-1,
                'serovar_shannon_index':-1,
                'serovar_simpson_index': -1,
                'serovar_simpson_index_e': -1,
                'serovar_chao1': 0,
                'num_serovars':0,
                'poportion_human':0,
                'taxa_entropy':-1,
                'taxa_shannon_index':-1,
                'taxa_simpson_index': -1,
                'taxa_simpson_index_e': -1,
                'taxa_chao1': -1,
            }
            summary[plasmid_id]['total_samples']+=1
            summary[plasmid_id]['samples'].append(sample_id)
            mobility = data['predicted_mobility']
            summary[plasmid_id]['mobility'][mobility] += 1

            rep = data['rep_type(s)'].split(",")
            for r in rep:
                if r not in summary[plasmid_id]['replicons']:
                    summary[plasmid_id]['replicons'][r] = 0
                summary[plasmid_id]['replicons'][r]+=1

            mob = data['relaxase_type(s)'].split(",")
            for m in mob:
                if m not in summary[plasmid_id]['relaxases']:
                    summary[plasmid_id]['relaxases'][m] = 0
                summary[plasmid_id]['relaxases'][m]+=1

            res_genes = data['resistance_genes']

            if len(res_genes) > 0:
                summary[plasmid_id]['num_resistant'] += 1
                for gene_id in res_genes:
                    if not gene_id in summary[plasmid_id]['resistance_genes']:
                        summary[plasmid_id]['resistance_genes'][gene_id] = 0
                    summary[plasmid_id]['resistance_genes'][gene_id] += res_genes[gene_id]


            if not 'metadata' in data:
                continue

            for field_id in data['metadata']:
                value = data['metadata'][field_id]

                if value == 'nan' or value == '':
                    value = 'unknown'

                if not field_id in summary[plasmid_id]:
                    continue

                if field_id == 'associated_taxa':
                    for v in value:
                        if v == '' or v == 'nan':
                            continue
                        if not v in summary[plasmid_id][field_id]:
                            summary[plasmid_id][field_id][v] = 0
                        summary[plasmid_id][field_id][v] += 1
                    continue
                if field_id in ('resistance_genes'):
                    continue

                if not value in summary[plasmid_id][field_id]:
                    summary[plasmid_id][field_id][value] = 0
                summary[plasmid_id][field_id][value] +=1

    for plasmid_id in summary:
        serovar_counts = list(summary[plasmid_id]['serovar'].values())
        if len(summary[plasmid_id]['year']) >0:
            summary[plasmid_id]['earliest_year'] = min(list(summary[plasmid_id]['year'].keys()))
        if 'human' in summary[plasmid_id]['primary_sample_category']:
            value = summary[plasmid_id]['primary_sample_category']['human']
        else:
            value = 0

        summary[plasmid_id]['poportion_human'] = value / summary[plasmid_id]['total_samples']

        summary[plasmid_id]['num_serovars'] = len(summary[plasmid_id]['serovar'])
        summary[plasmid_id]['proportion_resistant'] =  summary[plasmid_id]['num_resistant'] / summary[plasmid_id]['total_samples']

        summary[plasmid_id]['overall_mobility'] = max(summary[plasmid_id]['mobility'], key=summary[plasmid_id]['mobility'].get)
        if len(summary[plasmid_id]['serovar']) > 0:
            summary[plasmid_id]['overall_serovar'] = max(summary[plasmid_id]['serovar'],
                                                      key=summary[plasmid_id]['serovar'].get)

        if len(serovar_counts) > 0 and sum(serovar_counts) >= 10:
            summary[plasmid_id]['serovar_entropy'] = calc_shanon_entropy(serovar_counts)
            summary[plasmid_id]['serovar_shannon_index'] = alpha.shannon(serovar_counts)
            summary[plasmid_id]['serovar_simpson_index'] = alpha.simpson(serovar_counts)
            summary[plasmid_id]['serovar_simpson_index_e'] = alpha.simpson_e(serovar_counts)
            summary[plasmid_id]['serovar_chao1'] = alpha.chao1(serovar_counts)
        else:
            print("{}\t{}".format(plasmid_id,sum(serovar_counts)))
            print(summary[plasmid_id])
        human_removed_taxa = {}
        for taxon in summary[plasmid_id]['associated_taxa']:
            if taxon == 'homo sapiens':
                continue
            human_removed_taxa[taxon] = summary[plasmid_id]['associated_taxa'][taxon]

        taxa_counts = list(human_removed_taxa.values())
        if len(taxa_counts) > 0 and sum(taxa_counts) >= 10:
            summary[plasmid_id]['taxa_entropy'] = calc_shanon_entropy(taxa_counts)
            summary[plasmid_id]['taxa_shannon_index'] = alpha.shannon(taxa_counts)
            summary[plasmid_id]['taxa_simpson_index'] = alpha.simpson(taxa_counts)
            summary[plasmid_id]['taxa_simpson_index_e'] = alpha.simpson_e(taxa_counts)
            summary[plasmid_id]['taxa_chao1'] = alpha.chao1(taxa_counts)

    return summary

def serovar_metadata_summary(metadata):
    serovar_info = {}
    for sample_id in metadata:
        data = metadata[sample_id]
        serovar = data['serovar']
        if len(serovar) == 0:
            serovar = 'unknown'
        if not serovar in serovar_info:
            serovar_info[serovar] = {
            'continent':{},
            'country':{},
            'primary_sample_category':{},
            'secondary_sample_category':{},
            'associated_taxa':{},
            'year':{},
            'samples':[],
            'total_samples':0,
            'plasmids':{},
            'num_resistant': 0,
            'proportion_resistant': 0,
            'resistance_genes': {},
            'poportion_human': 0,
            'taxa_shannon_index': -1,
            'taxa_simpson_index': -1,
            'taxa_simpson_index_e': -1,
            'taxa_chao1': -1,
            'taxa_entropy':-1,
            'plasmid_entropy':-1,
            'plasmid_shannon_index': -1,
            'plasmid_simpson_index': -1,
            'plasmid_simpson_index_e': -1,
            'plasmid_chao1': -1,

        }
        serovar_info[serovar]['samples'].append(sample_id)
        serovar_info[serovar]['total_samples']+=1


        for field_id in data:
            value = data[field_id]
            if value == 'nan' or value == '':
                value = 'unknown'
            if not field_id in serovar_info[serovar]:
                continue

            if field_id == 'associated_taxa':
                for v in value:
                    if not v in serovar_info[serovar][field_id]:
                        serovar_info[serovar][field_id][v] = 0
                    serovar_info[serovar][field_id][v] += 1
                continue

            if not value in serovar_info[serovar][field_id]:
                serovar_info[serovar][field_id][value] = 0
            serovar_info[serovar][field_id][value] +=1

    for serovar in serovar_info:
        data = serovar_info[serovar]
        if 'human' in data['primary_sample_category']:
            value = data['primary_sample_category']['human']
        else:
            value = 0

        serovar_info[serovar]['poportion_human'] = value / serovar_info[serovar]['total_samples']

    return serovar_info

def add_resgenes_to_serovar_summary(abricate,serovar_data):

    for serovar in serovar_data:
        samples = serovar_data[serovar]['samples']
        for sample_id in samples:
            if not sample_id in abricate:
                continue
            for contig_id in abricate[sample_id]:
                for gene_id in abricate[sample_id][contig_id]:
                    if not gene_id in serovar_data[serovar]['resistance_genes']:
                        serovar_data[serovar]['resistance_genes'][gene_id] = 0
                    serovar_data[serovar]['resistance_genes'][gene_id] += 1
            if len(serovar_data[serovar]['resistance_genes']) > 0:
                serovar_data[serovar]['num_resistant']+=1
        serovar_data[serovar]['proportion_resistant'] = serovar_data[serovar]['num_resistant']/serovar_data[serovar]['total_samples']

    return serovar_data

def add_plasmids_to_serovar(mobtyper,serovar_data):
    for serovar in serovar_data:
        samples = serovar_data[serovar]['samples']
        for sample_id in samples:
            if not sample_id in mobtyper:
                continue
            plasmids = mobtyper[sample_id]
            for plasmid_id in plasmids:
                if not plasmid_id in serovar_data[serovar]['plasmids']:
                    serovar_data[serovar]['plasmids'][plasmid_id] = 0
                serovar_data[serovar]['plasmids'][plasmid_id]+=1
    return serovar_data

def add_diversity_data_to_serovar(serovar_data):
    for serovar in serovar_data:
        human_removed_taxa = {}
        data = serovar_data[serovar]
        for taxon in data['associated_taxa']:
            if taxon == 'homo sapien':
                continue
            human_removed_taxa[taxon] = data['associated_taxa'][taxon]

        taxa_counts = list(human_removed_taxa.values())
        if len(taxa_counts) > 0 and sum(taxa_counts) >= 10:
            serovar_data[serovar]['taxa_entropy'] = calc_shanon_entropy(taxa_counts)
            serovar_data[serovar]['taxa_shannon_index'] = alpha.shannon(taxa_counts)
            serovar_data[serovar]['taxa_simpson_index'] = alpha.simpson(taxa_counts)
            serovar_data[serovar]['taxa_simpson_index_e'] = alpha.simpson_e(taxa_counts)
            serovar_data[serovar]['taxa_chao1'] = alpha.chao1(taxa_counts)

        plasmid_counts = list(serovar_data[serovar]['plasmids'].values())
        if len(plasmid_counts) > 0 and sum(plasmid_counts) >= 10:
            serovar_data[serovar]['plasmid_entropy'] = calc_shanon_entropy(plasmid_counts)
            serovar_data[serovar]['plasmid_shannon_index'] = alpha.shannon(plasmid_counts)
            serovar_data[serovar]['plasmid_simpson_index'] = alpha.simpson(plasmid_counts)
            serovar_data[serovar]['plasmid_simpson_index_e'] = alpha.simpson_e(plasmid_counts)
            serovar_data[serovar]['plasmid_chao1'] = alpha.chao1(plasmid_counts)

    return  serovar_data

def write_gene_results(gene_plasmid_associations,molecule_type_gene_association,gene_serovar_associations):
    print(">>>>>>>>>>>>>>>>>>>>>")
    print("gene_id\tplasmid\tchromosome\ttotal\t"
          "proportion\tplasmid homogeneity\tserovar homogeneity\t"
          "plasmid entropy\tserovar entropy\thuman proportion")

    for gene_id in molecule_type_gene_association:
        plasmid_count = int(molecule_type_gene_association[gene_id]['plasmid'])
        chr_count = int(molecule_type_gene_association[gene_id]['chromosome'])
        total = plasmid_count + chr_count

        if len(gene_plasmid_associations[gene_id]['plasmid']) > 0:
            phomogeneity = calc_homogeneity(list(gene_plasmid_associations[gene_id]['plasmid'].values()))
            pentropy = calc_shanon_entropy(list(gene_plasmid_associations[gene_id]['plasmid'].values()))
        else:
            phomogeneity = -1
            pentropy = -1
        if gene_id in gene_serovar_associations and len(gene_serovar_associations[gene_id]['serovar']) > 0:
            shomogeneity = calc_homogeneity(list(gene_serovar_associations[gene_id]['serovar'].values()))
            human_proportion = gene_serovar_associations[gene_id]['sample_source']['human'] / \
                               (gene_serovar_associations[gene_id]['sample_source']['human'] +
                                gene_serovar_associations[gene_id]['sample_source']['non-human'])
            sentropy = calc_shanon_entropy(list(gene_serovar_associations[gene_id]['serovar'].values()))
        else:
            shomogeneity = -1
            human_proportion = -1
            sentropy = -1

        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(gene_id,
                                                              plasmid_count,
                                                              chr_count,
                                                              total,
                                                              float(plasmid_count) / total,
                                                              phomogeneity,
                                                              shomogeneity,
                                                              pentropy,
                                                              sentropy,
                                                              human_proportion))

def write_dict(data):

    header = ['key']
    for key in data:
        fields = data[key].keys()
        for f in fields:
            if f == 'samples' or f == 'plasmids':
                continue
            header.append(f)
        break
    print(">>>>>>>>>>>>>>>>>>>>>>>>>")
    print("{}".format("\t".join(header)))
    for key in data:
        data[key]['key'] = key


        row = []
        for f in header:
            if f in data[key]:

                row.append("{}".format(data[key][f]))
            else:
                row.append('nan')
        print("{}".format("\t".join(row)))
        #if key == 'Rubislaw':
         #   print (data[key])
          #  print(row)
           # sys.exit()

def calc_pairwise_jacard(data,field_id):
    sample_keys = list(data.keys())
    num_samples = len(sample_keys)
    distances = {}
    for i in range(0,num_samples):
        sample_id_1 = sample_keys[i]
        samples_1 = data[sample_id_1][field_id]
        if not sample_id_1 in distances:
            distances[sample_id_1] = {}
        for k in range(i,num_samples):
            sample_id_2 = sample_keys[k]
            samples_2 = data[sample_id_2][field_id]
            sample_union = set(samples_1).union(set(samples_2) )
            sample_intersection = set(samples_1).intersection(set(samples_2) )
            if len(sample_union) > 0:
                jaccard = 1 - len(sample_intersection)/len(sample_union)
            else:
                jaccard = 1
            distances[sample_id_1][sample_id_2] = jaccard
    return distances

def write_distance_matrix(data):
    sample_keys = list(data.keys())
    num_samples = len(sample_keys)
    print("\t{}".format("\t".join(sample_keys)))
    for i in range(0,num_samples):
        sample_id_1 = sample_keys[i]
        row = [""] * (i+1)
        for k in range(i,num_samples):
            sample_id_2 = sample_keys[k]
            distance = data[sample_id_1][sample_id_2]
            row.append(distance)
        print("{}{}".format(sample_id_1,'\t'.join(str(x) for x in row)))




def add_plasmid_proportion(serovar_data,mobtyper):
    for serovar in serovar_data:
        serovar_data[serovar]['count_plasmid_positive_samples'] = 0
    for sample_id in mobtyper:
        for serovar in serovar_data:
            if sample_id in serovar_data[serovar]['samples']:
                serovar_data[serovar]['count_plasmid_positive_samples']+=1
                break
    return serovar_data

def main():
    args = parse_args()
    mobtyper = process_mobtyper(args.mobtyper)
    metadata = process_metadata(args.metadata)
    mobtyper = add_metadata_to_mobtyper(mobtyper, metadata)


    outdir = args.outdir

    serovar_data = serovar_metadata_summary(metadata)

    abricate = process_abricate(args.abricate)

    contigs = process_contigs(args.contigs)

    mobtyper = add_resistance_genes(abricate, mobtyper, contigs)

    plasmid_summary = mobtyper_plasmid_summarize(mobtyper)

    serovar_data = add_diversity_data_to_serovar(add_plasmids_to_serovar(mobtyper,add_resgenes_to_serovar_summary(abricate,serovar_data)))
    serovar_data = add_plasmid_proportion(serovar_data,mobtyper)
    write_dict(plasmid_summary)
    write_dict(serovar_data)
    write_distance_matrix(calc_pairwise_jacard(plasmid_summary,'samples'))
    write_distance_matrix(calc_pairwise_jacard(serovar_data, 'plasmids'))


    gene_plasmid_associations = associate_genes_with_plasmids(abricate,contigs)
    molecule_type_gene_association = get_gene_molecule_associations(gene_plasmid_associations)
    gene_serovar_associations = associate_genes_serovars(abricate,metadata)

    write_gene_results(gene_plasmid_associations, molecule_type_gene_association, gene_serovar_associations)


    sys.exit()







main()

