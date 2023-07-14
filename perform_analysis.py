#!/usr/bin/env python
import os
import sys
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score,adjusted_rand_score
from scipy.stats import entropy
from skbio.diversity.alpha import simpson
import scipy.stats as stats
import skbio.diversity.alpha as alpha



def calc_MI(category_1,category_2):
    return adjusted_mutual_info_score(category_1,category_2,average_method='arithmetic')

def calc_fisherExact(contingency_table):
    return stats.fisher_exact(contingency_table)

def calc_shanon_entropy(value_list ,normalize=False):
    total = sum(value_list)
    values = []
    for v in value_list:
        values.append( v /total)
    e = entropy(values)
    if normalize:
        values = [1] * len(value_list)
        e = e / (entropy(values))

    return e


def calc_homogeneity(value_list):
    return 1 - simpson(value_list)

def calc_adjusted_rand(category_1,category_2):
    return adjusted_rand_score(category_1,category_2)

def parse_args():
    "Parse the input arguments, use '-h' for help"
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass
    parser = ArgumentParser(
        description="Process result files from different tools for visualization of plasmid dynamic",
        formatter_class=CustomFormatter)
    parser.add_argument('--abricate', type=str, required=True, help='Abricate report')
    parser.add_argument('--contigs', type=str, required=True, help='MOB-suite contig report')
    parser.add_argument('--mobtyper', type=str, required=True, help='MOB-suite mobtyper report')
    parser.add_argument('--metadata', type=str, required=True, help='Sample metadata')
    parser.add_argument('--gene_id', type=str, required=False, help='Sample metadata',default='blaCMY-2')
    parser.add_argument('--outdir', type=str, required=True, help='Output results to this directory')
    return parser.parse_args()


def read_data_into_df(file):
    return pd.read_csv(file,header=0, encoding = "UTF-8",sep="\t",low_memory=False)

def process_abricate(file,feat_id='blaCMY-2'):
    data = read_data_into_df(file)
    records = {}
    mapping = {}
    gene_samples = {}
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
        mapping[gene_id] = resistance
        if gene_id == feat_id:
            if not sample_id in gene_samples:
                gene_samples[sample_id] = []
            gene_samples[sample_id].append(contig_id)
    return {'records':records,'mapping':mapping,'target_samples':gene_samples}

def process_contigs(file):
    data = read_data_into_df(file)
    records = {}
    for index, row in data.iterrows():
        sample_id = row['sample_id']
        def_line = row['contig_id'].split("_")
        if len(def_line) == 1:
            def_line = row['contig_id'].split(" ")
        contig_id = str(def_line[0])

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
    if row['mash_neighbor_distance'] > 0.025:
        records[sample_id][plasmid_id]['secondary_cluster_id'] = '-'

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
            if not contig_id in contigs[sample_id][molecule_type]:
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
            if not plasmid_id in mobtyper[sample_id]:
                continue
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
                'plasmid_id':plasmid_id,
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

                if field_id in ['associated_taxa','resistance_genes']:
                    for v in value:
                        if v == '' or v == 'nan':
                            continue
                        if not v in summary[plasmid_id][field_id]:
                            summary[plasmid_id][field_id][v] = 0
                        summary[plasmid_id][field_id][v] += 1
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
            summary[plasmid_id]['serovar_entropy'] = calc_shanon_entropy(serovar_counts,True)
            summary[plasmid_id]['serovar_shannon_index'] = alpha.shannon(serovar_counts)
            summary[plasmid_id]['serovar_simpson_index'] = alpha.simpson(serovar_counts)
            summary[plasmid_id]['serovar_simpson_index_e'] = alpha.simpson_e(serovar_counts)
            summary[plasmid_id]['serovar_chao1'] = alpha.chao1(serovar_counts)

        human_removed_taxa = {}
        for taxon in summary[plasmid_id]['associated_taxa']:
            if taxon == 'homo sapiens':
                continue
            human_removed_taxa[taxon] = summary[plasmid_id]['associated_taxa'][taxon]

        taxa_counts = list(human_removed_taxa.values())
        if len(taxa_counts) > 0 and sum(taxa_counts) >= 10:
            summary[plasmid_id]['taxa_entropy'] = calc_shanon_entropy(taxa_counts,True)
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
            'serovar':serovar,
            'continent':{},
            'country':{},
            'primary_sample_category':{},
            'secondary_sample_category':{},
            'associated_taxa':{},
            'year':{},
            'samples':[],
            'total_samples':0,
            'count_plasmid_negative_samples':0,
            'count_plasmid_positive_samples':0,
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
            if field_id == 'serovar':
                continue
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

def add_plasmid_freq(serovar_data,mobtyper,metadata):
    for sample_id in mobtyper:
        if not sample_id in metadata:
            continue
        serovar = metadata[sample_id]['serovar']
        if serovar in serovar_data:
            serovar_data[serovar]['count_plasmid_positive_samples']+=1
    for serovar in serovar_data:
        serovar_data[serovar]['count_plasmid_negative_samples'] = serovar_data[serovar]['total_samples'] - serovar_data[serovar]['count_plasmid_positive_samples']
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
            serovar_data[serovar]['taxa_entropy'] = calc_shanon_entropy(taxa_counts, True)
            serovar_data[serovar]['taxa_shannon_index'] = alpha.shannon(taxa_counts)
            serovar_data[serovar]['taxa_simpson_index'] = alpha.simpson(taxa_counts)
            serovar_data[serovar]['taxa_simpson_index_e'] = alpha.simpson_e(taxa_counts)
            serovar_data[serovar]['taxa_chao1'] = alpha.chao1(taxa_counts)

        plasmid_counts = list(serovar_data[serovar]['plasmids'].values())
        if len(plasmid_counts) > 0 and sum(plasmid_counts) >= 10:
            serovar_data[serovar]['plasmid_entropy'] = calc_shanon_entropy(plasmid_counts, True)
            serovar_data[serovar]['plasmid_shannon_index'] = alpha.shannon(plasmid_counts)
            serovar_data[serovar]['plasmid_simpson_index'] = alpha.simpson(plasmid_counts)
            serovar_data[serovar]['plasmid_simpson_index_e'] = alpha.simpson_e(plasmid_counts)
            serovar_data[serovar]['plasmid_chao1'] = alpha.chao1(plasmid_counts)

    return  serovar_data

def write_gene_results(gene_plasmid_associations,molecule_type_gene_association,gene_serovar_associations,mapping,outfile):
    fh = open(outfile,'w')
    fh.write("gene_id\tplasmid\tchromosome\ttotal\t"
          "proportion\tplasmid homogeneity\tserovar homogeneity\t"
          "plasmid entropy\tserovar entropy\thuman proportion\tresistance\n")

    for gene_id in molecule_type_gene_association:
        resistance = mapping[gene_id]
        plasmid_count = int(molecule_type_gene_association[gene_id]['plasmid'])
        chr_count = int(molecule_type_gene_association[gene_id]['chromosome'])
        total = plasmid_count + chr_count

        if len(gene_plasmid_associations[gene_id]['plasmid']) > 0:
            phomogeneity = calc_homogeneity(list(gene_plasmid_associations[gene_id]['plasmid'].values()))
            pentropy = calc_shanon_entropy(list(gene_plasmid_associations[gene_id]['plasmid'].values()),True)
        else:
            phomogeneity = -1
            pentropy = -1
        if gene_id in gene_serovar_associations and len(gene_serovar_associations[gene_id]['serovar']) > 0:
            shomogeneity = calc_homogeneity(list(gene_serovar_associations[gene_id]['serovar'].values()))
            human_proportion = gene_serovar_associations[gene_id]['sample_source']['human'] / \
                               (gene_serovar_associations[gene_id]['sample_source']['human'] +
                                gene_serovar_associations[gene_id]['sample_source']['non-human'])
            sentropy = calc_shanon_entropy(list(gene_serovar_associations[gene_id]['serovar'].values()),True)
        else:
            shomogeneity = -1
            human_proportion = -1
            sentropy = -1

        fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene_id,
                                                              plasmid_count,
                                                              chr_count,
                                                              total,
                                                              float(plasmid_count) / total,
                                                              phomogeneity,
                                                              shomogeneity,
                                                              pentropy,
                                                              sentropy,
                                                              human_proportion,resistance))
    fh.close()

def write_dict(data,outfile):
    df = pd.DataFrame.from_dict(data,orient='index')
    df.to_csv(outfile,sep="\t",header=True,index=False)

def add_plasmid_proportion(serovar_data,mobtyper):
    for serovar in serovar_data:
        serovar_data[serovar]['count_plasmid_positive_samples'] = 0
    for sample_id in mobtyper:
        for serovar in serovar_data:
            if sample_id in serovar_data[serovar]['samples']:
                serovar_data[serovar]['count_plasmid_positive_samples']+=1
                break
    return serovar_data

def create_specific_gene_info(contigs,metadata,target_samples,outfile):
    data = {}
    for sample_id in target_samples:
        if not sample_id in contigs or not sample_id in metadata:
            continue
        target_contigs = target_samples[sample_id]
        serovar = metadata[sample_id]['serovar']
        if not serovar in data:
            data[serovar] = {
                'chromosome':0,
                'plasmid':{}
            }

        for contig_id in target_contigs:
            if not 'chromosome' in contigs[sample_id]:
                continue
            if contig_id in contigs[sample_id]['chromosome']:
                molecule_type = 'chromosome'
            else:
                molecule_type = 'plasmid'
            if molecule_type == 'plasmid':
                if contig_id not in contigs[sample_id][molecule_type]:
                    continue
                clust_id = contigs[sample_id][molecule_type][contig_id]['primary_cluster_id']
                sec_clust_id = contigs[sample_id][molecule_type][contig_id]['secondary_cluster_id']
                if not clust_id in data[serovar]['plasmid']:
                    data[serovar]['plasmid'][clust_id]= {}
                if not sec_clust_id in data[serovar]['plasmid'][clust_id]:
                    data[serovar]['plasmid'][clust_id][sec_clust_id] = 0
                data[serovar]['plasmid'][clust_id][sec_clust_id] += 1
            else:
                data[serovar]['chromosome']+=1
    fh = open(outfile,'w')
    fh.write("serovar\tmolecule_type\tprimary_id\tsecondary_id\tcount\n")
    for serovar in data:
        fh.write("{}\t{}\t{}\t{}\t{}\n".format(serovar,'chromosome','-','-',data[serovar]['chromosome']))
        for clust_id in data[serovar]['plasmid']:
            for secondary_clust_id in data[serovar]['plasmid'][clust_id]:
                fh.write("{}\t{}\t{}\t{}\t{}\n".format(serovar, 'plasmid', clust_id, secondary_clust_id, data[serovar]['plasmid'][clust_id][secondary_clust_id]))
    fh.close()
    return data

def write_replicon_relaxase(mobtyper,outfile):
    replicons = set()
    relaxases = set()
    for sample_id in mobtyper:
        plasmids = mobtyper[sample_id]
        for plasmid_id in plasmids:
            data = plasmids[plasmid_id]
            rep = set(data['rep_type(s)'].split(","))
            replicons = replicons | rep
            mob = set(data['relaxase_type(s)'].split(","))
            relaxases = relaxases | mob
    replicons = list(replicons)
    relxases = sorted(list(relaxases))
    counts = {}
    for r in replicons:
        counts[r] = {}
        for m in relaxases:
            counts[r][m] = 0
    for sample_id in mobtyper:
        plasmids = mobtyper[sample_id]
        for plasmid_id in plasmids:
            data = plasmids[plasmid_id]
            rep = data['rep_type(s)'].split(",")
            mob = data['relaxase_type(s)'].split(",")
            for r in rep:
                for m in mob:
                    counts[r][m]+=1

    fh = open(outfile,'w')

    fh.write("replicon\t{}\ttotal\n".format("\t".join(relaxases)))
    for replicon in counts:
        row = [replicon]
        total = 0
        for m in relaxases:
            row.append(counts[replicon][m])
            total += counts[replicon][m]
        row.append(total)
        fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()




def main():
    args = parse_args()

    outdir = args.outdir
    # initialize analysis directory
    if not os.path.isdir(args.outdir):
        print("Creating output directory {}".format(args.outdir))
        os.makedirs(args.outdir, 0o755)

    print("Reading sample metadata file {}".format(args.metadata))
    metadata = process_metadata(args.metadata)
    serovar_data = serovar_metadata_summary(metadata)

    print("Reading MOB-typer file {}".format(args.mobtyper))
    mobtyper = process_mobtyper(args.mobtyper)

    print("Merging MOB-typer with Sample metadata")
    mobtyper = add_metadata_to_mobtyper(mobtyper, metadata)

    print("Adding plasmid counts by serovar")
    serovar_data = add_plasmid_freq(serovar_data,mobtyper,metadata)

    print("Reading abricate file {}".format(args.abricate))
    abricate_results = process_abricate(args.abricate)
    abricate = abricate_results['records']
    mapping = abricate_results['mapping']


    print("Reading MOB-recon contigs file {}".format(args.contigs))
    contigs = process_contigs(args.contigs)

    print("Merging plasmid and abricate results")
    mobtyper = add_resistance_genes(abricate, mobtyper, contigs)
    plasmid_summary = mobtyper_plasmid_summarize(mobtyper)


    create_specific_gene_info(contigs, metadata, abricate_results['target_samples'], os.path.join(outdir,"specific.feature.data.txt"))

    print("Performing diversity estimates")
    serovar_data = add_diversity_data_to_serovar(add_plasmids_to_serovar(mobtyper,add_resgenes_to_serovar_summary(abricate,serovar_data)))
    serovar_data = add_plasmid_proportion(serovar_data,mobtyper)
    write_dict(plasmid_summary,os.path.join(outdir,"plasmid.summary.txt"))
    write_dict(serovar_data,os.path.join(outdir,"serovar.summary.txt"))
    gene_plasmid_associations = associate_genes_with_plasmids(abricate,contigs)
    molecule_type_gene_association = get_gene_molecule_associations(gene_plasmid_associations)
    gene_serovar_associations = associate_genes_serovars(abricate,metadata)
    write_gene_results(gene_plasmid_associations, molecule_type_gene_association, gene_serovar_associations,mapping,os.path.join(outdir,"genes.summary.txt"))

    write_replicon_relaxase(mobtyper, os.path.join(outdir,"replicon.table.txt"))

if __name__== '__main__':
    main()

