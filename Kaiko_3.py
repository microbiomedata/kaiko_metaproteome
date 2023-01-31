## This is originally ExtractTopRankTaxaFromUniRef.py

import pandas as pd
import numpy as np
import time


# @profile    
def run_diamond_tally(diamond_output, ntops, ncbi_taxa_folder, mode, fout, pident):

    filterby={}
    if pident: filterby['pident'] = float(pident)
    # #     if FLAGS.evalue: filterby['evalue'] = FLAGS.evalue
    # #     if FLAGS.mismatch: filterby['mismatch'] = FLAGS.mismatch
    print('filterby:', filterby)

    ############################################################
    # retrieve taxa lineage info such as taxa rank and lineage
    ############################################################
    print("Retrieving taxa lineage info...")
    taxa_table = read_ncbi_taxa_lineage(ncbi_taxa_folder / 'rankedlineage.dmp', 
                                        ncbi_taxa_folder / 'nodes.dmp')
    
    ############################################################
    # read diamond output file
    ############################################################
    print("Reading diamond output file...")
    dmd = read_dmd(diamond_output)

    ############################################################
    # filter by quality and taxa
    ############################################################
    print("Filtering by quality and taxa...")
    filtered_dmd = dmd_filter(dmd, filterby=filterby)
    filtered_dmd = collect_taxid(filtered_dmd)
    filtered_dmd['uniref_id'] = [i[0] for i in filtered_dmd.uniref_seq.str.split(" ",1)]

    ############################################################
    # garbage collection
    ############################################################
    print("Garbage collection...")
    del dmd
    import gc
    gc.collect()
    
    ############################################################
    # retrieve UniRef100 representative taxa and its members
    ############################################################
    print("Retrieving UniRef100 representative taxa and members...")
    unique_unirefs = set(filtered_dmd.uniref_id.drop_duplicates())
    taxa_member_tbl = get_taxa_members(ncbi_taxa_folder / 'uniref100_member_taxa_tbl.csv',
                                       unique_unirefs)
    
    ############################################################
    # merge
    ############################################################
    print("Adding taxa members...")
    merged = filtered_dmd.merge(taxa_member_tbl, left_on="uniref_id", right_on="uid", how="inner")
    print("  Final table size:{}".format(merged.shape))
    
    if merged.shape[0]!=filtered_dmd.shape[0]:
        print("[WARN] You might use the different version of UniRef100.fasta and .xml")
    
    # get the member taxa of each (scan, uid)
    unique_members = []
    scanids = merged.scans.tolist()
    uids = merged.uid.tolist()
    commons = merged.common_taxa.tolist()
    for ii, members in enumerate(merged.members.str.split(":").tolist()):
        for mm in members:
            unique_members.append({"scan":scanids[ii], "uid":uids[ii], "member_taxa":int(mm), "common_taxa":int(commons[ii])})
    print("  #members:{}".format(len(unique_members)))
    members_df = pd.DataFrame(unique_members)

    ############################################################
    # top-rank taxa
    ############################################################
    print("Filtering top-rank taxa by hits...")
    
    if mode=="member":
        taxa_col = 'member_taxa'
    elif mode=="common":
        taxa_col = 'common_taxa'
    else:
        print("  [ERR] please select the mode to be either 'member' or 'common'.")
        return
    
    unique_pepseq_taxa = members_df.drop_duplicates(subset=['scan',taxa_col])
    pepcount_taxid = unique_pepseq_taxa[taxa_col].value_counts()

    print('  unique peptides: {}'.format(members_df.scan.value_counts().shape[0]))
    print('  unique taxa: {}'.format(pepcount_taxid.shape[0]))

    def besthit(row):
        return pepcount_taxid[row[taxa_col]].nlargest(1, keep='all').index.tolist()

    besthits = []
    for besthit in unique_pepseq_taxa.groupby(by='scan').apply(besthit):
        besthits += besthit
    besthits = pd.Series(besthits).value_counts()
    
    if ntops > 0:
        _ntops = ntops
    else:
        _ntops = besthits.shape[0]
    
    top_taxa = besthits.nlargest(_ntops, keep='all').to_frame(name='hits')

    ############################################################
    # save top-rank taxa info
    ############################################################
    print("Saving top-rank taxa info...")
    df = pd.concat([taxa_table, top_taxa], join='inner', axis=1)
    df.index.name = 'tax_id'
    if fout: df.sort_values(by="hits", ascending=False).to_csv(fout)

# taxon info
def read_ncbi_taxa_lineage(rankedlineage_file, nodes_file):
    taxon_info = pd.read_csv(rankedlineage_file, sep='|', header=None,
                             names=['tax_id','tax_name','species','genus','family','order',
                                    'class','phylum','kingdom','superkingdom','null'])
    taxon_info = taxon_info.replace('\t\t', np.nan)
    for col in ['tax_name','species','genus','family','order','class','phylum','kingdom','superkingdom']:
        taxon_info[col] = taxon_info[col].str.extract(r'^\t(.*)\t$')[0].tolist()
        # taxon_info[col] = taxon_info[col].str.extract(r'^\t(.*)\t$').tolist()
    del taxon_info['null']
    print('  taxon_info:', taxon_info.shape, taxon_info.tax_id.drop_duplicates().shape)
    
    rank = pd.read_csv(nodes_file, sep='|', usecols=[0,2], header=None,
                     names=['tax_id','rank'])
    rank = rank.replace('\t\t', np.nan)
    rank['rank'] = rank['rank'].str.extract(r'^\t(.*)\t$')[0].tolist()
    # rank['rank'] = rank['rank'].str.extract(r'^\t(.*)\t$').tolist()
    print('  rank:', rank.shape, rank.tax_id.drop_duplicates().shape)
    
    final_df = taxon_info.merge(rank, left_on='tax_id', right_on='tax_id')
    print('  final_df:', final_df.shape, final_df.tax_id.drop_duplicates().shape)
    return final_df.set_index('tax_id')

def read_dmd(diamond_output):
    dmd_colnames = ['scans','uniref_seq','pident','evalue','mismatch']
    dmd = pd.read_csv(diamond_output, sep='\t', header=None, names=dmd_colnames)
    return dmd

def dmd_filter(dmd, filterby={}):
    if len(filterby) > 0:
        filterby_cond = None
        for col in filterby:
            if col == 'pident':
                tmp = dmd[col]>=filterby[col]
            elif col in ['evalue', 'mismatch']:
                tmp = dmd[col]<=filterby[col]
            else:
                print("  [ERR] no column found. For filtering, please use 'pident','evalue','mismatch'.")
                raise Exception

            if filterby_cond is None:
                filterby_cond = tmp
            else:
                filterby_cond &= tmp

        filtered_dmd = dmd[filterby_cond].copy()
        print("  org:{}, filtered:{}".format(dmd.shape[0], filtered_dmd.shape[0]))
        return filtered_dmd
    else:
        return dmd

def collect_taxid(filtered_dmd):
    filtered_dmd['taxid'] = filtered_dmd.uniref_seq.str.extract('^.+? TaxID=(\d*)\s?')
    # filtered_dmd['taxid'] = filtered_dmd.uniref_seq.str.extract('^.+? OX=(\d*)\s?')

    # drop the irrelevant
    # unknown
    # 1  # root
    # 9606  # Homo sapiens
    # 412755  # marine sediment metagenome
    # 408172  # marine metagenome
    drop_taxids = set(['1','2','9606','412755','408172'])

    _filtered_dmd = filtered_dmd[~((filtered_dmd.taxid=='')|(filtered_dmd.taxid.isin(drop_taxids)))].copy()
    _filtered_dmd.taxid = _filtered_dmd.taxid.astype(int)
    _filtered_dmd = _filtered_dmd.reset_index(drop=True)
    print('filter relevant taxa:{}'.format(_filtered_dmd.shape[0]))

    wrong_taxids = [444888, 55087, 210425]
    print("{} wrong taxids found in the results and corrected.".format(
        _filtered_dmd[_filtered_dmd.taxid.isin(wrong_taxids)].shape[0]))

    _filtered_dmd.loc[_filtered_dmd.taxid==444888, 'taxid'] = 629  # Yersinia
    _filtered_dmd.loc[_filtered_dmd.taxid==55087, 'taxid'] = 1386  # Bacillus
    _filtered_dmd.loc[_filtered_dmd.taxid==210425, 'taxid'] = 583  # Proteus
    return _filtered_dmd

def get_taxa_members(member_tbl_file, unique_unirefs):
    chunksize = 1000000

    stime = time.time()
    dfs = []
    num_iters = 0
    for chunk in pd.read_csv(member_tbl_file, chunksize=chunksize):
        tdf = chunk[chunk.uid.isin(unique_unirefs)]
        dfs.append(tdf)
        num_iters += 1

    df = pd.concat(dfs)
    print("  #Chunk:{}, Size:{}, {:.2f}min".format(len(dfs), df.shape, (time.time()-stime)/60))
    return df


