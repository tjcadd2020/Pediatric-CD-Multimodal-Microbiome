import pandas as pd
import numpy as np
from io import StringIO
import os
import re

#SV
###SGVFinder自带的reference对应的物种详情
species_reference = pd.read_csv('CD/SV/representatives.genomes.taxonomy.csv')
species_reference ['Species'] = species_reference .apply(lambda row: row['organism'] if pd.isna(row['Species']) else row['Species'], axis=1)
species_reference

###SGVFinder自带的reference对应的物种的基因组详情
gene_reference = open('CD/SV/representatives.genes.drepped.annotations.df','rb')
gene_reference=pd.read_pickle(gene_reference)
gene_reference[['dest', 'gene_name']] = gene_reference.index.tolist()
gene_reference[['Taxid', 'project_id']] = gene_reference['dest'].str.split('.', expand=True)
gene_reference['start_pos'] = gene_reference['start_pos']/1000
gene_reference['end_pos'] = gene_reference['end_pos']/1000
gene_reference = gene_reference.reset_index(drop=True)
gene_reference

###SV的名字分割出sv对应的物种taxid
merge_sv_col = pd.read_csv('CD/FDR/5.1 SV//train_merge_sv_index.txt',sep='\t',index_col=0)
merge_sv_col[['Taxid', 'Other']] = merge_sv_col['sv_name'].str.split('.',expand=True)
merge_sv_col['Taxid'] = merge_sv_col['Taxid'].astype(int)
merge_sv_col

###提取出差异的sv
diff_biomarker = pd.read_csv('CD/FDR/5.1 SV/4_differential_signature.txt',sep='\t',index_col=0)
diff_feature = merge_sv_col.iloc[diff_biomarker.columns]
diff_feature


###计算sv的平均长度
def extract_info(string):
    # 使用正则表达式提取信息
    match = re.match(r'(\d+)\.(PRJ\w+\d+):(\d+_\d+(?:;\d+_\d+)*)', string)
    
    if match:
        # 提取信息
        group1 = match.group(1)
        group2 = match.group(3)
        
        # 处理多组数字的情况
        if ';' in group2:
            # 计算数字差
            numbers = [int(x.split('_')[1]) - int(x.split('_')[0]) for x in group2.split(';')]
            count = len(numbers)
            # 计算平均值
            average_difference = sum(numbers) / count
        else:
            # 只有一组数字时，直接赋值
            count = 1
            average_difference = int(group2.split('_')[1]) - int(group2.split('_')[0])
        
        return group1, group2, count,average_difference
    else:
        # 无法匹配的情况，输出到控制台
        print(f"Cannot match string: {string}")
        return None

data_list = diff_feature['sv_name']
result = [extract_info(item) for item in data_list if extract_info(item) is not None]
diff_details = pd.DataFrame(result, columns=["Taxid", "Regions", "Region_count","mean_length"])
# 增加一列为原始的data_list
diff_details["sv_name"] = list(data_list)
diff_details['sv_type'] = list(diff_feature['sv_type'])
diff_details['Taxid'] = pd.to_numeric(diff_details['Taxid'], errors='coerce', downcast='integer')
diff_details

###每一个差异SV的物种情况
diff_genomes = pd.merge(diff_details,species_reference,on='Taxid')
print(pd.value_counts(diff_genomes['sv_type']))
print(len(pd.value_counts(diff_genomes['Species'])))
diff_genomes

####计算差异SV归属的各个物种的具体情况
counts_df = diff_genomes.groupby(['Species', 'organism','sv_type']).size().unstack(fill_value=0).reset_index()
avg_length_df = diff_genomes.groupby(['Species','sv_type'])['mean_length'].mean().unstack(fill_value=0).reset_index()
result_df = pd.merge(counts_df, avg_length_df, on='Species', suffixes=('_No', '_Average'))
result_df.columns = ['Species','organism', 'No.dSV', 'No.vSv', 'Average_length.dSV', 'Average_length.vSV']
result_df

###将每个sv按照其region进行分割
def process_sv_data(merge_sv_col):
    # 初始化空的DataFrame
    df = pd.DataFrame(columns=['sv_name', 'sv_type', 'Taxid', 'project_id', 'Species','start', 'end'])

    # 处理每个字符串
    for i in range(len(merge_sv_col['sv_name'])):
        s = merge_sv_col['sv_name'][i]
        sv_type = merge_sv_col['sv_type'][i]
        Species = merge_sv_col['Species'][i]
        parts = s.split(':')

        # 提取taxid和project_id
        taxid_project_id = parts[0].split('.')
        taxid = taxid_project_id[0]
        project_id = taxid_project_id[1]

        # 提取以冒号分隔的第二部分，并按分号分割成组
        groups = parts[1].split(';')

        # 处理每个组并添加到DataFrame
        data0 = []
        for group in groups:
            start_end = group.split('_')
            start = start_end[0]
            end = start_end[1]

            # 将提取的数据添加到临时列表
            data0.append({'sv_name': s, 'sv_type': sv_type, 'Taxid': taxid, 'project_id': project_id,'Species': Species,'start': start, 'end': end})

        # 使用concat将临时列表转换为DataFrame并附加到主DataFrame
        df = pd.concat([df, pd.DataFrame(data0)])

    # 重置索引
    df.reset_index(drop=True, inplace=True)

    return df

segments_details = process_sv_data(diff_genomes)
print(segments_details)
segments_details['start'] = segments_details['start'].astype(float)
segments_details['end'] = segments_details['end'].astype(float)

##将sv的region情况与SGVFinder的物种的基因组情况结合
all_segments_details = pd.merge(segments_details,gene_reference,on=['Taxid', 'project_id'],how='left')

###仅保留影响基因区域的regions
filtered_rows = all_segments_details[
    (all_segments_details['start'] <= all_segments_details['end_pos']) &
    (all_segments_details['end'] >= all_segments_details['start_pos'])
]
filtered_rows

###使用patric中的annotation
SV_pathway_reference = pd.read_table('CD/FDR/5.1 SV/SV_patric_reference/merged_file.csv',sep=',')
diff_pathway = pd.merge(filtered_rows,SV_pathway_reference,left_on='gene_name',right_on='Refseq Locus Tag',how='left')

KO_map = pd.read_csv('CD/pathway/KO_map.csv')
diff_pathway_map = pd.merge(diff_pathway[~diff_pathway['Pathway Name'].isna()],KO_map,left_on='Pathway Name',right_on='Pathway',how='left')
print(len(diff_pathway_map[['sv_name', 'start', 'end','sv_type']].drop_duplicates()))
print(len(diff_pathway_map[['sv_name','sv_type']].drop_duplicates()))


#InDel

###找出chrom和assembly的对应关系
directory = 'CD/indel/reference/genbank/'
all_reference_details = pd.DataFrame(columns=['Assembly_ID', 'CHROM'])

for filename in os.listdir(directory):
    if filename.endswith(".fna"):
        
        file_path = os.path.join(directory, filename)
        
        reference_details = pd.DataFrame(columns=['Assembly_ID', 'CHROM'])
        
        assembly_id = 'GCA_' + filename.split('_')[1]  # 根据文件名生成 Assembly_ID

        sequence_ids = []
        
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    sequence_id = line.split()[0][1:]
                    sequence_ids.append(sequence_id)
        for seq_id in sequence_ids:
            reference_details = pd.concat([reference_details,pd.DataFrame([[assembly_id,seq_id]], columns=reference_details.columns)])
        print(assembly_id)
        print(len(sequence_ids))
        print(len(reference_details))
    all_reference_details = pd.concat([all_reference_details, reference_details])
    print(len(all_reference_details))
print(all_reference_details)

###找出assembly的物种详情
with open('CD/indel/reference/Assembly_ID.txt', 'r') as file:
     lines = [line.strip() for line in file]
Assembly_summary = pd.read_csv('CD/indel/reference/Assembly_summary.txt',sep='\t')
Assembly_summary['Assembly_ID'] = lines
Assembly_summary

all_info = pd.merge(all_reference_details,Assembly_summary,on='Assembly_ID',how = 'outer')
all_info

###所有的符合标准Indels，得出所有indels的详细情况（归属物种）
details = pd.read_csv('CD/FDR/5.2 InDel//detail.txt',sep='\t',index_col=0)
df = pd.merge(details,all_info,on='CHROM',how='left')
df

###计算indel的长度和判断类型
df['Length_DIFF'] = df.apply(lambda row: len(row['REF']) - len(row['ALT']), axis=1)
df['Type'] = df['Length_DIFF'].apply(lambda diff: 'deletion' if diff > 0 else 'insertion' if diff < 0 else 'no_change')
df['Length_DIFF'] = abs(df['Length_DIFF'])
df

##差异InDel的feature
diff_feature = pd.read_csv('CD/FDR/5.2 InDel//4_differential_signature.txt',sep='\t',index_col=0)
diff_feature = df.iloc[diff_feature.columns.astype(int),:]
diff_feature

def process_chrom(chrom):
    # 如果字符串中包含下划线，则以下划线分割并返回第二部分
    if '_' in chrom:
        return chrom.split('_')[1]
    # 如果没有下划线，则直接返回原始字符串
    else:
        return chrom

# 将函数应用于'CHROM'列的每个元素
diff_feature['CHROM'] = diff_feature['CHROM'].str.split('.').str[0].apply(process_chrom)
diff_feature['Length_DIFF'] = diff_feature.apply(lambda row: len(row['REF']) - len(row['ALT']), axis=1)

diff_feature['Type'] = diff_feature['Length_DIFF'].apply(lambda diff: 'deletion' if diff > 0 else 'insertion' if diff < 0 else 'no_change')
diff_feature['Length_DIFF'] = abs(diff_feature['Length_DIFF'])
diff_feature

###结合patric的reference，需要下载每个物种对应的feature和pathway两个表作为reference
feature = pd.read_csv('CD/indel/indel_reference/merged_feature.csv',index_col=0)
feature = feature[['Accession','BRC ID','Start','End']]
pathway = pd.read_csv('CD/indel/indel_reference/merged_pathway.csv',index_col=0)
pathway = pathway[['BRC ID','Refseq Locus Tag','Gene','Product','Pathway Name','EC Description']]

merge_reference = pd.merge(feature,pathway,on='BRC ID').drop_duplicates()
merge_reference

###map差异indel和patric构建的reference
diff_pathway = pd.merge(diff_feature,merge_reference,left_on='CHROM',right_on='Accession',how='left').drop_duplicates()
diff_pathway['Indel_End'] = diff_pathway['POS'] + diff_pathway['Length_DIFF']
diff_pathway 

filtered_rows = diff_pathway[
    (diff_pathway['POS'] <= diff_pathway['End']) &
    (diff_pathway['Indel_End'] >= diff_pathway['Start'])
]
filtered_rows

KO_map = pd.read_csv('CD/pathway/KO_map.csv')
diff_pathway_map = pd.merge(filtered_rows[~filtered_rows['Pathway Name'].isna()],KO_map,left_on='Pathway Name',right_on='Pathway',how='left')
diff_pathway_map


#SNV (mainly for Bacteroides_vulgatus_57955)
snp = pd.read_csv('CD/SNV/train_snv/Bacteroides_vulgatus_57955/snps_info.txt',sep = '\t')
snp

rep = pd.read_csv('CD/SNV/train_snv/Bacteroides_vulgatus_57955/genome.features',sep = '\t')
rep

diff = pd.read_csv('CD/FDR/5.3 SNV//midas/4_differential_signature.txt',sep='\t',index_col=0)
diff = diff.columns
filtered_diff = [s for s in diff if s.startswith('Bacteroides_vulgatus_57955')]
print(len(filtered_diff))

data = [s.split('_ID') for s in diff]
df = pd.DataFrame(data, columns=['species', 'site_id'])
df['site_id'] = df['site_id'].astype(int)
df

bv_list = df[df['species']=='Bacteroides_vulgatus_57955']
bv_list['site_id'] = bv_list['site_id'].astype(int)
bv_list
bv_detail = pd.merge(bv_list,snp,how='left',on='site_id')
bv_detail

ref_pathway = pd.read_csv('CD/SNV/train_snv/Bacteroides_vulgatus_57955/BVBRC_pathways_genes.csv',sep=',')
ref_pathway['BRC ID'] = ref_pathway['BRC ID'].str.replace('^fig\|', '', regex=True)
print(pd.value_counts(ref_pathway['Pathway Name']))
ref_pathway

function = pd.merge(bv_detail,rep,on='gene_id',how='left')
diff_pathway = pd.merge(function,ref_pathway,how='left',left_on='gene_id',right_on='BRC ID')
diff_pathway

KO_map = pd.read_csv('CD/pathway/KO_map.csv')
print(KO_map)
diff_pathway_map = pd.merge(diff_pathway[~diff_pathway['Pathway Name'].isna()],KO_map,left_on='Pathway Name',right_on='Pathway',how='left')
print(diff_pathway_map)

####计算有义突变的数量
kska = diff_pathway[~diff_pathway['Pathway Name'].isna()]
kska[['amino_acid_1', 'amino_acid_2', 'amino_acid_3', 'amino_acid_4']] = kska['amino_acids'].str.split(',', expand=True)
# 创建一个新的索引列表
new_index = list(range(len(kska)))
kska = kska.rename(index=dict(zip(kska.index, new_index)))
kska
ks=0
for i in range(len(kska)):
    if ((kska.iloc[i,5] == 'A')|(kska.iloc[i,5] == 'C')) &((kska.iloc[i,6] == 'C')|(kska.iloc[i,6] == 'A')):
        if kska.loc[i,'amino_acid_1'] != kska.loc[i,'amino_acid_2']:
            ks +=1
    if ((kska.iloc[i,5] == 'A')|(kska.iloc[i,5] == 'G')) &((kska.iloc[i,6] == 'G')|(kska.iloc[i,6] == 'A')):
        if kska.loc[i,'amino_acid_1'] != kska.loc[i,'amino_acid_3']:
            ks +=1
    if ((kska.iloc[i,5] == 'A')|(kska.iloc[i,5] == 'T')) &((kska.iloc[i,6] == 'T')|(kska.iloc[i,6] == 'A')):
        if kska.loc[i,'amino_acid_1'] != kska.loc[i,'amino_acid_4']:
            ks +=1
    if ((kska.iloc[i,5] == 'C')|(kska.iloc[i,5] == 'G')) &((kska.iloc[i,6] == 'G')|(kska.iloc[i,6] == 'C')):
        if kska.loc[i,'amino_acid_2'] != kska.loc[i,'amino_acid_3']:
            ks +=1
    if ((kska.iloc[i,5] == 'C')|(kska.iloc[i,5] == 'T')) &((kska.iloc[i,6] == 'T')|(kska.iloc[i,6] == 'C')):
        if kska.loc[i,'amino_acid_2'] != kska.loc[i,'amino_acid_4']:
            ks +=1
    if ((kska.iloc[i,5] == 'G')|(kska.iloc[i,5] == 'T')) &((kska.iloc[i,6] == 'T')|(kska.iloc[i,6] == 'G')):
        if kska.loc[i,'amino_acid_3'] != kska.loc[i,'amino_acid_4']:
            ks +=1
print(ks)