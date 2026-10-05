from itertools import combinations
from numpy import log2, std, unique
from pandas import concat, isna, notna, read_excel, DataFrame
from pprint import pprint

# simple IDS-like JSON with easy-to-implement structure
ids = {'run_information': None,
       'results': {
           'quantification': {'amplification': None, 'cq': None},
           'gene_expression': {'original': None,
                               # 'pivoted': None,  # out of scope due to parameter variability
                               'm_values': None}
       }}

# Raw data

raw = read_excel(io='data/return/20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Quantification_Amplification_Results.xlsx',
                 sheet_name=None, engine='calamine')
pprint(raw)


raw_df = []
for sheet in raw:
    if not sheet == 'Run Information':
        tmp = raw[sheet].drop(columns='Unnamed: 0')
        tmp['Dye'] = sheet
        raw_df.append(tmp)
raw_df = concat(raw_df)
raw_df.info()
raw_df.head()

ids['results']['quantification']['amplification'] = raw_df.to_dict(orient='list')

# Cq

cq = read_excel(io='data/return/20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Quantification_Cq_Results.xlsx',
                sheet_name=0, engine='calamine')
cq.drop(columns='Unnamed: 0', inplace=True)
cq.info()
cq.head()
ids['results']['quantification']['cq'] = cq.to_dict(orient='list')


# Expression

# Use sheet "Run information" directly from Expression Excel file
run_info = read_excel(io='data/return/analyze_using_biological_group_sample_20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Gene_Expression_Results-Bar_Chart.xlsx',
                      sheet_name=1, engine='calamine', header=None)
mtdt = dict(zip(run_info[0].values.tolist(), run_info[1].values.tolist()))
ids['run_information'] = mtdt

expr = read_excel(io='data/return/analyze_using_biological_group_sample_20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Gene_Expression_Results-Bar_Chart.xlsx',
                  sheet_name=0, engine='calamine')
expr = read_excel(io='data/return/analyze_using_biological_groups_only_20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Gene_Expression_Results-Bar_Chart.xlsx',
                  sheet_name=0, engine='calamine')
expr = read_excel(io='data/return/analyze_using_sample_biological_group_20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Gene_Expression_Results-Bar_Chart.xlsx',
                  sheet_name=0, engine='calamine')
expr = read_excel(io='data/return/analyze_using_samples_only_20260310_SNB_7136_HCT116_fancm_KO_sirna_njd_MXTO-Gene_Expression_Results-Bar_Chart.xlsx',
                  sheet_name=0, engine='calamine')
expr.drop(columns='Unnamed: 0', inplace=True)
expr.info()
# expr.head()

analyze_using = expr.columns[2]

if analyze_using == 'Biological Group':
    raise RuntimeError('>> Output with setting `Analyze Using` = '
                       '`Biological Groups Only` is not supported')
else:
    expr['Analyze Using'] = expr.columns[2]
    expr.rename(columns={expr.columns[2]: 'Analyze Using Value'}, inplace=True)
    expr = expr.reindex(columns=expr.columns[0:2].tolist() + \
                        ['Analyze Using'] + \
                        expr.columns[2:-1].tolist())
expr.info()

# out of scope because of `Analyze using` variability
# pivoted = expr.pivot(columns=['Target'], index='Biological Group Sample',
#                      values=['Expression', 'Expression SEM'])
# pivoted.to_csv()

# view matching report export
# expr[["Target", "Biological Group Sample", "Control",
#       "Expression", "Expression SEM", "Corrected Expression SEM",
#       "Mean Cq", "Cq SEM", "P-Value"]]

# Calculating geNorm M-values

# If column expression is fully null, abort calculation
if int(isna(expr.Expression).sum()) == expr.shape[0]:
    raise RuntimeError('>> Column `Expression` is empty, have you selected '
                       'reference genes?')
else:
    data = expr[['Target', 'Expression', 'Mean Cq']]
    # Y targets have Expression = NA but Mean Cq != NA
    y_targets = data[data.Expression.isna() & data['Mean Cq'].notna()]\
        .Target.unique().tolist()

data = expr[['Analyze Using Value', 'Target', 'Mean Cq']]
data = data.pivot(columns=['Target'], index=['Analyze Using Value'],
                  values='Mean Cq')

m_value_threshold = 0.5

combs = list(combinations(iterable=y_targets, r=2))
pprint(combs)
res = []
for comb in combs:
    res.append({
        'num': comb[0],
        'denom': comb[1],
        'std_log2_cq_mean_ratio': std(log2(data[comb[0]] /
                                           data[comb[1]]))
    })
res = DataFrame.from_records(res)
print('res =', res, sep='\n')

m_values = []
for trgt in y_targets:
    m_value = res.query('num == @trgt | denom == @trgt')\
        ['std_log2_cq_mean_ratio'].mean()
    m_values.append({
        'Gene': trgt,
        'M_value': m_value,
        'M_value_threshold': m_value_threshold,
        'Selected': 'Yes' if m_value < m_value_threshold else 'No'
    })
m_values = DataFrame.from_records(m_values)
print('m_values =', m_values, sep='\n')

selected_ref_genes = m_values.query('Selected == "Yes"')['Gene'].tolist()
print(selected_ref_genes)

ids['results']['gene_expression']['original'] = expr.to_dict(orient='list')
# ids['results']['gene_expression']['pivoted'] = pivoted.to_dict(orient='list')  # TODO collapse colnames
ids['results']['gene_expression']['m_values'] = m_values.to_dict(orient='list')
pprint(ids)
