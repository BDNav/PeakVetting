'''
Written by Jonathan Monk for Naviaux Lab
Verion 1.27  4/7/2021
'''

from IPython import embed
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

GCPV_version = 1.27

# 1.27 adding error checks for dupe mets in AUC and that mets in RT match those in AUC



print('GCMS Peak Vetter version %0.2f'%GCPV_version)
import sys
version = sys.version_info
print('You are running python %i.%i'%(version.major, version.minor))
print('pandas version %s'%pd.__version__)

# choose a name for the study - the output file will reflect this:
study_name = "XXX"

# input AUC/SN filename (in excel, .xlsx):
filename = 'XXX.xlsx' 

# input RT filename (in excel, .xlsx):
RT_filename = 'XXX.xlsx'

# the name of the blank samples in the "Group" Column:
blank_name = 'Blank'

# choose case and control names for stats tab
CASE='case'
CONTROL='control'

#the name of the QC samples in the "Group" Column
qc_name = 'QC'

print('AUC File: %s'%(filename))
print('RT File: %s'%(RT_filename))
print('Blank name: %s'%(blank_name))
print('Case name: %s'%(CASE))
print('Control name: %s'%(CONTROL))
print('QC name: %s'%(qc_name))



mrm1_SN_cutoff = 0 # filter 1
mrm2_ratio_cutoff = 3 # filter 2
mrm1_ratio_cutoff = 10 # filter 3
mrm1_area_cutoff = 2000 # filter 4
RT_CUTOFF = 0.2 # filter 5
mrm2_SN_cutoff = 3 # filter 6

GROUP_COL_NAME=('Group','Group')
GROUP_COL = 1 # set the location of the group column


''' LOAD STUDY RETENTION TIMES '''
retention_times = pd.read_excel(RT_filename)
retention_times = retention_times.iloc[1:,:]
compounds = [c.replace(' Results','') for c in retention_times.columns]
to_rename = dict(zip(retention_times.columns.tolist(),compounds))
retention_times = retention_times.rename(columns=to_rename)

# ensure compounds = retention times:
if len(set(retention_times.columns.tolist())-set(compounds))!=0:
    print('Are you sure your retention times columns match those in your MRM file? This may cause an error')


for c in retention_times[retention_times.columns[2:]]:
    retention_times[c] = pd.to_numeric(retention_times[c])
retention_times.index=retention_times['Sample']
retention_time_names = retention_times.columns.tolist()[2:]


data = pd.read_excel(filename, header=[0,1], dtype={('Group','Group'):str})
data.index=data[('Sample','Name')]

blanks = data[data[GROUP_COL_NAME]==blank_name].index
qc_names = data[data[GROUP_COL_NAME]==qc_name].index
values = list(set(data.index.tolist()) - set(blanks.tolist()) - set(qc_names.tolist()))

if len(values)+len(qc_names)+len(blanks) != len(data):
    print('error - num of blanks + qc + values != total data size')


# REPLACE INFINITY SYMBOL WITH 20
inf = b'\xe2\x88\x9e'.decode('utf-8')
data = data.replace(inf, 20)
# REPLACE 0's with NA
data = data.replace(0, np.nan)

# fill blank areas with 1/10 of min measured val
print('filling blanks with 0 or NA to 1/10 of lowest measured val')
for c in data.columns:
    if c[1]=='Area':
        min = data.loc[blanks,c].min()
        replace_val = min*0.1
        new_vals = data.loc[blanks,c].fillna(replace_val).replace(0,replace_val)
        data.loc[blanks,c]=new_vals


outfile_name = study_name + "_vetted_peaks_py%i%i.xlsx"%(version.major, version.minor)
outfile = outfile_name

outfile = pd.ExcelWriter(outfile)

# 1 re-org sheet

mrm1_areas = {} # was areas
mrm2_areas = {} # was mrm3
mrm1_sns = {} # was sns
mrm2_sns = {}

groups = data[data.columns[GROUP_COL]].unique().tolist()
# groups = data[GROUP_COL_NAME].unique().tolist()
print('\nRunning...')
print('Found %i unique groups:'%len(groups))

print(groups)
groups.remove(blank_name)

# RESTRUCTURE THE DATAFRAME

unique_compounds = []

# NOW only want MRM1 for area and MRM2 for SN

for i in range(GROUP_COL+1,len(data.columns),4):
    name = data.columns[i][0].replace(' Results','')
    # print name
    unique_compounds.append(name)
    
    # script is dependent on structure: every 4 columns is a new metabolite
    
    area1 =  data.columns[i]
    area1 = data[area1]
    
    sn1 =  data.columns[i+1]
    sn1 = data[sn1]
    
    area2 = data.columns[i+2]
    area2 = data[area2]
    
    sn2 = data.columns[i+3]
    sn2 = data[sn2]
    
    mrm1_areas[name] = area1
    mrm2_areas[name] = area2
    mrm1_sns[name] = sn1
    mrm2_sns[name] = sn2
        
mrm1_areas = pd.DataFrame(mrm1_areas)
mrm2_areas = pd.DataFrame(mrm2_areas)
mrm2_sns = pd.DataFrame(mrm2_sns)
mrm1_sns = pd.DataFrame(mrm1_sns)

mrm1_areas.to_excel(outfile,'MRM1_AREA')
mrm1_sns.to_excel(outfile,'MRM1_SN')

mrm2_areas.to_excel(outfile,'MRM2_AREA')
mrm2_sns.to_excel(outfile,'MRM2_SN')


# CHECK for dupe compounds
tmp = pd.Series(unique_compounds)
cnts = tmp.value_counts()
dupes = cnts[cnts>1]
if len(dupes) > 0:
    print('ERROR you have duplicate columns in the AUC file for the following compounds:')
    print(dupes)
    embed()

# FILL IN BLANKS:
from scipy import stats
import numpy as np
import random

def fill_blanks(blank_areas):

    gmeans = {}
    for c in blank_areas.columns:
        measured_values = blank_areas[c].dropna()
        measured_values = measured_values[measured_values>0]
        if len(measured_values)==0:
            gmean = 200 # set analytes with no blanks to 200
        else:   
            gmean = stats.gmean(measured_values)
            std = blank_areas[c].std() # TO DO sample these - but keep all same for now... sample code below:
            
            mean = measured_values.mean()
            sigma = std/mean
            sample = np.random.lognormal(1, .2, 1000)*gmean
            gmeans[c]=gmean
        
        blank_areas[c] = blank_areas[c].apply(lambda v: gmean if (str(v)=='nan') else v) # for testing - so results are consistent
    return blank_areas
    
blank_areas_mrm1 = mrm1_areas.loc[blanks,:]
res = fill_blanks(blank_areas_mrm1)
mrm1_areas.loc[res.index,:] = res

blank_areas_mrm2 = mrm2_areas.loc[blanks,:]
res2 = fill_blanks(blank_areas_mrm2)
mrm2_areas.loc[res2.index,:] = res2

out = {}
mrm1_area_counts = mrm1_areas.transpose().count()
out[('Global Sums','MRM1_AREA_present')] = mrm1_area_counts

missing_mrm1_area_counts = len(mrm1_areas.transpose()) - mrm1_areas.transpose().count()
out[('Global Sums','MRM1_AREA_missing')] = missing_mrm1_area_counts

mrm2_area_counts = mrm2_areas.transpose().count()
out[('Global Sums','MRM2_AREA_present')] = mrm2_area_counts

missing_mrm2_area_counts = len(mrm2_areas.transpose()) - mrm2_areas.transpose().count()
out[('Global Sums','MRM2_AREA_missing')] = missing_mrm2_area_counts

mrm1_sns_counts = mrm1_sns.transpose().count()
out[('Global Sums','MRM1_SN_present')] = mrm1_sns_counts

missing_mrm1_sns_counts = len(mrm1_sns.transpose()) - mrm1_sns.transpose().count()
out[('Global Sums','MRM1_SN_missing')] = missing_mrm1_sns_counts

missing_mrm2_sns_counts = len(mrm2_sns.transpose()) - mrm2_sns.transpose().count()
out[('Global Sums','MRM2_SN_missing')] = missing_mrm2_sns_counts

rt_counts = retention_times.transpose().count()
out[('Global Sums','RT_present')] = rt_counts

missing_rt_counts = len(retention_times.transpose()) - retention_times.transpose().count()
out[('Global Sums','RT_missing')] = missing_rt_counts

mrm1_areas.to_excel(outfile, 'MRM1_AREAS_blanks_filled')
mrm2_areas.to_excel(outfile, 'MRM2_AREAS_blanks_filled')
# embed()

blank_means = mrm1_areas.loc[blanks,:].mean()
ratios = mrm1_areas/blank_means

ratios.to_excel(outfile, 'MRM1_AREA_RATIOS')

blank_means2 = mrm2_areas.loc[blanks,:].mean()
ratios2 = mrm2_areas/blank_means2

ratios2.to_excel(outfile, 'MRM2_AREA_RATIOS')



# 1. Filter 1 = MRM1 SN > SN cutoff
# 2. Filter 2 = AUC/Blank >3
# 3. Filter 3 = MRM1 AUC/Blank >10
# 4. Filter 4 = MRM1 AUC >2000
# 5. Filter 5 = RT > 0.2
# 6. Filter 6 = MRM2 > SN cutoff

# Filter 1 = MRM2 is not null
filter1_name = 'Filter 1: MRM1 S/N > %i'%mrm1_SN_cutoff
filter1 = mrm1_sns.where(mrm1_sns>mrm1_SN_cutoff)# 
filter1_exclude = mrm1_sns.where(mrm1_sns<=mrm1_SN_cutoff)

out[(filter1_name,'passed')] = filter1.transpose().count()
out[(filter1_name,'failed')] = mrm1_sns_counts - filter1.transpose().count()

# Filter 2 = MRM2 AUC/Blank>3
filter2_name = 'Filter 2 = MRM2 AUC/Blank >= %i'%mrm2_ratio_cutoff
filter2 = filter1.where(ratios2>=mrm2_ratio_cutoff)# 
filter2_out = ratios2.where(~filter2.isna())
filter2_exclude = filter1.where(ratios2<mrm2_ratio_cutoff)# 

out[(filter2_name,'passed')] = filter2.transpose().count()
out[(filter2_name,'failed')] = filter2_exclude.transpose().count()


# Filter 3 = MRM1 AUC/Blank > 10s
filter3_name = 'Filter 3 = MRM1 AUC/Blank >= %i'%mrm1_ratio_cutoff
filter3 = filter2.where(ratios>=mrm1_ratio_cutoff)# 
filter3_out = ratios.where(~filter3.isna())
filter3_exclude = filter2.where(ratios<mrm1_ratio_cutoff)# 

out[(filter3_name,'passed')] = filter3.transpose().count()
out[(filter3_name,'failed')] = filter3_exclude.transpose().count()


# Filter 4 = MRM1 > 2000
filter4_name = 'Filter 4 = MRM1 AUC >= %i'%mrm1_area_cutoff
filter4 = filter3.where(mrm1_areas>=mrm1_area_cutoff)# 
filter4_out = mrm1_areas.where(~filter4.isna())
filter4_exclude = filter3.where(mrm1_areas<mrm1_area_cutoff)# 

out[(filter4_name,'passed')] = filter4.transpose().count()
out[(filter4_name,'failed')] = filter4_exclude.transpose().count()


# Filter 5 = RT < 0.2
filter5_name = 'Filter 5 = RT < %0.2f'%RT_CUTOFF
filter5 = filter4.where(retention_times[unique_compounds]<RT_CUTOFF)# 
filter5_out = retention_times.where(~filter5.isna())
filter5_exclude = filter4.where(retention_times[unique_compounds]>=RT_CUTOFF)# 

out[(filter5_name,'passed')] = filter5.transpose().count()
out[(filter5_name,'failed')] = filter5_exclude.transpose().count()


# filter 6:
filter6_name = 'Filter 6 = MRM 2 S/N >= %i'%mrm2_SN_cutoff
filter6 = filter5.where(mrm2_sns>=mrm2_SN_cutoff)# 
filter6_out = mrm2_sns.where(~filter6.isna())
filter6_exclude = filter5.where(mrm2_sns<mrm2_SN_cutoff)# 

out[(filter6_name,'passed')] = filter6.transpose().count()
out[(filter6_name,'failed')] = filter6_exclude.transpose().count()



out = pd.DataFrame(out)

out = out[['Global Sums',filter1_name, filter2_name, filter3_name, filter4_name, filter5_name]]

filter1[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter1 = filter1[[GROUP_COL_NAME]+filter1.columns[:-1].tolist()]
filter2_out[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter2_out = filter2_out[[GROUP_COL_NAME]+filter2_out.columns[:-1].tolist()]
filter3_out[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter3_out = filter3_out[[GROUP_COL_NAME]+filter3_out.columns[:-1].tolist()]
filter4_out[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter4_out = filter4_out[[GROUP_COL_NAME]+filter4_out.columns[:-1].tolist()]

del filter5_out['Sample']
del filter5_out['Group']

filter5_out[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter5_out = filter5_out[filter4_out.columns.tolist()]

filter6_out[GROUP_COL_NAME] = data[GROUP_COL_NAME]
filter6_out = filter6_out[[GROUP_COL_NAME]+filter6_out.columns[:-1].tolist()]



out.to_excel(outfile, 'Troubleshooting')
filter1.to_excel(outfile, 'Filter 1 passed')
filter2_out.to_excel(outfile, 'Filter 2 passed')
filter3_out.to_excel(outfile, 'Filter 3 passed')
filter4_out.to_excel(outfile, 'Filter 4 passed')
filter5_out.to_excel(outfile, 'Filter 5 passed')
filter6_out.to_excel(outfile, 'Filter 6 passed')

''' adding new sheet for passed with AUCS'''
passed_AUCS = mrm1_areas.where(~filter6.isna())
passed_AUCS[GROUP_COL_NAME] = data[GROUP_COL_NAME]
passed_AUCS = passed_AUCS[[GROUP_COL_NAME]+passed_AUCS.columns[:-1].tolist()]

pos_counts = passed_AUCS.count()
to_keep = pos_counts[pos_counts>0]
passed_AUCS_only_pos = passed_AUCS[to_keep.index]
passed_AUCS_only_pos = passed_AUCS_only_pos.sort_values(GROUP_COL_NAME)
passed_AUCS_only_pos.to_excel(outfile, 'Positive MRM1 AUCs')

''' FALSE POSITVES TAB'''
false_positives = passed_AUCS.loc[blanks,:]
to_keep = false_positives.count()
to_keep = to_keep[to_keep>0]
false_positives = false_positives[to_keep.index]
false_positives.to_excel(outfile, 'False Positives')


''' QC TAB '''
qc_tab = passed_AUCS.loc[qc_names,:]
to_keep = qc_tab.count()
to_keep = to_keep[to_keep>0]
qc_tab = qc_tab[to_keep.index]
qc_tab.to_excel(outfile, 'QC')

count_positives = passed_AUCS

count_positives[GROUP_COL_NAME] = data[GROUP_COL_NAME]

all_positives = {}

count_positives = count_positives[count_positives[GROUP_COL_NAME]!=blank_name]

for c in unique_compounds:
    mrm1 = c
   
    if mrm1 in count_positives.columns:
        tmp = count_positives[[GROUP_COL_NAME,mrm1]].dropna()
        res = tmp.groupby(GROUP_COL_NAME).count()
        res_sum = res.sum()
        # print res.to_dict()
        # embed()
        all_positives[c] = res.to_dict()[mrm1]
        all_positives[c]['sum'] = res_sum[0]
all_positives = pd.DataFrame(all_positives).transpose().fillna(0)

all_positives = all_positives[all_positives['sum']>0]

vetted = mrm1_areas

# color areas according to whether the ratio is > ratio cutoff
def apply_color(x):
    def myfunc(x):
        from IPython import embed
        # embed()
        if type(x)==str or type(x)==unicode:
            return 'background-color: #ffffff'
        elif str(x) !='nan':
            return 'background-color: #90EE90'
        else:
            return 'background-color: #f5f5f5'
        
    
    return sn_masked_ratios_filtered.applymap(lambda val: myfunc(val))
    

all_positives.to_excel(outfile, 'count positives')


import datetime as dt

readme = {
    'Date of Analysis':	dt.datetime.today().strftime("%m/%d/%Y"),
    'Columns':'GCMS',
    'Groups':','.join(groups),
    'Study Name':study_name,
    'Input File Name':filename,
    'Filter 1: MRM1 S/N cutoff': mrm1_SN_cutoff,
    'Filter 2: MRM2 Area Ratio cutoff': mrm2_ratio_cutoff,
    'Filter 3: MRM1 Area Ratio cutoff': mrm1_ratio_cutoff,
    'Filter 4: MRM1 Area cutoff': mrm1_area_cutoff,
    'Filter 5: RT cutoff': RT_CUTOFF,
    'Filter 6: MRM2 S/N cutoff': mrm2_SN_cutoff,
    'Chemicals Targeted	':len(mrm2_sns.columns),
    'Chemicals Vetted ':len(all_positives),
    'Note': 'Any analytes with no blank values were set to 200',
    'GCMS PeakVet Version':GCPV_version
}
readme = pd.Series(readme)
readme.to_excel(outfile, 'readme')    

''' ADDING STATS info '''
data = passed_AUCS
targeted_chemicals = len(data.columns)-1

dt = data[data.columns[1:]].transpose() # transposed data
data['TML'] = dt.sum() 

data.loc[:,'hits'] = data[data.columns[1:]].transpose().count() # first col is group...


case_data = data[data[GROUP_COL_NAME]==CASE]
control_data = data[data[GROUP_COL_NAME]==CONTROL]

data = data[data[GROUP_COL_NAME].isin([CASE,CONTROL])]


# import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu, ttest_ind

from scipy.stats.mstats import gmean

# organize hit vs no hit data for fishers exact:
binary = data[data.columns[1:-1]].fillna(0)
binary[binary>0]=1
binary = binary.replace(0,'Not')
binary = binary.replace(1,'Hit')
binary[GROUP_COL_NAME] = data[GROUP_COL_NAME]

N = len(data)
CASE_N = len(case_data)
CONTROL_N = len(control_data)




out = []

# tml for all aucs
case_tml = case_data[case_data.columns[1:-1]].sum()
control_tml = control_data[control_data.columns[1:-1]].sum()

case_tml = case_tml[case_tml>0]
control_tml = control_tml[control_tml>0]

gmean_case_tml = gmean(case_tml)
gmean_control_tml = gmean(control_tml)

# case vs control tml
tml_stat, tml_mannw_pval = stats.mannwhitneyu(case_tml,control_tml, use_continuity=False, alternative='two-sided')

for c in data.columns[1:-1]:
    
    fisher_pval = np.nan
    ttest_pval = np.nan
    welch_pval = np.nan
    mannw_pval = np.nan
    
    case_counts = float(case_data[c].count())
    control_counts = float(control_data[c].count())
    
    percent_case = case_counts/CASE_N
    percent_control = control_counts/CONTROL_N
    
    if control_counts==0:
        count_ratio = np.nan
    else:    
        count_ratio = percent_case/percent_control # this is how bob calculated this

    case = case_data[c].dropna()
    control = control_data[c].dropna()
    
    case_gmean = gmean(case)
    control_gmean = gmean(control)
    
    gmean_ratio = case_gmean/control_gmean
    
    c_case_tml = case.sum() 
    c_control_tml = control.sum()
    
    
    tab = pd.crosstab(binary[GROUP_COL_NAME], binary[c]).transpose()
    if len(tab)>1:
        oddsratio, fisher_pval = fisher_exact(tab)
    rej, ttest_pval = ttest_ind(case, control) # students ttest
    rej, welch_pval = ttest_ind(case, control, equal_var=False) # welchs ttest
    try:
        stat, mannw_pval = stats.mannwhitneyu(case,control, use_continuity=False, alternative='two-sided') # NOTE: these parameters will make the MWpval=KWpval (for two samples)  UPDATED 10/26/2020 per Bob email
    except:
        print('MannWhitney error: case == control')
    
    if mannw_pval==0:
        mannw_pval = np.nan
    
    out.append({
        'MRM Name':c,
        '%s Counts'%CASE:case_counts,
        '%s Counts'%CONTROL:control_counts,
        'Percent %s Counts/%i'%(CASE, CASE_N):percent_case,
        'Percent %s Counts/%i'%(CONTROL, CONTROL_N):percent_control,
        'Count Ratio %s/%s'%(CASE, CONTROL):count_ratio,
        'Fisher\'s Exact test p value':fisher_pval,
        'Student\'s AUC ttest':ttest_pval,
        'Welch\'s AUC ttest p value':welch_pval,
        'Mann-Whitney AUC U-test value':mannw_pval,
        '%s Geomean'%CASE:case_gmean,
        '%s Geomean'%CONTROL:control_gmean,
        'Geomean Ratio %s/%s'%(CASE,CONTROL):gmean_ratio,
    })
      
    
out = pd.DataFrame(out)
out = out.sort_values('Fisher\'s Exact test p value')
out.index = out['MRM Name']

out_data = data.transpose()
out_data = pd.merge(out_data, out, left_index=True, right_index=True, how='left')

del out_data['MRM Name']

print('saving %s'%outfile_name)
out_data.to_excel(outfile,'stats')
outfile.save()