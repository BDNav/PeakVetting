'''
Written by Jonathan Monk for Naviaux Lab
Verion 7.1  4/7/2021
'''


import pandas as pd
from IPython import embed

version = '7.1';
print('LCMS peak vetter version %s'%version)

# 1. Edit the name of the input file below:
filename = 'XXX.csv'

# 2. Give your study a name (this is reflected in the name of the output file)
study_name = 'XXX

# 3. set the thresholds:
MRM1_vs_BLANK_THRESHOLD = 10 # MRM1 > this threshold
# SN_THRESHOLD

# 4. for statistics need to choose case and controls:
CASE='Mp'
CONTROL='Fp'
GROUP_COL_NAME='Group'


# 5. Run this program by typing in >python vet_peaks_LCMS.py

data = pd.read_csv(filename)

data['Group'] = data['Group'].replace('blank','Blank')

groups = data['Group'].unique().tolist()

exp_groups = list(set(groups)-set(['Blank']))

group_indices = {}
for group in exp_groups:
    idx = data[data['Group']==group].index
    group_indices[group] = idx
    
compounds = list(set(['_'.join(c.split('_')[:-1]) for c in data.columns]))
blanks = data[data['Group']=='Blank'].index
samples = data[data['Group']!='Blank'].index # assuming "sample" means anything besides a blank

mrm1s = [c for c in data.columns if '_MRM1' in c]
mrm2s = [c for c in data.columns if '_MRM2' in c]

mrm1_rename = [c.replace('_MRM1','') for c in mrm1s]
mrm2_rename = [c.replace('_MRM2','') for c in mrm2s]

mrm1 = data[mrm1s]
mrm2 = data[mrm2s]

mrm1 = mrm1.rename(columns=dict(zip(mrm1s,mrm1_rename)))
mrm2 = mrm2.rename(columns=dict(zip(mrm2s,mrm2_rename)))

blank_mrm1 = mrm1.loc[blanks,:]
blank_mrm2 = mrm2.loc[blanks,:]

from scipy import stats
import numpy as np
import random

def fill_blank_gmean(data):
    for c in data.columns:
        if len(data[c].dropna())==0: # if no blank values set all to 1000
                data.loc[:,c]=1000
        if len(data[c])!=len(data[c].dropna()): # otherwise add randomly sample gmean values:
            gmeans = stats.gmean(data[c].dropna())
            sigma1=0.2 
            sample1 = np.random.lognormal(0, sigma1, 10000)*gmeans
            for i in data.index:
                if str(data.loc[i,c])=='nan':
                    rnd = random.choice(sample1)
                    data.loc[i,c] = int(rnd)
    return data
blank_mrm1_filled = fill_blank_gmean(blank_mrm1)
blank_mrm2_filled = fill_blank_gmean(blank_mrm2)

mrm1_samples = mrm1.loc[samples,:]
mrm2_samples = mrm2.loc[samples,:]

mrm1_ratio = mrm1_samples/blank_mrm1_filled.mean()
mrm2_ratio = mrm2_samples/blank_mrm2_filled.mean()


out = pd.DataFrame()
out2 = pd.DataFrame()
ratios = pd.DataFrame()
ratios2 = pd.DataFrame()
blank_avgs = pd.DataFrame()
case_control_positive = {}

singles = []
for c in compounds[1:]:
    single=False
    try: # check to see if there are two MRMs
        data2 = data[[c+'_MRM1',c+'_MRM2']]
    except:
        data2 = data[[c+'_MRM1']]
        singles.append(c)
        single=True
    
    avg = data2.iloc[blanks].mean()
    orig = data2.iloc[blanks]
    new = data2.iloc[blanks].dropna()
    
    # additions of geo mean and distr for un measured blanks
    
    if len(orig)!=len(new) and len(new)!=0:
        # print orig
        gmeans = stats.gmean(orig.dropna())
        if len(gmeans)==2:
            
            stds = orig.std().tolist()
            sigma1=0.2
            sigma2=0.2
           
            sample1 = np.random.lognormal(0, sigma1, 10000)*gmeans[0]
            sample2 = np.random.lognormal(0, sigma2, 10000)*gmeans[1]
            for cix in range(0,2):
                c = orig.columns[cix]
                for i in orig.index:
                    if str(orig.loc[i,c])=='nan':
                        if cix==0:
                            rnd = random.choice(sample1)
                        else:
                            rnd = random.choice(sample2)
                        
                        orig.loc[i,c] = int(rnd)
                        data.loc[i,c] = int(rnd)
                    
        elif len(gmeans)==1:
            print(c + ' has only one MRM')
            sigma1=0.2
           
            sample1 = np.random.lognormal(0, sigma1, 10000)*gmeans[0]
            c = orig.columns[0]
            for i in orig.index:
                if str(orig.loc[i,c])=='nan':
                    rnd = random.choice(sample1)
                    orig.loc[i,c] = int(rnd)
                    data.loc[i,c] = int(rnd)
          
    if len(avg[avg.isnull()])>0:
        
        avg = avg.fillna(1000)
    data3=data2.iloc[samples]
    res = data3/avg
    
    if res[res>MRM1_vs_BLANK_THRESHOLD].sum().sum()>0:
        masked = data3
    else:
        masked = res[res>MRM1_vs_BLANK_THRESHOLD]
    
    if len(res[res>MRM1_vs_BLANK_THRESHOLD].dropna())>0: # now looks for BOTH > 3
        masked2 = data3
        res2 = res
        
        # check for how many case/control pass the threshold
        tmp=res[res>MRM1_vs_BLANK_THRESHOLD].dropna()
        case_control_positive[c] = {}
        for group in exp_groups:
            idx = group_indices[group]
            count = len(set(tmp.index.tolist()) & set(idx.tolist()))
            case_control_positive[c][group] = count
        
    else:
        masked2 = res[res>MRM1_vs_BLANK_THRESHOLD].dropna()
        res2 = res[res>MRM1_vs_BLANK_THRESHOLD].dropna()
    
    out = pd.concat([out,masked], axis=1)
    out2 = pd.concat([out2,masked2], axis=1)
    ratios = pd.concat([ratios,res], axis=1)
    ratios2 = pd.concat([ratios2,res2], axis=1)
    blank_avgs = pd.concat([blank_avgs,avg])

case_control_positive = pd.DataFrame(case_control_positive).transpose()
case_control_positive['sum'] = case_control_positive.transpose().sum()
    

def apply_color(x):
    # lambda x: np.nan if x < 90 else x
    def myfunc(x):
        if type(x)==str:
            return 'background-color: #ffffff'
        elif x < MRM1_vs_BLANK_THRESHOLD:
            return 'background-color: #f5f5f5'
        else:
            return 'background-color: #90EE90'
        
    
    return ratios2.fillna(0).applymap(lambda val: myfunc(val))
    
# add labels back to vetted peaks:

labels = data.loc[out2.index, data.columns[:2]]

out2 = pd.concat([labels, out2], axis=1)
ratios2 = pd.concat([labels, ratios2], axis=1)

positives = out2.where(ratios>MRM1_vs_BLANK_THRESHOLD)
pos_counts = positives.count()
to_keep = pos_counts[pos_counts>0]
passed_AUCs = positives[to_keep.index.tolist()]

passed_cols = [c[:-5] for c in passed_AUCs.columns] # get unique cols (filter MRM1 and MRM2)
passed_cols = list(set(passed_cols))  # need to drop compounds that don't pass BOTH MRM1 and MRM2 (dropnas from both)
final_passed_AUCs = pd.DataFrame()
for c in passed_cols:
    tmp = passed_AUCs[list(set([c+'_MRM1',c+'_MRM2'])&set(passed_AUCs.columns.tolist()))].dropna()
    # print tmp
    final_passed_AUCs = pd.concat([final_passed_AUCs, tmp], axis=1)

passed_AUCs = final_passed_AUCs


passed_AUCs[['Sample Name', 'Group']] = out2[['Sample Name', 'Group']]

new_cols = passed_AUCs.columns[:-2].tolist()
new_cols.sort()
new_cols = [c for c in new_cols if '_MRM1' in c] 
passed_AUCs = passed_AUCs[['Sample Name', 'Group']+new_cols]
passed_AUCs = passed_AUCs.sort_values('Group') # v6.5 export thest to excel
# embed()


''' calcs for stats tab '''


''' ADDING STATS info '''
def get_stats(data):
    targeted_chemicals = len(data.columns)-2

    dt = data[data.columns[2:]].transpose() # transposed data
    data['TML'] = dt.sum() 
    
    data.loc[:,'hits'] = data[data.columns[2:]].transpose().count() # first col is group...


    case_data = data[data[GROUP_COL_NAME]==CASE]
    control_data = data[data[GROUP_COL_NAME]==CONTROL]

    data = data[data[GROUP_COL_NAME].isin([CASE,CONTROL])]


    # import statsmodels.api as sm
    from scipy.stats import fisher_exact, mannwhitneyu, ttest_ind

    from scipy.stats.mstats import gmean

    # organize hit vs no hit data for fishers exact:
    binary = data[data.columns[2:-2]].fillna(0)
    binary[binary>0]=1
    binary = binary.replace(0,'Not')
    binary = binary.replace(1,'Hit')
    binary[GROUP_COL_NAME] = data[GROUP_COL_NAME]

    N = len(data)
    CASE_N = len(case_data)
    CONTROL_N = len(control_data)
    
    ''' calculating stats for stats tab '''

    out = []

    # tml for all aucs
    case_tml = case_data[case_data.columns[2:-2]].sum()
    control_tml = control_data[control_data.columns[2:-2]].sum()

    case_tml = case_tml[case_tml>0]
    control_tml = control_tml[control_tml>0]
    gmean_case_tml = gmean(case_tml)
    gmean_control_tml = gmean(control_tml)

    # case vs control tml
    tml_stat, tml_mannw_pval = stats.mannwhitneyu(case_tml,control_tml, use_continuity=False, alternative='two-sided')


    for c in data.columns[2:-2]:
        
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
            stat, mannw_pval = stats.mannwhitneyu(case,control, use_continuity=False, alternative='two-sided') # NOTE: these parameters will make the MWpval=KWpval (for two samples)
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

    data.index = data['Sample Name']
    del data['Sample Name']
    stats_out = data.transpose()
    stats_out = pd.merge(stats_out, out, left_index=True, right_index=True, how='left')

    del stats_out['MRM Name']
    
    return stats_out

stats_out = get_stats(passed_AUCs)

''' END STATS TAB '''

out2 = out2.style.apply(apply_color, axis=None)



    
outfile = study_name+'_vetted_peaks.xlsx'
writer = pd.ExcelWriter(outfile)

import datetime as dt

readme = {
    'Date of Analysis':	dt.datetime.today().strftime("%m/%d/%Y"),
    'MRM List':'Exogenous',
    'Columns':'Raptor + Shodex',
    'Study Name':study_name,
    'Input File Name':filename,
    'Input File Version Date':'7/22/2019',
    'Chemicals Targeted	':len(compounds),
    'Chemicals Detected (S/N >%i)'%MRM1_vs_BLANK_THRESHOLD:len(case_control_positive)
}

for group in exp_groups:
    readme['N (%s)'%group] = len(group_indices[group])


tmp = pd.DataFrame({'README':readme}).to_excel(writer,'README')
data.to_excel(writer,'original')
blank_avgs.to_excel(writer,'blank_averages')
ratios[data.columns[2:]].to_excel(writer,'ratios')
out[data.columns[2:]].to_excel(writer,'vetted')

out2.to_excel(writer,'vetted2')

passed_AUCs.to_excel(writer, 'Positive MRM1 AUCs')

case_control_positive.to_excel(writer,'count positives')

stats_out.to_excel(writer,'stats')

writer.save()
    
print('DONE, saved to: %s'%outfile)





