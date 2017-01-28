'''
(c) 2017 Alex Pai and Daniel Weiner

pTDT is a command line tool tool for analyzing the 
transmission of polygenic risk within families.


'''

import argparse 
import pandas as pd
import numpy as np
import scipy.stats as ss
import socket
import os
import time
import ptdt

try:
    x = pd.DataFrame({'A': ['1', '2', '3']})
    x = x.apply(pd.to_numeric, errors='coerce')
except AttributeError:
    raise ImportError('ptdt requires pandas version > 0.17.0')

    
def required_length(nmin,nmax):
    '''
    Raises error through argparse if length of optional
    arguments for an option is not within specified bounds.
    
    '''
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(f=self.dest,nmin=nmin,nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength

def parse(args):
    '''
    Takes arguments from parser and feeds into correct variables.
    
    '''
    ### Indices to parse PRS file, depending on column of ID and PRS 
    if len(args.prs) == 1:
        PRS = open(args.prs[0])
        PRS_iid = 1
        PRS_prs = 3
    elif len(args.prs) == 3:
        PRS = open(args.prs[0])
        PRS_iid = int(args.prs[1]) - 1
        PRS_prs = int(args.prs[2]) - 1
    else:
        print('Error: --prs incorrect number of arguments called')

    structure = args.structure
    subset = args.subset
    quad = args.quad
    table = args.table
    outname = args.out
    
    return PRS, PRS_iid, PRS_prs, structure, subset, quad, table, outname

def ptdt_analysis(PRS, PRS_iid, PRS_prs, structured):
    '''
    Takes PRS score file and strcture matrix and returns
    pTDT summary statistics.
    
    '''

    PRS_values = {}
    for line in PRS:
        try:
            PRS_values[line.split()[int(PRS_iid)]] = float(line.split()[int(PRS_prs)]) 
            #creating dict for each ID to call its PRS
        except ValueError:
            continue
            
    print('Creating pTDT matrix.', end="") #progress
    data = [] #initialize list to store each row of matched matrix

    for list_ in structured:
        row = [] #create empty list per row
        values = list_
        row.append(values[0])
        n = 4
        if quad == 'True': #Include sibling PRS if quad flag called
            n += 1
        for i in range(1,n): #fill each place in list with values from dict
            try:
                row.append(str(PRS_values[values[i]]))

            except KeyError: 
                row.append('NA') #return NA for missing PRS values
        data.append(row)

    print('.', end="")
    if quad == 'True': #Create a sibling column if quad flag called
        matrix = pd.DataFrame(data, columns=['ID','Proband','Father','Mother','Sibling'])
    else: 
        matrix = pd.DataFrame(data, columns=['ID','Proband','Father','Mother'])
    # convert dict values to numeric and force NaNs from missing
    IDs = matrix['ID']
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    matrix['ID'] = IDs #rescue non numeric IDs

    print('.', end="")

    # create output matrix
    output = pd.DataFrame(index = range(0, matrix.shape[0]), columns = ["FID","mp_PRS","pro_PRS","pro_pTDT"])

    output['FID'] = matrix[[0]]
    output['mp_PRS'] = (matrix.iloc[:,2]+ matrix.iloc[:,3])/2
    output['pro_PRS'] = matrix[[1]]
    sd = np.std(output['mp_PRS'], ddof=1)
    output['pro_pTDT'] = (output.iloc[:,2] - output.iloc[:,1])/sd
    if quad == 'True':
        output['sib_PRS'] = matrix[[4]]
        output['sib_pTDT'] = (output.iloc[:,4] - output.iloc[:,1])/sd
    print(' done.')

    # Quality control
    output1 = output.dropna(subset=['pro_PRS','mp_PRS'])
    corr1 = ss.pearsonr(output1['pro_PRS'],output1['mp_PRS'])[0]
    if quad == 'True':
        output2 = output.dropna(subset=['sib_PRS','mp_PRS'])
        corr2 = ss.pearsonr(output2['sib_PRS'],output2['mp_PRS'])[0]
    else:
        corr2 = 1
    if corr1 >= .2 and corr2 >= .2:
        print('QC pass.')
        log.write('QC pass.\n')
    else:
        print('WARNING: QC fail - Low correlation between mid-parent PRS and proband/sibling PRS.')
        log.write('WARNING: QC fail - Low correlation between mid-parent PRS and proband/sibling PRS.\n')
    
    # t-test
    x = output['pro_pTDT']
    x_totlength = len(x)
    x = x[~np.isnan(x)] # remove NaNs 
    x_usedlength = len(x)
    print('{0} probands used in pTDT analysis ({1} skipped due to missingness).' .format(x_usedlength,x_totlength-x_usedlength))
    log.write('{0} probands used in pTDT analysis ({1} skipped due to missingness).\n' .format(x_usedlength,x_totlength-x_usedlength))
    if x_usedlength == 0:
        raise RuntimeError('Column order in PRS file is incorrect')
    pTDT_mean = np.mean(x)
    pTDT_std = np.std(x,ddof=1)/np.sqrt(len(x))
    pTDT_pvalue = ss.ttest_1samp(x, popmean=0).pvalue

    values = [pTDT_mean,pTDT_std,pTDT_pvalue]

    if quad == 'True':
        y = output['sib_pTDT']
        y_totlength = len(y)
        y = y[~np.isnan(y)]
        y_usedlength = len(y)
        print('{0} siblings used in pTDT analysis ({1} skipped due to missingness).' .format(y_usedlength,y_totlength-y_usedlength))
        log.write('{0} siblings used in pTDT analysis ({1} skipped due to missingness).\n' .format(y_usedlength,y_totlength-y_usedlength))
        unaffected_mean = np.mean(y)
        unaffected_std = np.std(y,ddof=1)/np.sqrt(len(y))
        unaffected_pvalue = ss.ttest_1samp(y, popmean=0).pvalue
        values2 = [unaffected_mean,unaffected_std,unaffected_pvalue]
    
        values = values + values2
    
    if table == 'True':
        output['dad_PRS'] = matrix[[2]]
        output['mom_PRS'] = matrix[[3]]
        if quad == 'True':
            cols = ['FID','pro_PRS','sib_PRS','dad_PRS','mom_PRS','mp_PRS','pro_pTDT','sib_pTDT']
        else:
            cols = ['FID','pro_PRS','dad_PRS','mom_PRS','mp_PRS','pro_pTDT']
        output = output[cols]
        output = output.round(4)
        headers = list(output.columns.values)
        outfile = open(outname+'.ptdt.table','w+')
        outfile.flush()
        output.to_csv(outname+'.ptdt.table', na_rep = 'NA', header=headers, index=None, sep='\t', mode='a+')

    return values


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''
pTDT Script v1.0.0 (26 Jan 2017)
(c) 2017 Alex Pai and Daniel Weiner
Runs pTDT analysis for a specified family structure
See <https://github.com/ypaialex/ptdt> for more details''')
parser.add_argument('--prs', nargs ='+', help='\x1b[A\b\b\b\b\b\bX Y] \x1b[1CPRS score file [Family ID] [Individual ID] [PHENO] [PRS], X = [Individual ID] column number, Y = [PRS] column number, default X Y = 2 4',
                    metavar=('FILENAME', ''),required=True)
parser.add_argument('--structure', help='Family structure file (in this order): [Family ID] [Proband_IID] [Father_IID] [Mother_IID] ([Sibling_IID])', metavar='FILENAME',type=argparse.FileType('r'),required=True)
parser.add_argument('--subset',metavar='FILENAME',type=argparse.FileType('r'), 
                    help='Add subset of families to analyze, header optional',default=0)
parser.add_argument('--quad',help='Use if Family structure file includes unaffected siblings', 
                    action='store_const', const='True',default='False')
parser.add_argument('--print',dest='table',help='Prints table of pTDT values',action='store_const', const='True',default='False')
parser.add_argument('--out',help='Name of output files',required=True)
args = parser.parse_args()

if __name__ == "__main__":
    PRS, PRS_iid, PRS_prs, structure, subset, quad, table, outname = parse(args)

    print('pTDT Script v1.0.0 (26 Jan 2017)')
    print('(c) 2017 Alex Pai and Daniel Weiner')
    print('')
    print('Writing log file to '+outname+'.ptdt.log')
    print('Options invoked:')
    print('\t--PRS {0}' .format(PRS.name))
    print('\t--structure {0}' .format(structure.name))
    if subset != 0:
        print('\t--subset {0}' .format(subset.name))
    if quad == 'True':
        print('\t--quad')
    if table == 'True':
        print('\t--print')
    print('\t--out '+outname+'\n')
    print('')

    #log file
    log = open(outname+'.ptdt.log', 'w')
    log.write('pTDT Script v1.0.0 (26 Jan 2017)\n')
    log.write('Options invoked:\n')
    log.write('\t--PRS {0}\n' .format(PRS.name))
    log.write('\t--structure {0}\n' .format(structure.name))
    if subset != 0:
        log.write('\t--subset {0}\n' .format(subset.name))
    if quad == 'True':
        log.write('\t--quad\n')
    if table == 'True':
        log.write('\t--print\n')
    log.write('\t--out '+outname+'\n')
    log.write('\n')
    log.write('Hostname: {0}\n' .format(socket.gethostname()))
    name = os.getcwd()
    log.write('Working directory: {0}\n' .format(name))
    log.write('Start time: {0}\n\n' .format(time.strftime("%c")))

    if quad == 'True':
        first_line = structure.readline()
        try:
            sib_true = first_line.split()[4]
        except IndexError:
            raise RuntimeError('--quad invoked without full family structure')
    else:
        first_line = structure.readline()
        pass

    
    structured, miscount = ptdt.subset(args, log)
    values = ptdt_analysis(PRS, PRS_iid, PRS_prs, structured)

    #results
    print('')
    print('------------------------------------------')
    print('Proband analysis')
    print('pTDT mean: %.3E SD' % values[0])
    print('pTDT SE: %.3E SD' % values[1])
    print('pTDT pvalue: %.3E' % values[2])
    if quad == 'True':
        print('')
        print('Sibling analysis')
        print('pTDT mean: %.3E SD' % values[3])
        print('pTDT SE: %.3E SD' % values[4])
        print('pTDT pvalue: %.3E' % values[5])
    print('------------------------------------------')
    print('')

    if miscount >= 1:
        print('--subset: Missing families in structure file printed to ' + outname+'.missing')
        log.write('--subset: Missing families in structure file printed to ' + outname+'.missing\n')
    if table == 'True':
        log.write('--print: Table written to ' +outname+'.ptdt.table\n')
        print('--print: Table written to ' +outname+'.ptdt.table')
    results = open(outname+'.ptdt.result', 'w')
    results.write('Proband analysis\n')
    results.write('pTDT mean: %.3E SD\n' % values[0])
    results.write('pTDT SE: %.3E SD\n' % values[1])
    results.write('pTDT pvalue: %.3E\n' % values[2])
    if quad == 'True':
        results.write('\n')
        results.write('Sibiling analysis\n')
        results.write('pTDT mean: %.3E SD\n' % values[3])
        results.write('pTDT SE: %.3E SD\n' % values[4])
        results.write('pTDT pvalue: %.3E\n' % values[5])
    print('--output: Results written to '+outname+'.ptdt.result')
    log.write('--output: Results written to '+outname+'.ptdt.result\n\n')
    log.write('End time: {0}' .format(time.strftime("%c")))
    log.flush()
    results.flush()

