'''
(c) 2017 Alex Pai, Daniel Weiner, and Elise Robinson

Iterates through structure file and excludes families
found in subset file, if given.


'''

def subset(args, log):
    
    structure = args.structure
    subset = args.subset
    quad = args.quad
    outname = args.out
    
    if subset != 0:
        s_dict = {}
        structured = []
        count = 0
        for line in structure:
            count += 1
            fields = line.split()
            n = 4
            if quad == 'True': #Include extra sibling column if quad flag called
                n += 1
            s_dict[fields[0]] = fields[1:n] #create list of individual ID values based on family ID
        print('{0} families loaded from structure file.' .format(count))
        log.write('{0} families loaded from structure file.\n' .format(count))
        #now call dict values from subset values        
        count1 = 0
        miscount = 0
        missing = []
        for line in subset:
            count1 += 1
            num = line.split()[0]
            try:
                row = [num]+s_dict[str(num)]
                structured.append(row)
            except KeyError:
                miscount += 1
                missing.append(num)
                continue
        print('{0} families loaded from subset file ({1} families excluded from analysis).' .format(count1-miscount,count-count1+miscount))
        log.write('{0} families loaded from subset file ({1} families excluded from analysis).\n' .format(count1-miscount,count-count1+miscount))
        if miscount >= 1:
            missinglog = open(outname+'.missing', 'w+')
            print('\tWarning: {0} families not found in structure file' .format(miscount))
            log.write('\tWarning: {0} families not found in structure file.\n' .format(miscount))
            for item in missing:
                missinglog.write('{0}\n' .format(item))
            missinglog.flush()
    
    elif subset == 0:
        structured = []
        count = 0
        miscount = 0
        for line in structure:
            count += 1
            fields = line.split()
            n = 4
            if quad == 'True': #Include extra sibling column if quad flag called
                n += 1
            structured.append(fields[0:n])
        print('{0} families loaded from structure file.' .format(count))
        log.write('{0} families loaded from structure file.\n' .format(count))
        
    return structured, miscount

