# Any Python script in the scripts folder will be able to import from this module.

# Determine targets for optional rules
def provided(samplelist, condition):
    '''
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    '''
    if not condition:
        samplelist = []
    return samplelist
