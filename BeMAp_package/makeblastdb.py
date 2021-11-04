import os

gene = {1: 'aminoglycoside',
 2: 'beta-lactam',
 3: 'quinolone',
 4: 'macrolide',
 5: 'phenicol',
 6: 'tetracycline',
 7: 'trimethoprim',
 8: 'sulphonamide',
 9: 'rifampin',
 10: 'fosfomycin',
 11: 'polymxycin',
 12: 'nitrofuran'}

for i in gene.keys():
    AMR_db = 'makeblastdb -in database/ResFinder_db_detail/' + gene[i] + '_detail.fsa -dbtype nucl'
    os.system(AMR_db)

Inc_db = 'makeblastdb -in database/PlasmidFinder_db/enterobacteriaceae.fsa -dbtype nucl'
os.system(Inc_db)
