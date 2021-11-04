import pandas as pd

prop_color = {0: 'ffe6e6',
              1: 'e6ffff',
              2: 'ddeeff',
              3: 'ddffdd',
              4: 'fff2ff',
              5: 'f9ffe6',
              6: 'e7e6ff',
              7: 'efffe6',
              8: 'fff3e6',
              9: 'b3c7b6',
              10: 'bba9bc',
              11: 'fffbe6',
              12: 'ffe6ec',
              13: 'e6e7ff',
              14: 'd5e7d0',
              15: 'ddcfe6',
              16: 'e9d2de',
              17: 'cbd0e1',
              18: 'f5f5f5'}

def make_inc_color(summary):
    Inc = summary[summary['include?']=='yes']['Inc group']
    Inc_count = Inc.mask(Inc.str.contains(',', na=False)).value_counts()
    print(Inc_count)
    
    color_dict = pd.DataFrame()
    for i in range(len(Inc_count.index)):
        if Inc_count[Inc_count.index[i]] > 1:
            color_dict[i] = [Inc_count.index[i], prop_color[i]]
        else:
            color_dict[i] = [Inc_count.index[i], prop_color[18]]
    
    color_dict[len(Inc_count.index)] = ['Not identified',prop_color[18]]
    return color_dict.T
    
def make_country_color(summary):
    country = summary[summary['include?']=='yes']['country']
    country_colon = country.mask(country.str.contains(':', na=False)).value_counts()
    country_not_colon = country.str.extract(r'(.*?):')[0].value_counts()
    
    country_count = country_colon.add(country_not_colon,fill_value=0).sort_values(ascending=False)
    
    color_dict = pd.DataFrame()
    for i in range(len(country_count.index)):
        if country_count[country_count.index[i]] > 1:
            color_dict[i] = [country_count.index[i], prop_color[i%17]]
        else:
            color_dict[i] = [country_count.index[i], prop_color[18]]

    color_dict[len(country_count.index)] = ['Not identified',prop_color[18]]

    return color_dict.T
    
def make_organism_color(summary):
    organism = summary[summary['include?']=='yes']['organism']
    if len(organism.value_counts())<17:
        organism_count = organism.value_counts()
    else:
        organism_count = organism.str.extract('(.*?)\s')[0].value_counts()
    
    color_dict = pd.DataFrame()
    for i in range(len(organism_count.index)):
        if organism_count[organism_count.index[i]] > 1:
            color_dict[i] = [organism_count.index[i] + ' spp.', prop_color[i%17]]
        else:
            color_dict[i] = [organism_count.index[i] + ' spp.', prop_color[18]]

    color_dict[len(organism_count.index)] = ['Not identified', prop_color[18]]
    
    return color_dict.T

