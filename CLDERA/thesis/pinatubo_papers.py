import pdb
import numpy as np
import pandas as pd
from habanero import Crossref
import matplotlib.pyplot as plt

# Initialize Crossref API client
cr = Crossref()

def search_agu_articles(query, max_results=50):

    total_papers = []
    all_papers   = []
    dates        = []
    counts       = []
    years = np.arange(1991,2024,1)
    for year in years:
        print(f'searching {year}...')
        #dd = cr.works(query='pinatubo', filter={'member':'13',
        filters = {'has_abstract':True, 
                   'from-print-pub-date':'{}-01-01'.format(year), 
                   'until-print-pub-date':'{}-01-01'.format(year+1)}

        # --- get total count of articles
        totals = cr.works(query='climat', filter=filters, limit=1)['message']['total-results']
        totals += cr.works(query='atomspher', filter=filters, limit=1)['message']['total-results']
        total_papers.append(totals)

        # --- get results for query (up to 1000...)
        dd = cr.works(query=query, filter=filters, limit=1000,
                      select=['title', 'issued','reference', 'abstract'])
        items = dd['message']['items']
        abstracts = [i.pop('abstract') for i in items]
        mask = [(ai.find('climat')>-1 or ai.find('atmospher')>-1) for ai in abstracts]
        items = np.array(items)[mask]
        print('found {} papers'.format(len(items)))

        # --- gather returns
        counts.append(len(items))
        all_papers.extend(items)
        dates.extend([i['issued']['date-parts'][0] for i in items])

    return all_papers, dates, total_papers, counts, years

# Run the function with a search query
if(0):
    papers, dates, tot, counts, years = search_agu_articles('pinatubo')
    counts = np.array(counts) / (np.array(tot)/max(tot))
    plt.plot(years, counts, '-sk')
    plt.title('Pinatubo')

    plt.figure()
    papers, dates, tot, counts, years = search_agu_articles('volcanic forcing')
    counts = np.array(counts) / (np.array(tot)/max(tot))
    plt.plot(years, counts, '-sk')
    plt.title('Volcanic Forcing')

plt.figure()
papers, dates, tot, counts, years = search_agu_articles('volcanic aerosol')
counts = np.array(counts) / (np.array(tot)/max(tot))
plt.plot(years, counts, '-sk')
plt.title('Volcanic Aerosol')
plt.show()


