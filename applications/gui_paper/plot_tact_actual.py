#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

import countrydata

# XXX path to folder with output from `request_country.py`
datafolder = "."
df = countrydata.CollectCountryData(datafolder)
countrydata.AppendColor(df)
countrydata.AppendOfficalLockdown(df)
countrydata.AppendInferred(df, datafolder)

fig, ax = plt.subplots(1, 1, figsize=(5, 5.3))

yticks_y = []
yticks_label = []

myFmt = mdates.DateFormatter('%b %d')
ax.xaxis.set_major_formatter(myFmt)
ax.get_yaxis().set_visible(False)

for i,row in enumerate(df.sort_values(by='tint_mean').itertuples()):
    y = -i
    color = row.color
    official = row.official_lockdown
    ax.scatter(row.tint_mean,
               y,
               s=16,
               marker='o',
               facecolor='none' if pd.isnull(official) else color,
               edgecolor=color)
    ax.plot([row.tint_low, row.tint_high], [y, y],
            c="black",
            lw=3,
            alpha=0.15,
            zorder=-10)

    yticks_y.append(y)
    yticks_label.append(row.displayname)

    if not pd.isnull(official):
        ax.scatter(official, y, s=16, c=color, marker='|')
        ax.plot([official, row.tint_mean], [y, y], c=color)
    ax.annotate(row.displayname,
                xy=(row.tint_mean, y),
                fontsize=7,
                xytext=(15, 0),
                textcoords='offset points',
                va='center')

ax.set_xticks(
    list(
        map(pd.to_datetime, [
            "2020-03-01",
            "2020-03-15",
            "2020-04-01",
            "2020-04-15",
            "2020-05-01",
            "2020-05-15",
        ])))

ax.set_yticks(yticks_y)
ax.set_yticklabels(yticks_label, fontsize=7)

fig.tight_layout()
fig.savefig("scatter_tint_actual.pdf")

