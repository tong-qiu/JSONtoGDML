import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pickle
from matplotlib import rc, rcParams
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches
import copy
from itertools import cycle

def errobar_cat(cat, height, labels, **kwargs):
    settings = {
        "xlabel" : "",
        "ylabel": 'Number of Steps',
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        "title2": r"$\mathit{\sqrt{s}=13.6\:TeV}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "",
        "filename": "ratio",
        "log_y":False,
        "upper_y": 1.7, 
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    
    matplotlib.rcParams['font.sans-serif'] = "FreeSans"
    matplotlib.rcParams['font.family'] = "sans-serif"

    cycol = cycle(["#3f90da","#ffa90e","#bd1f01","#94a4a2","#832db6","#a96b59","#e76300","#b9ac70","#717581","#92dadd"])


    fig, ax = plt.subplots(figsize=(10,8))

    datadic = {}
    for i in range(len(labels)):
        datadic[labels[i]] = height[i]

    x = np.arange(len(cat))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0


    for attribute, measurement in datadic.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute, color=next(cycol))
        # ax.bar_label(rects, padding=3)
        multiplier += 1
    ax.set_xticks(x + width, cat)
    plt.xticks(rotation=45, rotation_mode="anchor", ha='right',)
    # for i in range(len(height)):
    #     ax.bar(cat[i], height[i], label=labels[i], color=next(cycol), alpha=0.7, edgecolor='black', linewidth=1.2)
    ax.legend(loc='upper right',prop={'size': 18}, frameon=False)

    ymin, ymax = ax.get_ylim()
    ax.set_ylim([0, ymax* settings["upper_y"]])
    ax.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=25, transform=ax.transAxes)
    ax.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
    ax.set_ylabel(settings['ylabel'], fontsize=20)
    if settings['log_y']:
        ax.set_yscale('log')
        ax.set_ylim([0.0001, 100**(ymax * settings["upper_y"])])
        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,numticks=100))
        ax.minorticks_on()


    fig.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.25)
    plt.close(fig)

def get_cats(dicts):
    cats = []
    for each_dict in dicts:
        for each_key in each_dict:
            cats.append(each_key)
    cats = list(set(cats))

    values = []
    for each_dict in dicts:
        each_values = []
        for each_cat in cats:
            if each_cat in each_dict:
                each_values.append(each_dict[each_cat])
            else:
                each_values.append(0)
        values.append(each_values)
    return cats, values

def processfile(filename):
    steps = {}
    steps_noworld = {}
    with open(filename) as f:
        for each_line in f:
            if "G4Worker_0 >" in each_line:
                components = each_line.split()
                if len(components) == 12:
                    if components[2].isdigit():
                        if components[11] not in steps:
                            steps[components[11]] = 1
                            steps_noworld[components[11]] = 0
                            if "world" not in components[10].lower():
                                steps_noworld[components[11]] = 1
                        else:
                            steps[components[11]] += 1
                            if "world" not in components[10].lower():
                                steps_noworld[components[11]] += 1
    cats, values = get_cats([steps, steps_noworld])
    return steps, steps_noworld

def main():
    steps_TIK, steps_noworld_ITK = processfile("ITK.log")
    steps_sim, steps_noworld_sim = processfile("fsl_simplified.log")
    cats, values = get_cats([steps_TIK, steps_noworld_ITK, steps_sim, steps_noworld_sim])

    errobar_cat(cats, values, ["Full geo", "Full geo no world", "simplified" , "simplified no world"], filename="everything", title3="Steps count for each G4 process")
    errobar_cat(cats, [values[0], values[2]], ["Full geo", "simplified"], filename="full", title3="Steps count for each G4 process")
    errobar_cat(cats, [values[1], values[3]], ["Full geo no world", "simplified no world"], filename="noworld", title3="Steps count for each G4 process")

if __name__ == "__main__":
    main()
