#---------------------------------------------------------------------------
# pourbaix diagram plot script
#---------------------------------------------------------------------------
# Hark Lee 2016

#---------------------------------------------------------------------------
# title
#---------------------------------------------------------------------------
#elements_or_formulas.sort()
title_el=""
title_comp=" [ "
for i,el in enumerate(elements_or_formulas):
    if i==len(elements_or_formulas)-1:
        title_el+=str(el)
        if comp_dict=={}:
            title_comp+=str(round(1.0/len(elements_or_formulas),2))+" ] Pourbaix Diagram"
        else:
            title_comp+=str(round(comp_dict[el],2))+" ] Pourbaix Diagram"
    else:
        title_el+=str(el)+"-"
        if comp_dict=={}:
            title_comp+=str(round(1.0/len(elements_or_formulas),2))+" "
        else:
            title_comp+=str(round(comp_dict[el],2))+" "

title=title_el+title_comp
limits=[[-2,16],[-3,3]]
label_domains=True
passive=False

#---------------------------------------------------------------------------
# import modules
#---------------------------------------------------------------------------
import platform
import matplotlib
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FormatStrFormatter
import six
from six.moves import map
from six.moves import zip
import numpy as np
import re
import collections
from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer
from pymatgen.analysis.pourbaix.maker import PREFAC
from pymatgen.analysis.pourbaix.entry import MultiEntry
from pymatgen.phasediagram.plotter import uniquelines
from pymatgen.util.string_utils import latexify
from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.util.coord_utils import in_coord_list
from matplotlib.patches import Polygon
from pymatgen import Element
from itertools import chain
import operator

#---------------------------------------------------------------------------
# plot
#---------------------------------------------------------------------------
def latexify_ion(formula):
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", formula)

plt = get_publication_quality_plot(24, 16)

(stable, unstable) = plotter.pourbaix_plot_data(limits)

optim_colors = ['#0000FF', '#FF0000', '#00FF00', '#FFFF00', '#FF00FF',
                '#FF8080', '#DCDCDC', '#800000', '#FF8000']
optim_font_colors = ['#FFC000', '#00FFFF', '#FF00FF', '#0000FF', '#00FF00',
                     '#007F7F', '#232323', '#7FFFFF', '#007FFF']
mark_passive = {key: 0 for key in stable.keys()}

if plotter._pd._elt_comp:
    maxval = max(six.iteritems(plotter._pd._elt_comp), key=operator.itemgetter(1))[1]
    key = [k for k, v in plotter._pd._elt_comp.items() if v == maxval]
#passive_entry = key[0]
passive_entry = 'Mg'

def list_elts(entry):
    elts_list = set()
    if isinstance(entry, MultiEntry):
        for el in chain.from_iterable([[el for el in e.composition.elements] for e in entry.entrylist]):
            elts_list.add(el)
    else:
        elts_list = entry.composition.elements
    return elts_list

for entry in stable:
    #if passive_entry + str("(s)") in entry.name:
    #    mark_passive[entry] = 2
    #    continue
    if "(s)" not in entry.name:
        continue
    #elif len(set([Element("O"), Element("H")]).intersection(set(list_elts(entry)))) > 0:
    else:
        mark_passive[entry] = 1
        for e in entry.entrylist:
            if e.name==passive_entry+str("(s)"):
                mark_passive[entry] = 2
        
if limits:
    xlim = limits[0]
    ylim = limits[1]
else:
    xlim = plotter._analyzer.chempot_limits[0]
    ylim = plotter._analyzer.chempot_limits[1]

h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                       [xlim[1], -xlim[1] * PREFAC]])
o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                       [xlim[1], -xlim[1] * PREFAC + 1.23]])
neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

ax = plt.gca()
ax.set_xlim(xlim)
ax.set_ylim(ylim)

label_list=[]
c_list=[]
fontcolor_list=[]
plot_passive_list=[]
for entry in stable.keys():
    xy = plotter.domain_vertices(entry)
    c = plotter.get_center(stable[entry])
    c_list.append(c)

    if mark_passive[entry] == 1:
        color = optim_colors[0]
        fontcolor = optim_font_colors[0]
        colorfill = True
        plot_passive_list.append(True)
    elif mark_passive[entry] == 2:
        color = optim_colors[1]
        fontcolor = optim_font_colors[1]
        colorfill = True
        plot_passive_list.append(True)
    else:
        color = "w"
        colorfill = False
        fontcolor = "k"
        plot_passive_list.append(False)
    fontcolor_list.append(fontcolor)

    patch = Polygon(xy, facecolor=color, closed=True, lw=3.0, fill=colorfill)
    ax.add_patch(patch)
    if label_domains:
        str_name = ""
        if isinstance(entry, MultiEntry):
            for e in entry.entrylist:
                str_name += latexify_ion(latexify(e.name)) + " + "
            label = str_name[:-3]
        else:
            label = latexify_ion(latexify(entry.name))
        label_list.append(label)
        #plt.annotate(plotter.print_name(e), c, color=fontcolor, fontsize=20)

lw = 3
plt.plot(h_line[0], h_line[1], "r--", linewidth=lw)
plt.plot(o_line[0], o_line[1], "r--", linewidth=lw)
plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

#---------------------------------------------------------------------------
# labels
#---------------------------------------------------------------------------
# labels
counter=0
xy_correction_list=[]
xshift=0.04
for xy,label in zip(c_list,label_list):
    #plt.annotate(counter, xy, fontsize=15, color="b")
    print str(counter)+" - "+str(label)+" - "+str(len(label))
    counter+=1
    xy_correction_list.append([-xshift*len(label), 0])

# apply label corrections

#xy_correction_list[7][0]+=1


counter=0
visual_counter=1
for xy,xy_correction,label,fontcolor,plot_passive in zip(c_list,xy_correction_list,label_list,fontcolor_list,plot_passive_list):
    x=xy[0]+xy_correction[0]
    y=xy[1]+xy_correction[1]
    counter_plot=counter
    
    if(passive):
        if(plot_passive):
            plt.annotate(label+" ["+str(counter_plot)+"]",[x,y], fontsize=15, color=fontcolor)
            #plt.annotate(label,[x,y], fontsize=15, color="b")    
            visual_counter+=1
    else:
        plt.annotate(label+" ["+str(counter)+"]",[x,y], fontsize=15, color=fontcolor)
        
    counter+=1

#---------------------------------------------------------------------------
# plot settings
#---------------------------------------------------------------------------
plt.xlabel("pH", fontsize=30)
plt.ylabel("E (V)", fontsize=30)
t=plt.title(title, fontsize=30)
t.set_y(1.02)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

#---------------------------------------------------------------------------
# display plot
#---------------------------------------------------------------------------
plt.show()
