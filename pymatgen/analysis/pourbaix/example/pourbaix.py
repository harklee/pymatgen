#---------------------------------------------------------------------------
# pourbaix diagram script
# (in conjunction with plot.py)
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# import modules
#---------------------------------------------------------------------------
# warnings
import warnings

# I/O with databases
from matgendb import QueryEngine
from pymatgen.matproj.rest import MPRester

# element, composition
from pymatgen import Composition, Element

# phase diagrams
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

# pourbaix diagrams
from pymatgen.analysis.pourbaix.entry import PourbaixEntry
from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.plotter import PourbaixPlotter
#from pymatgen.analysis.pourbaix.analyzer import PourbaixAnalyzer

# compatibility
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility, MaterialsProjectCompatibility, AqueousCorrection

# ion / reference solid dictionaries 
from pymatpro.pourbaix_tools.IonIO import IonIO
from pymatpro.pourbaix_tools.reference_energies import IonReference


#---------------------------------------------------------------------------
# elements and compositions
#---------------------------------------------------------------------------
elements_or_formulas=["Mg","Sn","Zr","Sb"]
#comp_dict={"Mg":0.9,"Ca":0.05,"Sr":0.05}
comp_dict={'Mg':0.90,'Sn':0.04,'Zr':0.03,'Sb':0.03}
conc_dict={}

#---------------------------------------------------------------------------
# data switches
#---------------------------------------------------------------------------
inter=True
ss=True

#---------------------------------------------------------------------------
# hydride screen switch
#---------------------------------------------------------------------------
nohydride=True

#---------------------------------------------------------------------------
# Materials Project database Rest API
#---------------------------------------------------------------------------
mpr=MPRester(api_key="")

#---------------------------------------------------------------------------
# intermetallic database
#---------------------------------------------------------------------------
if inter:
    qe_inter=QueryEngine(host="mongodb03.nersc.gov",user="admin_hark_inter",password="",database="vasp_hark_inter",port=27017,collection="tasks")

#---------------------------------------------------------------------------
# solid solution oxide database
#---------------------------------------------------------------------------
if ss:
    qe_oxide=QueryEngine(host="mongodb03.nersc.gov",user="admin_hark_mg",password="",database="vasp_hark_mg",port=27017,collection="tasks")

#---------------------------------------------------------------------------
# solid solution metal database
#---------------------------------------------------------------------------
if ss:
    qe_metal=QueryEngine(host="mongodb03.nersc.gov",user="admin_hark_vasp",password="",database="vasp_hark_vasp",port=27017,collection="tasks")

#---------------------------------------------------------------------------
# compatibility module
#---------------------------------------------------------------------------
compat=MaterialsProjectAqueousCompatibility("Advanced")
#compat=MaterialsProjectCompatibility("Advanced",correct_peroxide=True)
#aqcompat = AqueousCorrection("/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/pymatgen-3.4.0-py2.7-macosx-10.6-intel.egg/pymatgen/entries/MPCompatibility.yaml")

#---------------------------------------------------------------------------
# elements
#---------------------------------------------------------------------------
elements = set()
pbx_elts = list(elements_or_formulas)
for formula in pbx_elts:
    try:
        comp = Composition(formula)
        elements.update(comp.elements)
    except Exception:
        raise ValueError(str(formula) + " is not a valid element or formula!")

elements_or_formulas = [elf.symbol for elf in elements if elf.symbol not in ["O", "H", ""]]

#---------------------------------------------------------------------------
# function to check duplicates
#---------------------------------------------------------------------------
def contains_entry(entry_list, entry):
    for e in entry_list:
        if e.entry_id == entry.entry_id or (abs(entry.energy_per_atom - e.energy_per_atom) < 1e-6 and entry.composition.reduced_formula == e.composition.reduced_formula):
            return True

#---------------------------------------------------------------------------
# function to check duplicates
#---------------------------------------------------------------------------
def hydride(x):
    el=x.composition.elements
    if Element('H') in el:
        if len(el)>=2:
            if Element('O') in el:
                return False
            else:
                return True
        else:
            return False
    else:
        return False

#---------------------------------------------------------------------------
# elements + O + H
#---------------------------------------------------------------------------
elements_to_query = list(elements_or_formulas)
if "O" not in elements_to_query:
    elements_to_query.append("O")
if "H" not in elements_to_query:
    elements_to_query.append("H")

#---------------------------------------------------------------------------
# load intermetallic database
#---------------------------------------------------------------------------
#results = qe_inter.get_entries_in_system(elements_to_query, optional_data=["oxide_type"])
print "----------------------------------------------------------------------------"
print " Accessing Intermetallic Database"
print "----------------------------------------------------------------------------"
entries_inter=qe_inter.get_entries_in_system(elements_to_query,inc_structure=False,optional_data=["oxide_type"],additional_criteria={"vaspinputset_name":"DictVaspInputSet"})
print "----------------------------------------------------------------------------"
print " Number of Intermetallics = ",len(entries_inter)
print "----------------------------------------------------------------------------"
print " "
    
#---------------------------------------------------------------------------
# load solid solution oxide database
#---------------------------------------------------------------------------
print "----------------------------------------------------------------------------"
print " Accessing SS Oxide (hark_mg) Database"
print "----------------------------------------------------------------------------"
entries_oxide=qe_oxide.get_entries_in_system(elements_to_query,inc_structure=False,optional_data=["oxide_type"],additional_criteria={"vaspinputset_name":"DictVaspInputSet","nelements":4})
#entries_oxide=qe_metal.get_entries_in_system(elements_to_query,inc_structure=True,optional_data=["oxide_type"],additional_criteria={"vaspinputset_name":"DictVaspInputSet","nelements":4})
print "----------------------------------------------------------------------------"
print " Number of SS Oxides = ",len(entries_oxide)
print "----------------------------------------------------------------------------"
print " "
        
#---------------------------------------------------------------------------
# load solid solution metal database
#---------------------------------------------------------------------------
print "----------------------------------------------------------------------------"
print " Accessing SS Metal (hark_vasp) Database"
print "----------------------------------------------------------------------------"
entries_metal=qe_metal.get_entries_in_system(elements_to_query,inc_structure=False,optional_data=["oxide_type"],additional_criteria={"vaspinputset_name":"DictVaspInputSet","nelements":3})
print "----------------------------------------------------------------------------"
print " Number of SS Metals = ",len(entries_metal)
print "----------------------------------------------------------------------------"
print " "

#---------------------------------------------------------------------------
# load materials project database
#---------------------------------------------------------------------------
entries_mp = mpr.get_entries_in_chemsys(elements_to_query)
print "----------------------------------------------------------------------------"
print " Number of MP Compounds = ",len(entries_mp)
print "----------------------------------------------------------------------------"
print " "

#---------------------------------------------------------------------------
# combine entries
#---------------------------------------------------------------------------
results=list()
for x in entries_mp:
    results.append(x)
if(inter):
    for y in entries_inter:
        results.append(y)
if(ss):
    for z in entries_metal:
        results.append(z)
    for z in entries_oxide:
        results.append(z)

#---------------------------------------------------------------------------
# apply corrections (at once)
#---------------------------------------------------------------------------
#results_compat = compat.process_entries(results)

#---------------------------------------------------------------------------
# duplicate check and apply corrections (individually)
#---------------------------------------------------------------------------
results_unique=list()
for entry in results:
    entry_compat = compat.process_entry(entry)
    if not contains_entry(results_unique,entry_compat):
        results_unique.append(entry_compat)
print "----------------------------------------------------------------------------"
print " After Checking Duplicates = ",len(results_unique)
print "----------------------------------------------------------------------------"
print " "

#---------------------------------------------------------------------------
# remove hydrides
#---------------------------------------------------------------------------p
print "----------------------------------------------------------------------------"
print " Removing Hydrides"
print "----------------------------------------------------------------------------"
print " Number Before Removal = ",len(results_unique)
print "----------------------------------------------------------------------------"
results_final=list()
if nohydride:
    for entry in results_unique:
        if hydride(entry):
            print entry.composition.reduced_formula+"  -  "+entry.entry_id
        else:
            results_final.append(entry)
else:
    results_final=results_unique
    
print "----------------------------------------------------------------------------"
print " Number After Removal = ",len(results_final)
print "----------------------------------------------------------------------------"            
print " "

#---------------------------------------------------------------------------
# generate phase diagram and find stable solids
#---------------------------------------------------------------------------
pd = PhaseDiagram(results_final)

if(len(elements_or_formulas)<3):
    pdp = PDPlotter(pd)
    pdp.get_plot()

stable_solids = pd.stable_entries
stable_solids_minus_h2o = [entry for entry in stable_solids if 
                           entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]

#---------------------------------------------------------------------------
# construct solid entries
#---------------------------------------------------------------------------
solid_entries = list()
for entry in stable_solids_minus_h2o:
    pbx_entry = PourbaixEntry(entry)
    pbx_entry.g0_replace(pd.get_form_energy(entry))
        
    pbx_entry.reduced_entry()
    solid_entries.append(pbx_entry)

#---------------------------------------------------------------------------
# construct ion entries
#---------------------------------------------------------------------------
ion_entries = []
ion_io = IonIO()
all_in = lambda a, b: bool(set(a) <= set(b))
experimental_data_lacking_element = []
for el_f in elements_or_formulas:
    # if the ion dictionary doesn't have the element
    ions = ion_io.get_ions(el_f)
    if len(ions) == 0:
        experimental_data_lacking_element.append(el_f)

    # find reference solids and check if they have other elements
    ref_solid_elt_list = set()
    for ion in ions:
        for elt in Composition(ion.ref_solid).elements:
            ref_solid_elt_list.add(elt.symbol)
    if all_in(ref_solid_elt_list, elements_to_query):
        pd_ion = pd
    else:
        print "Error : Missing Terminal Elements"
        print "Reference Solid List : "+ref_solid_elt_list
        print "Elements to Query : "+elements_to_query

    # set ion concentrations
    if el_f in conc_dict:
        set_conc = conc_dict[el_f]
    else:
        set_conc = 1.0e-8

    # construct ion entries
    for ion in ions:
        ir = IonReference(ion.el[0], ion.ion, ion.exp_energy, ion.ref_solid, ion.exp_ref_solid_energy, pd_ion, ion.name)
        pbx_entry = ir.ion_pourbaix_entry
        pbx_entry.entry_id = ion.source
        pbx_entry.conc = set_conc
        ion_entries.append(pbx_entry)

# if any missing ions
if len(experimental_data_lacking_element) >= 1:
    eles = experimental_data_lacking_element
    message1 = "element" if len(eles) == 1 else "elements"
    message2 = eles[0] if len(eles) == 1 else "{} and {}".format(", ".join(eles[:-1]), eles[-1])
    warnings.warn("No aqueous ion experimental data for {} {}.".format(message1, message2))

#---------------------------------------------------------------------------
# all entries for Pourbaix diagrams
#---------------------------------------------------------------------------
all_entries = solid_entries + ion_entries

#---------------------------------------------------------------------------
# automatically set compositions if not given
#---------------------------------------------------------------------------
elements_or_formulas.sort(reverse=True)
if (len(elements_or_formulas) > 1) and (comp_dict == {}):
    comp = [1.0 / len(elements_or_formulas)] * len(elements_or_formulas)
    comp_dict = dict(zip(elements_or_formulas, comp))

#---------------------------------------------------------------------------
# construct Pourbaix diagram
#---------------------------------------------------------------------------
pbx = PourbaixDiagram(all_entries, comp_dict)

#---------------------------------------------------------------------------
# plot Pourbaix diagram
#---------------------------------------------------------------------------
plotter = PourbaixPlotter(pbx)
execfile("plot.py")

