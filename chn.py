# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:10:23 2019
CHN Program for Formula Weights and Elemental Analysis
@author: Mark Runyan
"""
from decimal import Decimal
from datetime import datetime
from atomicWeightsDecimal import atomicWeightsDecimal
from prettytable import PrettyTable

def getList(dict):
    """return list of keys in dicionary """
    return dict.keys()

class Sample():
    """Chemcial Sample"""
    def __init__(self, name, mw=0, formula={}, exactmass=0, hydrate=0):
        self.name = name
        self.mw = mw
        self.formula = {}
        self.exactmass = exactmass
        self.hydrate = hydrate
        
    def input_sample_name(self):
        sample_name = input("Sample Name: ")
        return(sample_name)

    def calc_exactmass(self, formula):
        """Cacluate Monoisotopic Mass from formula"""
        em = 0
        for element in formula:
            em = em + formula[element]*atomicWeightsDecimal[element]["abundant"]
        return(em)
        
    def calc_mw(self, formula):
        """Cacluate Molecular Weight from formula"""
        mw = 0
        for element in formula:
            mw = mw + Decimal(formula[element])*atomicWeightsDecimal[element]["standard"]
        return(mw)
    
    def calc_percents(self, formula, mw, pct={}):
        """Calculate elemental percentages from formula and mw or fw"""
        if "C" in pct:
            pass
        else:
            # Use formula keys to generate nested dictionary for theoretical and experimental results 
            elements = getList(formula)
            pctcat = {'theoretical':'','experimental':''}
            # create nested dictionary
            pct = {k: dict.fromkeys(pctcat,0) for k in elements}
        # calcualate theoretical elemental percentages
        for element in pct:
            pctx =  (Decimal(formula[element])*atomicWeightsDecimal[element]["standard"]/mw)*100
            pct[element]['theoretical'] = "{0:0.2f}".format(pctx)
        return(pct)

    def create_hydrate_formula(self, formula, H2Ox):
        # new dictionary for hydrate formula
        hydrate_formula = formula.copy()
        if "O" in hydrate_formula:
            # udpate
            hydrate_formula["O"] = hydrate_formula.get("O") + H2Ox
        else:
            hydrate_formula["O"] = H2Ox
        if "H" in hydrate_formula:
            # udpate
            hydrate_formula["H"] = hydrate_formula.get("H") + H2Ox*2
        else:
            hydrate_formula["H"] = H2Ox*2
        return(hydrate_formula)

    def format_formula(self, formula):
        """Extract MF from dictionary for output"""
        formula_string = ""
        for e in formula:
            # for elements with one atom, remove number 1
            formula_string = formula_string + e + str.replace(str(formula[e]),"1","",10) + " "
        return(formula_string)
        
    def print_chn_data(self, name, formula, elem_percents, em, mw, fw, hydrate=0, filename=None):
        """ Print CHN Analysis Results """
        tab = PrettyTable()
        header = []
        now = datetime.now()
        gxp_date = now.strftime("%d-%b-%Y %I:%M %p")
        rpt_title = "CHN Analysis Report - "+gxp_date
        header.append(rpt_title)
        header.append(name)
        if hydrate == 0:
            header.append(self.format_formula(formula))
            header.append("MW: "+"{0:0.2f}".format(mw))
        else:
            header.append(self.format_formula(formula)+" * "+ str(hydrate)+"H2O")
            header.append("MW: "+"{0:0.2f}".format(mw)+" (free)")
            header.append("FW: "+"{0:0.2f}".format(fw))
        header.append("Monoisotopic Mass: "+"{0:0.5f}".format(em))
        
        if elem_percents['C']['experimental'] > 0:
            tab.field_names = ["Element","Theoretical Percentage","Experimental","Difference"]
            for element in elem_percents:
                
                if elem_percents[element]['experimental'] > 0:
                    diff = float(elem_percents[element]['experimental'])-float(elem_percents[element]['theoretical'])
                    diff = "{0:0.2f}".format(diff)
                    tab.add_row([element,elem_percents[element]['theoretical'],elem_percents[element]['experimental'],diff])
                else:
                    tab.add_row([element,elem_percents[element]['theoretical'],'',''])
        else:
            tab.field_names = ["Element","Theoretical Percentage"]
            for element in elem_percents:
                row_list = [element,elem_percents[element]['theoretical']]
                tab.add_row(row_list)
        if filename:
            f = open(filename, "a")
            print("", file=f)
            print("\n".join(header), file=f)
            print(tab, file=f)
            f.close()
        else:
            print("")
            print("\n".join(header))
            print(tab)
            
    def input_mf(self):
        common_elements = {'C':'Carbons','H':'Hydrogens','N':'Nitrogens','O':'Oxygens','Cl':'Chlorines','F':'Fluorines'}
        print("Molecular formula Input")
        # formula = {} replace with class instance
        for ce in common_elements:
            while True:
                try:
                    ce_n = input("How many "+common_elements.get(ce) + " ? ")
                    if not ce_n:
                        ce_n = 0
                    ce_n = int(ce_n)
                except ValueError:
                    print("Please input a number, or leave blank")
                    continue
                self.formula[ce] = ce_n
                break
        #print("Common elements : "+str(formula))    
        res = input("Other Elements? [N] ")
        res = res or 'N'
        all_elements = getList(atomicWeightsDecimal)
        # might want to remove common elements already entered
        while res.upper() == 'Y':
            symb = input("Please input element symbol : ")
            if not symb:
                res = input("You entered nothing.  Are you finished? [Y]")
                res = res or 'Y'
                if res == 'Y':
                    res = 'N' # reversed logic, finished do not continue
                    continue        
                if res != 'Y': # not finished continue
                    res = 'Y'
                    continue
            if symb in all_elements:
                while res == 'Y':
                    try:
                        symb_n = input("How many "+symb+" ? ")
                        self.formula[symb] = int(symb_n)
                    except ValueError:
                        print("Please input a number.")
                        continue
                    break
                res = input("More Elements? [Y] ")
                res = res or 'Y'
                continue
            else:
                print("Symbol "+symb+" not recognized")
                continue
        # remove elements from dictionary with n of zero
        self.formula = {e:n for e,n in self.formula.items() if n!=0}
        return(self.formula)

    def print_menu(self):
        print("")
        print("1)Enter new MF\n2)Enter experimental data\n3)Add Water of hydration\n4)Print to file\n5)Exit Program")
# end of class Sample
                
print("Welcome to the CHN Program for Elemental Analysis")
print("")

# sample.formula = formula
i = 0
menu_choice = '1'
while True:
    if menu_choice == '1':
        sample_name = input("Sample Name: ")
        i =+ 1
        sample_no="sample_no_"+str(i)
        # instantiate Sample          
        sample_no = Sample(sample_name)    
        formula = sample_no.input_mf()
        mw = sample_no.calc_mw(formula)
        fw = mw
        em = sample_no.calc_exactmass(formula)
        elem_percents = sample_no.calc_percents(formula, mw)
        sample_no.print_chn_data(sample_name, formula, elem_percents, em, mw, fw)    
    
    sample_no.print_menu()
    #1)Enter new MF\n2)Enter experimental data\n3)Add Water of hydration\n4)Print to file\n5)Exit Program
    
    menu_choice = input("Input a menu option number : ")
    if menu_choice == '1':
        continue

    if menu_choice == '2':
        combustion_elements = ['C','H','N','S']
        # remove elements without measureable combusion products
        pct_experimental = {e: formula[e] for e in set(combustion_elements) & set(formula.keys())}
        print("Enter experimental determined percentages")
        for element in sorted(pct_experimental):
            while True:
                try:
                    elem_percents[element]['experimental'] = float(input("%"+element+": "))
                except ValueError:
                    print("Please input a decimal number.")
                    continue
                break
        try:
            hydrate
        except NameError:
            hydrate = 0
        sample_no.print_chn_data(sample_name, formula, elem_percents, em, mw, fw, hydrate)
        continue

    if menu_choice == '3':
        H2Ox = float(input("Molar Ratio : "))
        # new dictionary for hydrate formula
        hydrate_formula = sample_no.create_hydrate_formula(formula, H2Ox)
        hydrate = H2Ox
        fw = sample_no.calc_mw(hydrate_formula)
        elem_percents = sample_no.calc_percents(hydrate_formula, fw, elem_percents)
        sample_no.print_chn_data(sample_name, formula, elem_percents, em, mw, fw, hydrate)
        continue
    
    if menu_choice == '4':
        default = 'chn_results.txt'
        filename = input('Filename: [%s]' % default+' ')
        filename = filename or default
        try:
            hydrate
        except NameError:
            hydrate = 0
        sample_no.print_chn_data(sample_name, formula, elem_percents, em, mw, fw, hydrate, filename)
        continue
       
    if menu_choice == '5':
        print('Done.')
        break
