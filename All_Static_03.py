## Import libraries for analysis

import os
import re
import pandas as pd
from collections import Counter
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET

## Initialize lists for relevent data

list_of_file_names = []
list_of_nct_id = []
list_of_study_type = []
list_of_study_first_submitted_qc = []
list_of_phase = []
list_of_overall_status = []
list_of_lead_sponsor = []
set_of_all_data_items = set()
list_of_data_items = []
list_of_data_items_with_count = []

## Extract data field names and their counts and children from each xml file

#for path, dirs, files in os.walk("D:/AILON/data/ClinicalTrials/NCT0000xxxx/"):  
for path, dirs, files in os.walk("D:/AILON/data/ClinicalTrials/"):  
    for filename in files:
        #print(filename)
        list_of_file_names.append(filename)
        list_of_nct_id.append("")
        list_of_lead_sponsor.append("")
        list_of_study_type.append("")
        list_of_phase.append("")
        list_of_overall_status.append("")
        list_of_study_first_submitted_qc.append("")
        
        file = path + "/" + filename
        tree = ET.parse(file)
        root = tree.getroot()
        
        l = []
        #set_of_data_items = set()
        for element in root:
            l.append(element.tag)
            set_of_all_data_items.add(element.tag)
            
            if element.tag == 'study_type':
                list_of_study_type.pop(len(list_of_study_type)-1)
                list_of_study_type.append(element.text)
            
            if element.tag == 'phase':
                list_of_phase.pop(len(list_of_phase)-1)
                list_of_phase.append(element.text)
            
            if element.tag == 'overall_status':
                list_of_overall_status.pop(len(list_of_overall_status)-1)
                list_of_overall_status.append(element.text)
            
            if element.tag == 'study_first_submitted_qc':
                list_of_study_first_submitted_qc.pop(len(list_of_study_first_submitted_qc)-1)
                list_of_study_first_submitted_qc.append(element.text)
            
            if len(element.getchildren()) > 0:
                for subelement in element:
                    l.append(element.tag+'>'+subelement.tag)
                    set_of_all_data_items.add(element.tag+'>'+subelement.tag)
                    
                    if subelement.tag == 'nct_id':
                        list_of_nct_id.pop(len(list_of_nct_id)-1)
                        list_of_nct_id.append(subelement.text)
                    
                    if len(subelement.getchildren()) > 0:
                        for subsubelement in subelement:
                            l.append(element.tag+'>'+subelement.tag+'>'+subsubelement.tag)
                            set_of_all_data_items.add(element.tag+'>'+subelement.tag+'>'+subsubelement.tag)
                            
                            if (element.tag == 'sponsors' and subelement.tag == 'lead_sponsor' and subsubelement.tag == 'agency'):
                                list_of_lead_sponsor.pop(len(list_of_lead_sponsor)-1)
                                list_of_lead_sponsor.append(subsubelement.text)
                        
        #list_of_data_items.append(set_of_data_items)
        list_of_data_items_with_count.append(Counter(l))

df = pd.DataFrame(0,columns=set_of_all_data_items, index=list_of_nct_id)

i = 0
for file_data_items in list_of_data_items_with_count:
    for item in file_data_items:
        df.iloc[i][item] = list_of_data_items_with_count[i][item]
    i += 1

df['status_nm'] = list_of_overall_status
df['lead_sponsor_nm'] = list_of_lead_sponsor
df['study_type_nm'] = list_of_study_type
df['phase_nm'] = list_of_phase
df['study_first_submitted_qc_nm'] = list_of_study_first_submitted_qc

df.to_csv('occurrence_table.csv')