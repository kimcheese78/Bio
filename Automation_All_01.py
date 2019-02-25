import os
import pandas as pd
from biopandas.mol2 import PandasMol2
from biopandas.mol2 import split_multimol2
import numpy as np

default_DUDE_data_path = 'D:\\AILON\\data\\DUDE'
default_PDB_data_path = 'D:\\AILON\\data\\scPDB'
default_output_path = 'D:\\AILON\\data\\AILONET_DATA'

atlist = ['C.ani', 'C.ar','C.cat', 'C.1', 'C.2', 'C.3', 'N.am','N.pl3','N.ar','N.1','N.2','N.3','N.4', 'O.3',  'O.2','O.co2', 'S.2', 'S.3', 'S.o2', 'S.o','H','Mg','Ca','Fe','Zn','Si','I','P.3','Br','Cl','F']
atnum=[12,12,12,12,12,12,14,14,14,14,14,14,14,16,16,16,32,32,32,32,1,24,40,56,65,28,127,31,80,35,19]

# Get the list of Proteins in DUDE
os.chdir(default_DUDE_data_path)
proteins = os.listdir()

# Get the list of Proteins in PDB
os.chdir(default_PDB_data_path)


scPDB = os.listdir()

# Get the relation table between DUDE code and scPBD code
relation_table = pd.read_csv('D:\\AILON\\data\\PDB_DUDE_table_result.csv')

for protein in proteins:
    # Retrieve PDB code
    pdb_code = relation_table[relation_table['Target_name']==protein.upper()]['PDB_code'].to_string(index=False).lower().strip()
    # Target Directory Generation
    os.chdir(default_output_path)
    protein_directory = protein.upper() + '_' + pdb_code.upper()
    os.makedirs(protein_directory + '/' + protein_directory + '_Active')
    os.makedirs(protein_directory + '/' + protein_directory + '_Decoy')
    
    # Change directory to Data
    os.chdir(default_DUDE_data_path + '\\' + protein)
    
    ## Process Actives
    # Change sdf into pdbqt and Separate actives into their own files
    os.system("7z e actives_final.sdf.gz")
    os.system('obabel actives_final.sdf -O actives_filt.pdbqt --unique title')
    os.system('obabel actives_filt.pdbqt -O ' + default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Active\\actives_mole.pdbqt -m')
    
    # Change the file name with CHEMBL id
    os.chdir(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Active')
    files = os.listdir()
    for file in files:
        file_name = next(open(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Active\\' + file)).split(' ')[4].rstrip()
        os.rename(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Active\\' + file, default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Active\\' + file_name + '.pdbqt')
    
    # Change directory to Data
    os.chdir(default_DUDE_data_path + '\\' + protein)
    
    ## Process Decoys
    # Change sdf into pdbqt and Separate decoys into their own files
    os.system("7z e decoys_final.sdf.gz")
    os.system('obabel decoys_final.sdf -O decoys_filt.pdbqt --unique title')
    os.system('obabel decoys_filt.pdbqt -O ' + default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Decoy\\decoys_mole.pdbqt -m')
    # Change the file name with ZINC id
    os.chdir(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Decoy')
    files = os.listdir()
    for file in files:
        file_name = next(open(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Decoy\\' + file)).split(' ')[4].rstrip()
        os.rename(default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Decoy\\' + file, default_output_path + '\\' + protein_directory + '\\' + protein_directory + '_Decoy\\' + file_name + '.pdbqt')
    
    
    ## Process Config file
    df_count = relation_table[relation_table['Target_name']==protein.upper()]['count']
    df_count.reset_index()
    count = df_count.iloc[0]
    if count > 0:
        # Search folder in scPDB for protein
        pdb_code = relation_table[relation_table['Target_name']==protein.upper()]['PDB_code'].to_string(index=False).lower().strip()
        scPDB_folder_names = [name for name in scPDB if pdb_code in name]
        os.chdir(default_PDB_data_path + '\\' + scPDB_folder_names[0])
        
        pmol = PandasMol2().read_mol2('site.mol2')
        x = pmol.df['x']
        y = pmol.df['y']
        z = pmol.df['z']
        
        weight=[]
        
        for atomtypelist in range(len(pmol.df['atom_type'])):
            weight.append(atnum[atlist.index(pmol.df['atom_type'][atomtypelist])])
        
        cmx = sum(np.multiply(x, weight)) / sum(weight)
        cmy = sum(np.multiply(y, weight)) / sum(weight)
        cmz = sum(np.multiply(z, weight)) / sum(weight)
        
        file_name = pdb_code + '_config.txt'
        f = open(default_output_path + '\\' + protein_directory + '\\' + file_name,'w+')
        
        f.write('receptor = ' + pdb_code + '.pdbqt\n\n')
        f.write('center_x = %18.15f\n' % cmx)
        f.write('center_y = %18.15f\n' % cmy)
        f.write('center_z = %18.15f\n\n' % cmz)
        f.write('size_x = 20\n')
        f.write('size_y = 20\n')
        f.write('size_z = 20\n\n')
        f.write('exhaustiveness = 8\n\n')
        f.write('log = ' + pdb_code + '_log.txt\n\n')
        f.write('num_mode = 1')
        
        f.close()