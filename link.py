from pathlib import Path
import json
import sys
from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit.Chem import AllChem as rdkit
from collections import defaultdict
import py3Dmol
from IPython.display import Image
import matplotlib.pyplot as plt
import stk
import stko
import numpy as np
import shutil
import os
import string
from collections import OrderedDict
from ase import Atoms
from ase.io import read
from collections import defaultdict
import networkx as nx
import pulp
import logging
from ase.units import Hartree, mol, kJ, Bohr

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QTextEdit
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWebEngineCore import QWebEngineSettings 
from PyQt6.QtCore import QUrl, QTimer
from PyQt6.QtGui import QImage, QPixmap, QIcon, QFont, QTextCursor
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QLabel
import io
from rdkit.Chem import Draw

class MoleculeViewer(QWidget):
    def __init__(self, rdkit_mol):
        super().__init__()
        self.rdkit_mol = rdkit_mol
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Molecule Visualization')
        self.setGeometry(150, 150, 600, 600)

        layout = QVBoxLayout()
        img = Draw.MolToImage(self.rdkit_mol, size=(500, 500))
        
        # Convert PIL Image to QPixmap
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        qimg = QImage.fromData(buffer.getvalue())
        pixmap = QPixmap.fromImage(qimg)

        # Create QLabel and set the pixmap
        label = QLabel()
        label.setPixmap(pixmap)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        layout.addWidget(label)
        self.setLayout(layout)

def check_type_and_shape(value, name, expected_type, expected_shape):
    value = np.asarray(value)
    if value.size == 0:
        sys.exit(f'ERROR: No data found in the QM Hessian output file(s) for "{name}".'
                    ' Exiting...\n\n')
    elif value.dtype != np.dtype(expected_type):
        raise TypeError(f'"{name}" property expected a type of "{np.dtype(expected_type)}",'
                        f' but got "{value.dtype}" for the QM Hessian output.')
    elif value.shape != expected_shape:
        sys.exit(f'ERROR: "{name}" property expected a shape of "{expected_shape}", but got '
                    f'"{value.shape}" for the QM Hessian output. Exiting...\n\n')
    return value

def show_stk_mol(stk_mol):
    data = rdkit.MolToMolBlock(stk_mol.to_rdkit_mol())
    p = py3Dmol.view(
        data=data,
        style={'stick':{'colorscheme':'cyanCarbon'}}, 
        width=400,
        height=400,
    )
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    p.show()



class link():
    def __init__(self, molecule_name, BB_name, BB_connector):
        self.molecule_name = molecule_name #分子名稱
        self.BB_names = BB_name #BB名稱
        self.BB_connector = BB_connector #BB的接頭       
        self.deleter_num = []
        self.dihedral_number_list = [] #dihedral編號                
        self.dihedral_number_pdb = "" #dihedral編號存成文字檔用來放到.pdb檔內         
        self.molecule_smile = {}  #分子的smiles
        self.new_BBname = []  #BB的新名稱       
        self.name_counter = {} #BB名字有沒有重複
        self.atom_name = []  # 讀取rtp檔裡的各個原子名稱照rtp順序存成list
        self.dihedral_rtp = []   # dihedral以原子名稱編號
        self.dict_atomid = {}  #BB編號對應分子的原子編號 
        self.atom_type = []
        self.pdb_information = []
        self.count = [0]
        self.number = 0
        self.atom_name_dic = {}
        self.BB_dic = {}
        self.qforce_atomtype = []
        self.dihedral_dic = {} #用atomtype找dihedral
   
    
        self.create_FF_folders() #建立FF所需的資料夾
        self.STK_create_building_blocks() #讀取BB使用stk定義為一個一個building block
        self.build_molecule() #將各個BB使用stk連在一起
        self.show_molecule() #顯示分子
        self.renameBB_copy() #複製BB的rtp檔並重新命名
        self.replace_BBname_in_rtp() #替換rtp檔裡的BB名稱
        self.pdb_file_information() #pdb檔的資訊
        self.rtp_add_bonds() #加入rtp資訊
        #self.update_dihedral_number_pdb_create_SCAN() #找出對稱的dihedral並建立dihedral.txt以加入在.pdb檔上(需計算完Guassian後再執行一次)
        #self.establish_atp() #建立atp檔
        #self.map_atomnames_to_dihedral() #建立dihedral的rtp檔
        #self.rtp_add_impropers() #加入rtp資訊
        #self.rtp_add_dihedrals() #加入rtp資訊
        #self.arrange_rtp() #將rtp檔排整齊
        #self.forcefield_itp() #建立forcefield.itp，新增相連接的項目 (需辨認出類似的項還沒完成)
        #self.forcefield_itp_atomtype() #填入itp atomtype項  
        #self.DDEC6_charge() #計算DDEC6 charge (需計算完電荷後再執行一次)
        #self.improper_number() #improper的編號
        


     
           
    def DDEC6_charge(self):
       
        n_atoms = self.molecule_stk.get_num_atoms()
        with open(f'{self.molecule_folder_path}/QM_GASPHASE/{self.molecule_name}.log', 'r') as file:
            for index, line in enumerate(file):
                if "Input orientation" in line:
                    start = index + 5
                    end = index + n_atoms + 5
                    break

            elements = []
            coords = []
            file.seek(0)  # 將檔案指標重設為開始位置
            for i, line in enumerate(file):
                if start <= i < end:
                    element = line.split()[1]
                    elements.append(int(element))
                    coord = line.split()[3:6]
                    coords.append([float(i) for i in coord])


        with open(f'{self.molecule_folder_path}/QM_GASPHASE/xtb/wbo','r') as file:

            b_orders = [[0] * n_atoms for _ in range(n_atoms)]  # 初始化为0
            file.seek(0)
            lines = file.readlines()
            for line in lines:
                if not line.split():
                    break
                atom_1 = int(line.split()[0])-1
                atom_2 = int(line.split()[1])-1
                bond_order = float(line.split()[2])

                
                b_orders[atom_1][atom_2] = bond_order

            elements = check_type_and_shape(elements, 'elements', int, (n_atoms,))
            coords = check_type_and_shape(coords, 'coords', float, (n_atoms, 3))
            b_orders = check_type_and_shape(b_orders, 'b_orders', float, (n_atoms, n_atoms))
            n_bonds = b_orders.sum(axis=1).round().astype(int)

            elements = np.array(elements)
            coords = np.asfarray(coords, float)
            coords = np.reshape(coords, (-1, 3))
            coords = np.array(coords)* Bohr


        graph = nx.Graph()
        for i_idx, i_elem in enumerate(elements):
            graph.add_node(i_idx,  elem=i_elem ,n_bonds=n_bonds[i_idx],
                                coords=coords[i_idx])
            # add bonds
            for j_idx, j_elem in enumerate(elements):
                b_order = b_orders[i_idx, j_idx]
                if b_order > 0.3:
                    id1, id2 = sorted([i_elem, j_elem])
                    b_order_half_rounded = np.round(b_order*2)/2
                    vec = coords[i_idx] - coords[j_idx]
                    dist = np.sqrt((vec**2).sum())
                    graph.add_edge(i_idx, j_idx, vector=vec, length=dist, order=b_order,
                                        type=f'{id1}({b_order_half_rounded}){id2}', n_rings=0)


        rings = nx.minimum_cycle_basis(graph) #循環
        rings3 = [r for r in rings if len(r) == 3]

        for i in range(n_atoms):
            graph.nodes[i]['n_ring'] = sum([i in ring for ring in rings]) #儲存原子的環數
        #
        for atoms in graph.edges: #全部的鍵結
            ring_members = [set(atoms).issubset(set(ring)) for ring in rings]
            graph.edges[atoms[0],atoms[1]]['n_rings'] = sum(ring_members) #鍵結再幾個環裡面
            graph.edges[atoms[0],atoms[1]]['in_ring'] = any(ring_members)
            graph.edges[atoms[0],atoms[1]]['in_ring3'] = any(set(atoms).issubset(set(ring))for ring in rings3)

        atom_ids = [[] for _ in range(n_atoms)]

        bonds, angles, dihedrals = [], [], []
        
        #辨認topology
        for i in range(n_atoms):
            neighbors = nx.bfs_tree(graph, source=i,  depth_limit=3).nodes
            for n in neighbors:
                paths = nx.all_simple_paths(graph, i, n, cutoff=3)
                for path in paths:

                    if len(path) == 2 and path[0] < path[1]: #後面的條件是保證不會重複辨認同樣的情況
                        bonds.append(path)
                    elif len(path) == 3 and path[0] < path[2]:
                        angles.append(path)
                    elif len(path) == 4 and path[1] < path[2]:
                        dihedrals.append(path)

        dihedrals.sort(key=lambda x: [x[1], x[2], x[0], x[3]]) #dihedral的第二項最小為優先順序，再來是第三項，以此類推

        #self.bonds = bonds
        #self.angles = angles
        #self.dihedrals = dihedrals
        #print(dihedrals)
        
        for i in range(n_atoms):
            neighbors = nx.bfs_tree(graph, i, depth_limit=4).nodes
            for n in neighbors:
                paths = nx.all_simple_paths(graph, i, n, cutoff=4)
                for path in map(nx.utils.pairwise, paths):
                    types = [graph.edges[edge[0],edge[1]]['type'] for edge in path]
                    atom_ids[i].append("-".join(types))
                atom_ids[i].sort()
        if 4 < 0:
            atom_ids = [i for i in range(n_atoms)]

        list = []

        for n in range(n_atoms):
            if n in [item for sub in list for item in sub]:
                continue
            eq = [i for i, a_id in enumerate(atom_ids) if atom_ids[n] == a_id]
            list.append(eq)  

        atomtype_dic = {}
        for n,atoms in enumerate(list):
            for atom in atoms:
                atomtype_dic[str(atom+1)] = n
                
        dihedral = []
        dihedral_atomtype = []
        dihedral_BB = []
        for d in self.dihedral_number_list:
            a1 = []
            a4 = []
            for i in self.molecule_stk.get_bonds():
                link1 = d[1]-1
                link2 = d[2]-1
                a2 = i.get_atom1().get_id()+1
                a3 = i.get_atom2().get_id()+1
                cond1 = (i.get_atom1().get_id() != link1 or i.get_atom2().get_id() != link2)
                cond2 = (i.get_atom1().get_id() != link2 or i.get_atom2().get_id() != link1)
                cond3 = i.get_atom1().get_id() == link1 or i.get_atom2().get_id() == link1
                cond4 = i.get_atom1().get_id() == link2 or i.get_atom2().get_id() == link2
                if cond1 and cond2:
                    if cond3 or cond4:
                        if i.get_atom1().get_id()+1 == link1+1:
                            a1.append(f'{i.get_atom2().get_id()+1}')
                        if i.get_atom2().get_id()+1 == link1+1:
                            a1.append(f'{i.get_atom1().get_id()+1}')
                        if i.get_atom1().get_id()+1 == link2+1:
                            a4.append(f'{i.get_atom2().get_id()+1}')
                        if i.get_atom2().get_id()+1 == link2+1:
                            a4.append(f'{i.get_atom1().get_id()+1}')
                            
            for atom1 in a1:
                for atom4 in a4:
                    #print(f'{atom1} {d[1]} {d[2]} {atom4} ')
                    #print(f'atom type : {atomtype_dic[atom1]} {atomtype_dic[str(d[1])]} {atomtype_dic[str(d[2])]} {atomtype_dic[atom4]}')
                    #print(f'{self.atom_name_dic[str(atom1)]}  {self.atom_name_dic[str(d[1])]}  {self.atom_name_dic[str(d[2])]}  {self.atom_name_dic[str(atom4)]}')
                    #print(f'{self.BB_dic[str(atom1)]}  {self.BB_dic[str(d[1])]}  {self.BB_dic[str(d[2])]}  {self.BB_dic[str(atom4)]}')
                    dihedral.append(f'{atom1} {d[1]} {d[2]} {atom4}')
                    dihedral_atomtype.append(f'{atomtype_dic[atom1]} {atomtype_dic[str(d[1])]} {atomtype_dic[str(d[2])]} {atomtype_dic[atom4]}')
                    dihedral_BB.append(f'{self.BB_dic[str(atom1)]}  {self.BB_dic[str(d[1])]}  {self.BB_dic[str(d[2])]}  {self.BB_dic[str(atom4)]}')

        dihedral_dic = {}
        for a,b,c in zip(dihedral,dihedral_atomtype,dihedral_BB):
            dihedral_dic[a] = {'atomtype' : [b] , 'BB' : [c]}



        for atomnumber, value in dihedral_dic.items():
            atomtype = value['atomtype']

            atomtype_verse = [' '.join(item.split()[::-1]) for item in atomtype]
            
            BB = value.get('BB') 
            for atomnumber1,value1 in dihedral_dic.items():
                atomtype1 = value1['atomtype']
                BB1 = value1.get('BB')
                if (atomtype1 == atomtype_verse or atomtype1 == atomtype) and BB1 not in [v['BB'] for v in self.dihedral_dic.values()]:
                    self.dihedral_dic[atomnumber1] = {'atomtype': atomtype1, 'BB': BB1}   

        #print(f'dihedral_dic:{self.dihedral_dic}')


        with open (f'{self.molecule_folder_path}/QM_GASPHASE/DDEC6_even_tempered_net_atomic_charges.xyz','r')as f:
            lines = f.readlines()
            q = []
            for line in lines[2:]:
                if not line.split():
                    break
                q.append(line.split()[-1])

        q = np.array(q, dtype=float)
        q = np.reshape(q,(1,n_atoms))   

        avg_terms = []
        for term in q:
            term = np.array(term)
            for tlist in list:

                term[tlist] = term[tlist].mean().round(5)
            avg_terms.append(term)

        q_avg = avg_terms[0]        

        total = q_avg.sum()
        q_integer = int(round(total))
        extra = int(round(100000 * round(total - q_integer, 5)))

        if extra != 0:
            if extra > 0:
                sign = 1
            else:
                sign = -1
                extra = - extra

            n_eq = [len(tlist) for tlist in list]
            eq_no = [f"{i:05d}" for i, _ in enumerate(n_eq)]

            var = pulp.LpVariable.dicts("x", eq_no, lowBound=0, cat='Integer')
            prob = pulp.LpProblem('prob', pulp.LpMinimize)
            prob += pulp.lpSum([var[n] for n in eq_no])
            prob += pulp.lpSum([eq * var[eq_no[i]] for i, eq in enumerate(n_eq)]) == extra
            prob.solve()
        
            if prob.status == 1:
                for i, v in enumerate(prob.variables()):
                    q_avg[list[i]] -= sign * v.varValue / 100000
            else:
                print('Failed to equate total of charges to the total charge of '
                        'the system. Do so manually')

        number = 0
        for n,rtp in enumerate(self.new_BBname):
            with open(f"{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/{rtp}.rtp",'r+') as f:
                lines = f.readlines()
                for index,line in enumerate(lines):
                    if ' atoms ' in line:
                        start = index + 1
                    elif ' bonds ' in line:
                        end = index
        # 将文件指针移动到文件末尾
                f.seek(0)  # 将文件指针移动到文件开头
                for i in lines[:start]:
                    f.write(i)  # 将原始内容写回文件
                
                for i in lines[start:end]:
                    if i.strip() != '':
                        formatted_line = "     {:<16} {:<14.5f} {:<6} \n".format(i.strip(), q_avg[number], str(n))
                        f.write(formatted_line)
                        number += 1
                f.write('\n')
                for i in lines[end:]:
                    f.write(i)  # 将剩余内容写回文件

    def forcefield_itp_atomtype(self):
        if not os.path.exists(f'../molecule/{self.molecule_name}/QM_GASPHASE/DDEC_atomic_Rcubed_moments.xyz'):
            print(f"無法填入charge和LJ參數，沒有DDEC6的結果，請先計算DDEC6 charge")
            raise SystemExit
        else:
            molecule = read(f'../molecule/{self.molecule_name}/QM_GASPHASE/DDEC_atomic_Rcubed_moments.xyz')

        B2nm = 0.0529177
        H2kJ = 2625.5
        Vf = {'H':7.842800,'C':34.939018,'N':25.922105,'O':22.1}
        Bf = {'H':6.5,'C':46.6,'N':24.2,'O':15.6}
        Rf = {'H':3.099151,'C':3.93063,'N':3.250329,'O':3.02}

        f = open(f'../molecule/{self.molecule_name}/QM_GASPHASE/DDEC_atomic_Rcubed_moments.xyz','r')
        lines = f.readlines()

        natom = int(lines[0].split()[0])
        SYMB = []
        Vi = []
        for i in range(2,natom+2):
            SYMB.append(lines[i].split()[0])
            Vi.append(float(lines[i].split()[-1]))
        Bi = []
        for i, V in enumerate(Vi):
            Bi.append((V/Vf[SYMB[i]])**2.0*Bf[SYMB[i]])
        Ai = []
        for i, B in enumerate(Bi):
            Ri = (Vi[i]/Vf[SYMB[i]])**(1.0/3.0)*Rf[SYMB[i]]
            Ai.append(0.5*B*(2*Ri)**6.0)

        sigma = []
        eps = []
        for i in range(len(Bi)):
            sigma.append((Ai[i]/Bi[i])**(1.0/6.0)*B2nm)
            eps.append(Bi[i]**2.0/(4.0*Ai[i])*H2kJ)
            #print (i+1,SYMB[i],(Ai[i]/Bi[i])**(1.0/6.0)*B2nm,Bi[i]**2.0/(4.0*Ai[i])*H2kJ)


        rtp_path = []
        with open(f"../molecule/{self.molecule_name}/STRUCTURE/pdb.txt","r") as f :
            lines = f.readlines()
            for line in lines:
                tmp = line.split()
                if (f"../molecule/{self.molecule_name}/FORCEFIELD/{self.molecule_name}.ff/{tmp[1]}.rtp") not in rtp_path:
                    rtp_path.append(f"../molecule/{self.molecule_name}/FORCEFIELD/{self.molecule_name}.ff/{tmp[1]}.rtp")
                                        
        from collections import defaultdict

        # 使用defaultdict来初始化字典，避免了检查键是否存在的逻辑
        sigma_dict = defaultdict(list)
        eps_dict = defaultdict(list)

        # 合并循环，直接填充字典
        for n, name in enumerate(self.atom_type):
            sigma_dict[name].append(sigma[n])
            eps_dict[name].append(eps[n])

        # 计算平均值
        sum_sigma_dict = {key: sum(value)/len(value) for key, value in sigma_dict.items()}
        sum_eps_dict = {key: sum(value)/len(value) for key, value in eps_dict.items()}

        with open(f"{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/forcefield.itp", 'r+') as file:
            content = file.read()
            file.seek(0, 0)  # 将文件指针移动到文件开头
            file.write("[ defaults ]\n")
            file.write("; nbfunc      comb-rule       gen-pairs        fudgeLJ  fudgeQQ\n")
            file.write("    1             3             yes              0.5      0.5\n\n")
            file.write("[ atomtypes ]\n")
            file.write(";     name      mass      charge       ptype     sigma            eps\n")
            for a in sum_sigma_dict:
                if '_H' in a:
                    file.write(f'   {a:>8}  {1.0080:>8}  {0.000:>8}   {"A":>8}  {sum_sigma_dict[a]:>14.5e}  {sum_eps_dict[a]:>14.5e}\n')
                elif '_C' in a:
                    file.write(f'   {a:>8}  {12.011:>8}  {0.000:>8}   {"A":>8}  {sum_sigma_dict[a]:>14.5e}  {sum_eps_dict[a]:>14.5e}\n')
                elif '_N' in a:
                    file.write(f'   {a:>8}  {14.007:>8}  {0.000:>8}   {"A":>8}  {sum_sigma_dict[a]:>14.5e}  {sum_eps_dict[a]:>14.5e}\n')
                elif '_O' in a:
                    file.write(f'   {a:>8}  {15.999:>8}  {0.000:>8}   {"A":>8}  {sum_sigma_dict[a]:>14.5e}  {sum_eps_dict[a]:>14.5e}\n')
                elif '_S' in a:
                    file.write(f'   {a:>8}  {32.060:>8}  {0.000:>8}   {"A":>8}  {sum_sigma_dict[a]:>14.5e}  {sum_eps_dict[a]:>14.5e}\n')
            file.write(f'\n{content}')  # 将原始内容写回文件            
            

    def forcefield_itp(self):
        if not os.path.exists(f'../molecule/{self.molecule_name}/QM_GASPHASE/{self.molecule_name}.xyz'):
            print(f"無法建立forcefield.itp，沒有optimize後的{self.molecule_name}.xyz檔案在QM_GASPHASE資料夾中，請先optimize分子")
            raise SystemExit
        else:
            molecule = read(f'../molecule/{self.molecule_name}/QM_GASPHASE/{self.molecule_name}.xyz')


        bond = []
        bond_number = []
        for n in range(len(self.dihedral_number_list)):
            for i in self.molecule_stk.get_bonds():
                link1 = self.dihedral_number_list[n][1]-1
                link2 = self.dihedral_number_list[n][2]-1
                cond1 = (i.get_atom1().get_id() == link1) and (i.get_atom2().get_id() == link2)
                cond2 = (i.get_atom1().get_id() == link2) and (i.get_atom2().get_id() == link1)
                if cond1 or cond2:
                    if i.get_atom1().get_id() == link1:
                        bond_number.append((f'{link1+1}',f'{link2+1}'))
                        bond.append((self.atomtype_dic[f'{link1+1}'],self.atomtype_dic[f'{link2+1}']))
                    if i.get_atom2().get_id() == link1:
                        bond_number.append((f'{link2+1}',f'{link1+1}'))
                        bond.append((self.atomtype_dic[f'{link2+1}'],self.atomtype_dic[f'{link1+1}']))

        bond_set = []
        for i in bond:
            if (i[0],i[1]) not in bond_set and (i[1],i[0]) not in bond_set:
                bond_set.append((i[0],i[1]))

        Angle = []
        Angle_number = []
        for n in range(len(self.dihedral_number_list)):
            for i in self.molecule_stk.get_bonds():
                link1 = self.dihedral_number_list[n][1]-1
                link2 = self.dihedral_number_list[n][2]-1
                cond1 = (i.get_atom1().get_id() != link1 or i.get_atom2().get_id() != link2)
                cond2 = (i.get_atom1().get_id() != link2 or i.get_atom2().get_id() != link1)
                cond3 = i.get_atom1().get_id() == link1 or i.get_atom2().get_id() == link1
                cond4 = i.get_atom1().get_id() == link2 or i.get_atom2().get_id() == link2
                if cond1 and cond2:
                    if cond3 or cond4:
                        if i.get_atom1().get_id()+1 == link1+1:
                            Angle_number.append((f'{link2+1}',f'{link1+1}',f'{i.get_atom2().get_id()+1}'))
                            Angle.append((self.atomtype_dic[f'{link2+1}'],self.atomtype_dic[f'{link1+1}'], self.atomtype_dic[f'{i.get_atom2().get_id() + 1}']))
                        if i.get_atom2().get_id()+1 == link1+1:
                            Angle_number.append((f'{link2+1}',f'{link1+1}',f'{i.get_atom1().get_id()+1}'))
                            Angle.append((self.atomtype_dic[f'{link2+1}'],self.atomtype_dic[f'{link1+1}'], self.atomtype_dic[f'{i.get_atom1().get_id() + 1}']))
                        if i.get_atom1().get_id()+1 == link2+1:
                            Angle_number.append((f'{link1+1}',f'{link2+1}',f'{i.get_atom2().get_id()+1}'))
                            Angle.append((self.atomtype_dic[f'{link1+1}'],self.atomtype_dic[f'{link2+1}'], self.atomtype_dic[f'{i.get_atom2().get_id() + 1}']))
                        if i.get_atom2().get_id()+1 == link2+1:
                            Angle_number.append((f'{link1+1}',f'{link2+1}',f'{i.get_atom1().get_id()+1}'))
                            Angle.append((self.atomtype_dic[f'{link1+1}'],self.atomtype_dic[f'{link2+1}'], self.atomtype_dic[f'{i.get_atom1().get_id() + 1}']))

        Angle_set = []
        for i in Angle:
            if (i[0],i[1],i[2]) not in Angle_set and (i[2],i[1],i[0]) not in Angle_set:
                Angle_set.append((i[0],i[1],i[2]))
        
        proper = []
        proper_number = []
        for n in range(len(self.dihedral_number_list)):
            Angle_proper = []
            for i in self.molecule_stk.get_bonds():
                link1 = self.dihedral_number_list[n][1]-1
                link2 = self.dihedral_number_list[n][2]-1
                cond1 = (i.get_atom1().get_id() != link1 or i.get_atom2().get_id() != link2)
                cond2 = (i.get_atom1().get_id() != link2 or i.get_atom2().get_id() != link1)
                cond3 = i.get_atom1().get_id() == link1 or i.get_atom2().get_id() == link1
                cond4 = i.get_atom1().get_id() == link2 or i.get_atom2().get_id() == link2

                if cond1 and cond2:
                    if cond3 or cond4:
                        if i.get_atom1().get_id()+1 == link1+1:
                            Angle_proper.append((f'{link2+1}',f'{link1+1}', f'{i.get_atom2().get_id() + 1}'))
                        if i.get_atom2().get_id()+1 == link1+1:
                            Angle_proper.append((f'{link2+1}',f'{link1+1}', f'{i.get_atom1().get_id() + 1}'))
                        if i.get_atom1().get_id()+1 == link2+1:
                            Angle_proper.append((f'{link1+1}',f'{link2+1}', f'{i.get_atom2().get_id() + 1}'))
                        if i.get_atom2().get_id()+1 == link2+1:
                            Angle_proper.append((f'{link1+1}',f'{link2+1}', f'{i.get_atom1().get_id() + 1}'))
                            
            for n,atom in enumerate(Angle_proper):
                for i in self.molecule_stk.get_bonds():
                    link0 = atom[0]
                    link1 = atom[1]
                    link2 = atom[2]
                    cond1 = i.get_atom1().get_id() != int(link1)-1 and i.get_atom2().get_id() != int(link1)-1
                    cond2 = i.get_atom1().get_id() == int(link2)-1 or i.get_atom2().get_id() == int(link2)-1
                    if cond1:
                        if cond2:
                            if i.get_atom1().get_id()+1 == int(link2):
                                proper_number.append((f'{link0}',f'{link1}',f'{link2}', f'{i.get_atom2().get_id() + 1}'))
                                proper.append((self.atomtype_dic[f'{link0}'],self.atomtype_dic[f'{link1}'],self.atomtype_dic[f'{link2}'], self.atomtype_dic[f'{i.get_atom2().get_id() + 1}']))
                            if i.get_atom2().get_id()+1 == int(link2):
                                proper_number.append((f'{link0}',f'{link1}',f'{link2}', f'{i.get_atom1().get_id() + 1}'))
                                proper.append((self.atomtype_dic[f'{link0}'],self.atomtype_dic[f'{link1}'],self.atomtype_dic[f'{link2}'], self.atomtype_dic[f'{i.get_atom1().get_id() + 1}']))

                
        proper_set = []
        for i in proper:
            if (i[0],i[1],i[2],i[3]) not in proper_set and (i[3],i[2],i[1],i[0]) not in proper_set:
                proper_set.append((i[0],i[1],i[2],i[3]))

        '''
        dihedral_atomtype = [] #dihedral的atomtype
        dihedral_atomtype_set = [] #dihedral的atomtype的集合(重複幾種atomtype的組合)
        dihedral_atomtype_dic = {} #dihedral的atomtype的字典
        for a in self.dihedral_dic:
            a1,a2,a3,a4 = a.split()[0] , a.split()[1], a.split()[2], a.split()[3]
            
            atom1, atom2, atom3, atom4 =  self.atomtype_dic[f'{a1}'] , self.atomtype_dic[f'{a2}'], self.atomtype_dic[f'{a3}'], self.atomtype_dic[f'{a4}']
            atom_type_tuple = (atom1,atom2,atom3,atom4)
            dihedral_atomtype.append(atom_type_tuple)
            


            if atom_type_tuple in dihedral_atomtype_dic:
                if not isinstance(dihedral_atomtype_dic[atom_type_tuple], list):
                    dihedral_atomtype_dic[atom_type_tuple] = [dihedral_atomtype_dic[atom_type_tuple]]
                dihedral_atomtype_dic[atom_type_tuple].append((atom1,atom2,atom3,atom4))
            else:
                dihedral_atomtype_dic[atom_type_tuple] = (atom1,atom2,atom3,atom4)
       

        for i in dihedral_atomtype:
            if (atom1,atom2,atom3,atom4) not in dihedral_atomtype_set and (atom4,atom3,atom2,atom1)not in dihedral_atomtype_set:
                dihedral_atomtype_set.append((atom1,atom2,atom3,atom4))
        
        print(dihedral_atomtype_set)
        print(dihedral_atomtype_dic)
        '''

        # 读取XYZ文件
        bond_distance = []
        bond_info = {}
        for a,n in zip(bond,bond_number):
            distance = molecule.get_distance(int(n[0])-1, int(n[1])-1)
            bond_distance.append(distance)


            key = tuple(sorted(a))  # 使用排序后的元组作为键，包含键和原子编号

            if key in bond_info:
                bond_info[key].append((distance, (n[0], n[1])))  # 将距离和原子编号作为元组添加到列表中
            else:
                bond_info[key] = [(distance, (n[0], n[1]))] 

            '''
            key = tuple(sorted(a))  # 排序元组 (a[0], a[1])，确保键的顺序一致
            if key in bond_info:
                bond_info[key].append(distance)
            else:
                bond_info[key] = [distance]
            '''
            '''   
            if a in bond_info: #順序是正的
                bond_info[a].append(distance) #順序反改為正
            if (a[1],a[0]) in bond_info:
                bond_info[(a[1],a[0])].append(distance)
            else:
                bond_info[a] = [distance]
            ''' 



        bond_itp = []
        for key,value in bond_info.items():
            distances = [item[0] for item in value]
            atom_numbers = [item[1] for item in value]  
            average_distance = sum(distances) / len(distances) 

            bond_itp.append((key, "{:.5f}".format(average_distance * 0.1), atom_numbers))


        Angle_ase = []
        Angle_info = {}

        for a,n in zip(Angle,Angle_number):
            angle_ase = molecule.get_angle(int(n[0])-1, int(n[1])-1, int(n[2])-1)
            Angle_ase.append(angle_ase)
            
            angle_key = frozenset(a)

            if a in Angle_info:
                Angle_info[a].append((angle_ase, (n)))
            else:
                Angle_info[a] = [(angle_ase, (n))]


        Angle_itp = []
        for key,value in Angle_info.items():
            angle = [item[0] for item in value]
            atom_numbers = [item[1] for item in value]
            average_angle = sum(angle) / len(angle)

            Angle_itp.append((key , "{:.5f}".format(average_angle ),atom_numbers))


        Proper_ase = []
        Proper_info = {}
        for a,n in zip(proper,proper_number):
            proper_ase = molecule.get_dihedral(int(n[0])-1, int(n[1])-1, int(n[2])-1,int(n[3])-1)
            if 160 < proper_ase < 200 :
                proper_ase = 180
            if 0 < proper_ase < 20 or 330 < proper_ase <380:
                proper_ase = 0

            Proper_ase.append(proper_ase)
            proper_key = frozenset(a)

            if a in Proper_info:
                Proper_info[a].append((proper_ase, (n)))
            else:
                Proper_info[a] = [(proper_ase, (n))]

        proper_itp = []
        for key,value in Proper_info.items():
            angle = [item[0] for item in value]
            atom_numbers = [item[1] for item in value]
            average_angle = sum(angle) / len(angle)
            proper_itp.append((key , "{:.3f}".format(average_angle),atom_numbers))

        for name in set(self.BB_names):
            source_file = f"qforce_itp/{name}.itp"
            #source_file = f"C://Users//egg45//OneDrive - NTHU//resreach//Data base//qforce_RTP//improper\{name}.rtp"
            destination_folder = f"{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff"
            destination_file = f"{destination_folder}/{name}.itp"
            os.makedirs(destination_folder, exist_ok=True)    
            print(f"Copying {source_file} to {destination_file}")        
            shutil.copy(source_file, destination_file)



        with open(f'{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/forcefield.itp', 'w') as file:
            file.write('[ bondtypes ]\n') 
            BB_values = {}
            for bondtype in bond_itp:
                BB_values['BB1'] = bondtype[0][0].split('_')[0]
                BB_values['BB2'] = bondtype[0][1].split('_')[0]
                atom_numbers = ','.join(f'({a[0]},{a[1]})' for a in bondtype[2])

                for key in BB_values:
                    if BB_values[key] in ('1BA', '1BB', '1BC', '1BD', '1BE', '1BF', '1BG', '1BH','TLA','TLB','TLC','BUA'):
                        BB_values[key] = '1BA'

                BB1 = BB_values['BB1']
                BB2 = BB_values['BB2']
                with open(f'../force_constant_database/bond.txt') as f:
                    lines = f.readlines()
                    found_in_database = False
                    for line in lines:
                        tmp = line.split()
                        try:
                            #有再database中找到
                            cond1 = (BB1 in tmp[0] and BB2 in tmp[1]) 
                            cond2 = (BB2 in tmp[0] and BB1 in tmp[1])
                            if cond1 or cond2:
                                file.write(f'     {bondtype[0][0]:<9} {bondtype[0][1]:<38} {1:<5} {bondtype[1]:<18} {tmp[2]}     ;{atom_numbers}\n')
                                found_in_database = True
                                break
                        except IndexError:
                            continue

                if not found_in_database:
                    for line in lines:
                        tmp = line.split()
                        try:
                            # 查找自己BB的CC、CN force constant
                            cond3 = (bondtype[0][0]in tmp[0] and bondtype[0][1][4] in tmp[1][4]) or (bondtype[0][0]in tmp[1] and bondtype[0][1][4] in tmp[0][4])
                            cond4 = (bondtype[0][1] in tmp[0] and bondtype[0][0][4] in tmp[1][4]) or (bondtype[0][1] in tmp[1] and bondtype[0][0][4] in tmp[0][4])
                            if cond3 or cond4:
                                file.write(f'     {bondtype[0][0]:<9} {bondtype[0][1]:<38} {1:<5} {bondtype[1]:<18} {tmp[2]}     ;{atom_numbers} ;BBdata CC、CN\n')
                                break
                        except IndexError:
                            continue
               
            file.write('\n[ angletypes ]\n')
            for angletype in Angle_itp:
                BB1 = angletype[0][0].split('_')[0]
                BB2 = angletype[0][1].split('_')[0]
                BB3 = angletype[0][2].split('_')[0]
                a1 = angletype[0][0].split('_')[1]
                a2 = angletype[0][1].split('_')[1]
                a3 = angletype[0][2].split('_')[1]
                BB_values = {'BB1': BB1, 'BB2': BB2, 'BB3': BB3, 'a1': a1, 'a2': a2, 'a3': a3}
                #print(BB_values)
                atom_numbers = ','.join(f'({a[0]},{a[1]},{a[2]})' for a in angletype[2])

                for n in range(1, 4):
                    # 在这里进行变量的更新
                    if BB_values[f'BB{n}'] in ('1BA', '1BB', '1BC', '1BD', '1BE', '1BF', '1BG', '1BH','TLA','TLB','TLC','BUA'):
                        BB_values[f'BB{n}'] = f'1BA_{BB_values[f"a{n}"][0]}'
                    else:
                        BB_values[f'BB{n}'] = f'{BB_values[f"BB{n}"]}_{BB_values[f"a{n}"]}'

                BB1 = BB_values['BB1']
                BB2 = BB_values['BB2']
                BB3 = BB_values['BB3']    

                with open(f'../force_constant_database/angle.txt') as f:
                    lines = f.readlines()
                    found_in_database = False
                    for line in lines:
                        tmp = line.split()            
                        try:
                            cond1 = (BB1 in tmp[0] and BB2 in tmp[1] and BB3 in tmp[2])
                            cond2 = (BB3 in tmp[0] and BB2 in tmp[1] and BB1 in tmp[2])
                            if cond1 or cond2:
                                file.write(f'     {angletype[0][0]:<9} {angletype[0][1]:<9} {angletype[0][2]:<28} {1:<5} {angletype[1]:<18} {tmp[3]}    ;{atom_numbers}\n')
                                found_in_database = True
                                break
                        except IndexError:
                            continue

                    if not found_in_database:
                        for line in lines:
                            tmp = line.split()
                            try:
                                cond3 = (angletype[0][0] in tmp[0] and angletype[0][1] in tmp[1] and angletype[0][2][4] in tmp[2][4]) or (angletype[0][0] in tmp[2] and angletype[0][1] in tmp[1] and angletype[0][2][4] in tmp[0][4])
                                cond4 = (angletype[0][2] in tmp[0] and angletype[0][1] in tmp[1] and angletype[0][0][4] in tmp[2][4]) or (angletype[0][2] in tmp[2] and angletype[0][1] in tmp[1] and angletype[0][0][4] in tmp[0][4])
                                if cond3 or cond4:
                                    file.write(f'     {angletype[0][0]:<9} {angletype[0][1]:<9} {angletype[0][2]:<28} {1:<5} {angletype[1]:<18} {tmp[3]}    ;{atom_numbers} ;BBdata CCC、CNC...\n')
                                    break

                            except IndexError:
                                continue
                       
            file.write('\n[ dihedraltypes ]\n;dihedral\n')

            arepeat = [] #同樣的atomtype

            for a,b in self.dihedral_dic.items():
                a1,a2,a3,a4 = a.split()[0] , a.split()[1], a.split()[2], a.split()[3]
                atom1, atom2, atom3, atom4 =  self.atomtype_dic[f'{a1}'] , self.atomtype_dic[f'{a2}'], self.atomtype_dic[f'{a3}'], self.atomtype_dic[f'{a4}']
                atomtype = b.get('atomtype')
                atomtype_verse = [' '.join(item.split()[::-1]) for item in atomtype]

                if atomtype_verse in arepeat or atomtype in arepeat:
                    type = f';type{len(arepeat)//2}'
                    self.dihedral_dic[a]['type'] = type
                    self.dihedral_dic[a]['datomtype'] =  atom1, atom2, atom3, atom4
                    #file.write(f'     {atom1:<7} {atom2:<7} {atom3:<7} {atom4:<7} {3:<5} 0.0 0.0 0.0 0.0 0.0 0.0  {type:<10} ({a})')

                if atomtype_verse not in arepeat and atomtype not in arepeat:
                    arepeat.append(atomtype)
                    arepeat.append(atomtype_verse)
                    type = f';type{len(arepeat)//2}'
                    self.dihedral_dic[a]['type'] = type
                    self.dihedral_dic[a]['datomtype'] =  atom1, atom2, atom3, atom4
                    #file.write(f'     {atom1:<7} {atom2:<7} {atom3:<7} {atom4:<7} {3:<5} 0.0 0.0 0.0 0.0 0.0 0.0  {type:<10} ({a})')


            itp_dihedral_dic = {}
            datomtyperepeat = []

            for atomnum , value in self.dihedral_dic.items():
                da1 ,da2 ,da3 ,da4 = value['datomtype'][0] ,value['datomtype'][1] ,value['datomtype'][2] ,value['datomtype'][3]

                if value['datomtype']  in datomtyperepeat:
                    try:
                        itp_dihedral_dic[value['datomtype']]['atomnumber'].append(atomnum)
                    except KeyError:
                        itp_dihedral_dic[(da4,da3,da2,da1)]['atomnumber'].append(atomnum)

                if value['datomtype'] not in datomtyperepeat:
                    itp_dihedral_dic[value['datomtype']] = {'atomnumber' : [atomnum], 'type' : value['type']}
                    datomtyperepeat.append((da1,da2,da3,da4))
                    datomtyperepeat.append((da4,da3,da2,da1))

            for atomtpye , info in itp_dihedral_dic.items():
                atom_numbers_list = info['atomnumber']
                atom_numbers_string = ','.join(['(' + ','.join(item.split()) + ')' for item in atom_numbers_list])
                #print(atom_numbers_string)
                file.write(f'     {atomtpye[0]:<7} {atomtpye[1]:<7} {atomtpye[2]:<7} {atomtpye[3]:<7} {3:<5} 0.0 0.0 0.0 0.0 0.0 0.0  {info["type"]:<10} {atom_numbers_string:<10}\n')
            


            file.write('\n;proper\n')       
            for propertype in proper_itp:
                BB_values = {}
                for i in range(4):
                    BB_values[f'BB{i+1}'] = propertype[0][i].split('_')[0]
                    BB_values[f'a{i+1}'] = propertype[0][i].split('_')[1]
                #原本的BB+atomtype
                BB1_O , BB2_O, BB3_O, BB4_O = propertype[0][0], propertype[0][1], propertype[0][2], propertype[0][3]
                #原本的BB
                BB1_B , BB2_B, BB3_B, BB4_B = BB_values['BB1'], BB_values['BB2'], BB_values['BB3'], BB_values['BB4']

                #原本BB的園子是甚麼
                BB1_a , BB2_a, BB3_a, BB4_a = BB_values['a1'][0], BB_values['a2'][0], BB_values['a3'][0], BB_values['a4'][0]

                #替換成Benzene的名字對應database
                for n in range(1, 5):
                    BB = BB_values[f'BB{n}']
                    a = BB_values[f'a{n}']
                    if BB in ('1BA', '1BB', '1BC', '1BD', '1BE', '1BF', '1BG', '1BH','TLA','TLB','TLC','BUA'):
                        BB_values[f'BB{n}'] = f'1BA_{a[0]}'
                    else:
                        BB_values[f'BB{n}'] = f'{BB_values[f"BB{n}"]}_{BB_values[f"a{n}"]}'

                #替換成1BA的BB+atomtype
                BB1, BB2, BB3, BB4 = BB_values['BB1'], BB_values['BB2'], BB_values['BB3'], BB_values['BB4']
                atom_numbers = ','.join(f'({a[0]},{a[1]},{a[2]},{a[3]})' for a in propertype[2])

                for a in propertype[2]:
                    #proper各自屬於甚麼BB
                    molecule_BB_1 , molecule_BB_2, molecule_BB_3, molecule_BB_4 = self.BB_dic[a[0]], self.BB_dic[a[1]], self.BB_dic[a[2]], self.BB_dic[a[3]]
                #print(BB1_O , BB2_O, BB3_O, BB4_O)
                #print(molecule_BB_1 , molecule_BB_2, molecule_BB_3, molecule_BB_4)

                with open(f'../force_constant_database/proper.txt') as f:
                    lines = f.readlines()
                    atom_numbers_annotation = ';' + atom_numbers
                    for line in lines:
                        tmp = line.split()
                        found_in_database = False
                        #自己BB內para的條件

                        try:
                            cond1 = BB2_B in ('TLB','TLC','BNB') and BB3_B in ('TLB','TLC','BNB')
                            cond2 = BB1_a != 'H' and BB4_a != 'H' and propertype[1] == '0.000' 
                            cond3 = ((molecule_BB_1 == molecule_BB_2 == molecule_BB_3) or (molecule_BB_2 == molecule_BB_3 == molecule_BB_4)) 
                            cond4 = (BB2_B in tmp[1] and BB3_B in tmp[2]) 
                            if  cond1 and cond2 and cond3 and cond4:
                                file.write(f'     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18} {2:<5} {propertype[1]:<18} {tmp[5]}    {atom_numbers_annotation:<10} ;BB inside para position\n')
                                found_in_database = True
                                break
                        except IndexError:
                            continue


                    #BB相接para的條件
                    if not found_in_database:
                        for line in lines:
                            tmp = line.split()        
                            if molecule_BB_1 != molecule_BB_4 and molecule_BB_2 == molecule_BB_3 and molecule_BB_1 != molecule_BB_4 and molecule_BB_4 != molecule_BB_2 and molecule_BB_4 != molecule_BB_3:
                                try: 
                                    condpara1 = (BB2_O in tmp[1] and BB3_O in tmp[2]) or (BB3_O in tmp[1] and BB2_O in tmp[2])
                                    condpara2 = '1BA_C' in tmp[0] and  '1BA_C' in tmp[3]
                                    if condpara1 and condpara2:
                                        file.write(f'     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18} {2:<5} {propertype[1]:<18} {tmp[5]:<10}    {atom_numbers_annotation} ;para position\n')
                                        found_in_database = True
                                        break
                                except IndexError:
                                    continue

                    #換成Benzene找
                    if not found_in_database:
                        for line in lines:
                            tmp = line.split()
                            try:
                                cond1 = (BB1 in tmp[0] and BB2 in tmp[1] and BB3 in tmp[2] and BB4 in tmp[3])
                                cond2 = (BB4 in tmp[0] and BB3 in tmp[1] and BB2 in tmp[2] and BB1 in tmp[3])
                                if cond1 or cond2:
                                    file.write(f'     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18} {2:<5} {propertype[1]:<18} {tmp[5]:<10}    {atom_numbers_annotation} \n')
                                    found_in_database = True
                                    break
                            except IndexError:
                                continue
                        '''                                          
                        try:
                            cond1 = (BB1 in tmp[0] and BB2 in tmp[1] and BB3 in tmp[2] and BB4 in tmp[3])
                            cond2 = (BB4 in tmp[0] and BB3 in tmp[1] and BB2 in tmp[2] and BB1 in tmp[3])
                            if cond1 or cond2:
                                file.write(f'     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18} {1:<5} {propertype[1]:<18} {tmp[5]}    ;{atom_numbers}\n')
                                found_in_database = True
                                break
                        except IndexError:
                            continue
                        ''' 
                    #找CCCC、CNCC
                    if not found_in_database:
                        for line in lines:
                            tmp = line.split()
                            try:
                                cond3 = (propertype[0][0] in tmp[0] and propertype[0][1] in tmp[1] and propertype[0][2] in tmp[2] and propertype[0][3][4] in tmp[3][4]) or (propertype[0][0] in tmp[3] and propertype[0][1] in tmp[2] and propertype[0][2] in tmp[1] and propertype[0][3][4] in tmp[0][4]) 
                                cond4 = (propertype[0][3] in tmp[0] and propertype[0][2] in tmp[1] and propertype[0][1] in tmp[2] and propertype[0][0][4] in tmp[3][4]) or (propertype[0][3] in tmp[3] and propertype[0][2] in tmp[2] and propertype[0][1] in tmp[1] and propertype[0][0][4] in tmp[0][4])
                                if cond3 or cond4 :                                 
                                    file.write(f'     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18} {2:<5} {propertype[1]:<18} {tmp[5]:<10}    {atom_numbers_annotation} ;BBdata CCCC、CNCC...\n')
                                    found_in_database = True
                                    break

                            except IndexError:
                                continue
                    
                    if not found_in_database:
                        file.write(f'\n     {propertype[0][0]:<9} {propertype[0][1]:<9} {propertype[0][2]:<9} {propertype[0][3]:<18}   ;{atom_numbers:<10} ;not found in database\n\n')

            file.write('\n')
            for BB_name in set(self.BB_names):
                file.write(f'#include "./itp/{BB_name}.itp"\n')

    def arrange_rtp(self):
        for n in range(0, len(self.new_BBname)):
            with open(f"{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.new_BBname[n]}.rtp", 'r') as rtp_file:
                lines = rtp_file.readlines()
                for index, line in enumerate(lines):
                    if '[ atoms ]' in line:
                        start = index + 1 
                    if '[ bonds ]' in line:
                        end = index 
                        break

                # Iterate through the lines again to align them
                for i in range(start, end):
                    # If the line contains two specific items, align it
                    if '[ atoms ]' not in lines[i] and '[ bonds ]' not in lines[i]:
                        items = lines[i].split()
                        if len(items) >= 2:  # Ensure there are at least two items
                            items[0] = f"    {items[0]:<6}"  # Left align first item to width < 5 with leading spaces
                            items[1] = f"{items[1]:<6}\n"  # Left align second item to width < 5
                            lines[i] = ''.join(items)  # Join the items back into a line


                impropers_found = False
                end = len(lines)      
                for index, line in enumerate(lines):
                    if '[ bonds ]' in line:
                        start = index + 1 
                    if '[ impropers ]' in line:
                        end = index
                        impropers_found = True 
                    if '[ dihedrals ]' in line and not impropers_found:
                        end = index

                for i in range(start, end):
                    
                    # If the line contains two specific items, align it
                    if  '[ bonds ]' not in lines[i] and '[ impropers ]' not in lines[i]:
                        items = lines[i].split()

                        if len(items) >= 2:  
                            items[0] = f"    {items[0]:<6}"  
                            items[1] = f"{items[1]:<6}\n"
                            lines[i] = ''.join(items)    

                start = 0
                end = 0
                impropers_found = False
                for index, line in enumerate(lines):
                    if '[ impropers ]' in line:
                        start = index + 1
                        end = len(lines)
                        impropers_found = True
                    if '[ dihedrals ]' in line and not impropers_found:
                        start = index + 1
                        end = len(lines)
                        break

                for i in range(start, end):
                    # If the line contains two specific items, align it
                    if  '[ impropers ]' not in lines[i] and '[ dihedrals ]' not in lines[i]:
                        items = lines[i].split()
                        if len(items) >= 4:  
                            items[0] = f"    {items[0]:<6}"  
                            items[1] = f"{items[1]:<6}"
                            items[2] = f"{items[2]:<6}"
                            items[3] = f"{items[3]:<6}\n" 
                            lines[i] = ''.join(items)     
            
            with open(f"{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.new_BBname[n]}.rtp", 'w') as rtp_file:
                rtp_file.writelines(lines)

    def  rtp_add_dihedrals(self):
        for atom_number,value in self.dihedral_dic.items():
            a1 , a2 , a3 , a4 = atom_number.split()[0] , atom_number.split()[1] , atom_number.split()[2] , atom_number.split()[3]
            rtp_name = self.BB_dic[str(a2)]
            atom1 = self.atom_name_dic[str(a1)]
            atom2 = self.atom_name_dic[str(a2)]
            atom3 = self.atom_name_dic[str(a3)]
            atom4 = self.atom_name_dic[str(a4)]
            with open(f'{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/{rtp_name}.rtp', '+r') as f:
                lines = f.readlines()
                found_dihedrals = False
                for n, line in enumerate(lines):
                    if '[ dihedrals ]' in line:
                        found_dihedrals = True     
                if not found_dihedrals :
                    lines.append('\n  [ dihedrals ]\n')
                # Calculate the indentation of start + 3 line
                end = len(lines)

                lines.insert(end , f"{' ' *4} {atom1:<6} {atom2:<8} +{atom3:<8} +{atom4:<6}\n")
            
                f.seek(0)
                f.writelines(lines)

    def rtp_add_impropers(self):
        proper_atomname = []
        proper_atomtype = []
        proper_number = []
        for n in range(len(self.dihedral_number_list)):
            Angle_proper = []
            for i in self.molecule_stk.get_bonds():
                link1 = self.dihedral_number_list[n][1]-1
                link2 = self.dihedral_number_list[n][2]-1
                cond1 = (i.get_atom1().get_id() != link1 or i.get_atom2().get_id() != link2)
                cond2 = (i.get_atom1().get_id() != link2 or i.get_atom2().get_id() != link1)
                cond3 = i.get_atom1().get_id() == link1 or i.get_atom2().get_id() == link1
                cond4 = i.get_atom1().get_id() == link2 or i.get_atom2().get_id() == link2

                if cond1 and cond2:
                    if cond3 or cond4:
                        if i.get_atom1().get_id()+1 == link1+1:
                            Angle_proper.append((f'{link2+1}',f'{link1+1}', f'{i.get_atom2().get_id() + 1}'))
                        if i.get_atom2().get_id()+1 == link1+1:
                            Angle_proper.append((f'{link2+1}',f'{link1+1}', f'{i.get_atom1().get_id() + 1}'))
                        if i.get_atom1().get_id()+1 == link2+1:
                            Angle_proper.append((f'{link1+1}',f'{link2+1}', f'{i.get_atom2().get_id() + 1}'))
                        if i.get_atom2().get_id()+1 == link2+1:
                            Angle_proper.append((f'{link1+1}',f'{link2+1}', f'{i.get_atom1().get_id() + 1}'))

            for n,atom in enumerate(Angle_proper):
                for i in self.molecule_stk.get_bonds():
                    link0 = atom[0]
                    link1 = atom[1]
                    link2 = atom[2]
                    cond1 = i.get_atom1().get_id() != int(link1)-1 and i.get_atom2().get_id() != int(link1)-1
                    cond2 = i.get_atom1().get_id() == int(link2)-1 or i.get_atom2().get_id() == int(link2)-1
                    if cond1:
                        if cond2:
                            if i.get_atom1().get_id()+1 == int(link2):
                                proper_number.append((int(link0), int(link1), int(link2), i.get_atom2().get_id() + 1))
                                proper_atomtype.append((self.atomtype_dic[f'{link0}'],self.atomtype_dic[f'{link1}'],self.atomtype_dic[f'{link2}'], self.atomtype_dic[f'{i.get_atom2().get_id() + 1}']))
                                proper_atomname.append((f"{self.BB_dic[link0]}_{self.atom_name_dic[link0]}",f"{self.BB_dic[link1]}_{self.atom_name_dic[link1]}",f"{self.BB_dic[link2]}_{self.atom_name_dic[link2]}", f"{self.BB_dic[str(i.get_atom2().get_id() + 1)]}_{self.atom_name_dic[str(i.get_atom2().get_id() + 1)]}"))
                            if i.get_atom2().get_id()+1 == int(link2):
                                proper_number.append((int(link0), int(link1), int(link2), i.get_atom1().get_id() + 1))
                                proper_atomtype.append((self.atomtype_dic[f'{link0}'],self.atomtype_dic[f'{link1}'],self.atomtype_dic[f'{link2}'], self.atomtype_dic[f'{i.get_atom1().get_id() + 1}']))
                                proper_atomname.append((f"{self.BB_dic[link0]}_{self.atom_name_dic[link0]}",f"{self.BB_dic[link1]}_{self.atom_name_dic[link1]}",f"{self.BB_dic[link2]}_{self.atom_name_dic[link2]}", f"{self.BB_dic[str(i.get_atom1().get_id() + 1)]}_{self.atom_name_dic[str(i.get_atom1().get_id() + 1)]}"))

        for rtp,number in zip(proper_atomname,proper_number):
            atom1 = rtp[0].split('_')[1]
            atom2 = rtp[1].split('_')[1]
            atom3 = rtp[2].split('_')[1]
            atom4 = rtp[3].split('_')[1]
            rtpname1 = rtp[0].split('_')[0]
            rtpname2 = rtp[1].split('_')[0]
            rtpname3 = rtp[2].split('_')[0]
            rtpname4 = rtp[3].split('_')[0]
        
            with open(f'{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/{rtpname2}.rtp','r+') as f:
                lines = f.readlines()
                found_impropers = False  # 添加一个标志来指示是否找到了 '[ impropers ]'
                for n, line in enumerate(lines):
                    if '[ impropers ]' in line:
                        found_impropers = True
                        start = n+1
                if not found_impropers:  # 如果没有找到 '[ impropers ]'
                    lines.append('\n[ impropers ]\n')  # 在列表末尾添加 '[ impropers ]\n'
                    start = len(lines)  # 设置开始行为列表的最后一行
                end = len(lines)

                start_plus_3_indentation = len(lines[start + 1]) - len(lines[start + 1].lstrip()) 
                

                if number[0] < number[1] and rtpname2 == rtpname4:
                    lines.insert(end - 1, f"{' ' * start_plus_3_indentation} -{atom1:<6} {atom2:<8} {atom3:<8} {atom4:<6} \n")
                if number[0] > number[1] and rtpname2 == rtpname4:
                    lines.insert(end - 1, f"{' ' * start_plus_3_indentation} +{atom1:<6} {atom2:<8} {atom3:<8} {atom4:<6} \n")
                if rtpname2 != rtpname4:
                    if number[0] < number[1] and number[2] < number[3]:
                        lines.insert(end - 1, f"{' ' * start_plus_3_indentation} -{atom1:<6} {atom2:<8} {atom3:<8} +{atom4:<6} \n")
                    if number[0] > number[1] and number[2] > number[3]:
                        lines.insert(end - 1, f"{' ' * start_plus_3_indentation} +{atom1:<6} {atom2:<8} {atom3:<8} -{atom4:<6} \n")
                    if number[0] < number[1] and number[2] > number[3]:
                        lines.insert(end - 1, f"{' ' * start_plus_3_indentation} -{atom1:<6} {atom2:<8} {atom3:<8} -{atom4:<6} \n")
                    if number[0] > number[1] and number[2] < number[3]:
                        lines.insert(end - 1, f"{' ' * start_plus_3_indentation} +{atom1:<6} {atom2:<8} {atom3:<8} +{atom4:<6} \n")
                
                f.seek(0)
                f.writelines(lines) 

    def rtp_add_bonds(self):
        #-的
        for dihedral in self.dihedral_number_list:
            d1 = dihedral[1]
            d2 = dihedral[2]

            with open(f'{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.BB_dic[str(d1)]}.rtp', 'r+') as f:
                lines = f.readlines()
                for n, line in enumerate(lines):
                    if '[ bonds ]' in line:
                        start = n + 1
                    elif '[ impropers ]' in line:
                        end = n
                        break  # Assuming you want to stop loop once you found the end
                else:
                    end = len(lines)
                # Calculate the indentation of start + 3 line
                start_plus_3_indentation = len(lines[start + 1]) - len(lines[start + 1].lstrip())

                # Insert print statement with preceding lines at the desired location
                lines.insert(end - 1, f"{' ' * start_plus_3_indentation}{self.atom_name_dic[str(d1)]:<5} +{self.atom_name_dic[str(d2)]:<5}\n")
                # Now, write the modified content back to the file
                f.seek(0)
                f.writelines(lines)

        for dihedral in self.dihedral_number_list:
            d1 = dihedral[2]
            d2 = dihedral[1]

            with open(f'{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.BB_dic[str(d1)]}.rtp', 'r+') as f:
                lines = f.readlines()
                for n, line in enumerate(lines):
                    if '[ bonds ]' in line:
                        start = n + 1
                    elif '[ impropers ]' in line:
                        end = n
                        break  # Assuming you want to stop loop once you found the end
                else:
                    end = len(lines)
                # Calculate the indentation of start + 3 line
                start_plus_3_indentation = len(lines[start + 1]) - len(lines[start + 1].lstrip())

                # Insert print statement with preceding lines at the desired location
                lines.insert(end - 2, f"{' ' * start_plus_3_indentation}{self.atom_name_dic[str(d1)]:<5} -{self.atom_name_dic[str(d2)]:<5}\n")
                # Now, write the modified content back to the file
                f.seek(0)
                f.writelines(lines)




    def create_FF_folders(self):  #建立FF所需的資料夾
        folders = ["FORCEFIELD", "QM_GASPHASE/xtb", "SCANS", "STRUCTURE"] 
        for folder in folders:
            if folder == "FORCEFIELD":
                path = f"molecule/{self.molecule_name}/{folder}/{self.molecule_name}.ff"
            else:
                path = f"molecule/{self.molecule_name}/{folder}"
            os.makedirs(path, exist_ok=True)
        self.molecule_folder_path = f"molecule/{self.molecule_name}"
        with open (f'{self.molecule_folder_path}/STRUCTURE/BB_names_connector.txt','w') as file:
            file.write(f'molecule_name = "{self.molecule_name}"\n')
            file.write(f'BB_name = {self.BB_names}\n')
            file.write(f'BB_connector = [')
            for index,factory in enumerate(self.BB_connector):
                if 'Bromo' in str(factory):
                    file.write('Br')
                if 'Iodo' in str(factory):
                    file.write('I')
                if 'Fluoro' in str(factory):
                    file.write('F')
                
                # 如果不是最后一个元素，则写入逗号
                if index != len(self.BB_connector) - 1:
                    file.write(',')
            
            # 写入右方括号
            file.write(']\n')

        
    def STK_create_building_blocks(self):
        self.BB_stk_all = []
        with open('BB_info.json') as f:
            BB_info = json.load(f)
        for BB_name_n,BB_name in enumerate(self.BB_names):
            for molecule_info in BB_info.values():
                for BB in molecule_info["BB"]:
                    if BB == BB_name:
                        smiles =  molecule_info["smiles"]
                        if 'Br' in  molecule_info['BB'][BB_name]:
                            Br_num = molecule_info['BB'][BB_name]['Br']

                        else:
                            Br_num = None

                        if 'I' in  molecule_info['BB'][BB_name]:
                            I_num = molecule_info['BB'][BB_name]['I']

                        else:
                            I_num = None

                        if 'F' in  molecule_info['BB'][BB_name]:
                            F_num = molecule_info['BB'][BB_name]['F']

                        else:
                            F_num = None
                        
                        building_block = stk.BuildingBlock(smiles)
                        mol_rdkit =  building_block.to_rdkit_mol()
                        atom = list(mol_rdkit.GetAtoms())
                        for n, a in enumerate(atom):
                            if Br_num is not None and n+1 in Br_num: #變數
                                a.SetAtomicNum(35) #Br
                            if I_num is not None and n+1 in I_num:
                                a.SetAtomicNum(53) #I 
                            if F_num is not None and n+1 in F_num:
                                a.SetAtomicNum(9) #F

                        atom_num = [ a.GetAtomicNum() for a in atom]

                        BB_stk = stk.BuildingBlock.init_from_rdkit_mol(mol_rdkit,self.BB_connector[BB_name_n])
                        self.BB_stk_all.append(BB_stk)

    
    def create(self,B0, B1):
       
        #建立分子
        c1 = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(B0, B1),
                repeating_unit='AB',
                num_repeating_units=1,
                optimizer=stk.MCHammer(),
            ),
        )
        
        connect_atom = []
        bond_infos = [ x for x in c1.get_bond_infos()]

        for b in bond_infos:
            BBinfo = b.get_building_block_id()
            if BBinfo is None:
                bondinfo = b.get_bond()
                #print(bondinfo)
                connect_atom.append(bondinfo.get_atom1())
                connect_atom.append(bondinfo.get_atom2())        

        for i    in range(0,len(connect_atom)-1,2):
            if i+1 < len(connect_atom):
                for b in bond_infos:
                    bond = b.get_bond()
                    cond1 = bond.get_atom1() == connect_atom[i] #or  bond.get_atom2() == connect_atom[i]
                    cond2 = bond.get_atom1() != connect_atom[i+1] and  bond.get_atom2() != connect_atom[i+1]
                    if cond1 and cond2:
                        n1 = bond.get_atom2().get_id()+1 #python從0開始，編號要+1
                        n2 = bond.get_atom1().get_id()+1
                        
                    cond3 =  bond.get_atom1() == connect_atom[i+1] #or  bond.get_atom2() == connect_atom[i+1]
                    cond4 = bond.get_atom1() != connect_atom[i] and  bond.get_atom2() != connect_atom[i]
                    if cond3 and cond4:
                        n3 = bond.get_atom1().get_id()+1 #python從0開始，編號要+1
                        n4 = bond.get_atom2().get_id()+1   

        self.dihedral_number_list.append((n1,n2,n3,n4))

                    
        rdkit_c1 = c1.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_c1)
        rdkit.MMFFOptimizeMolecule(rdkit_c1,maxIters=500)
        smiles1 = rdkit.MolToSmiles(c1.to_rdkit_mol(),allHsExplicit=True)

        self.molecule_smile.update({''.join(self.molecule_name):{'smiles':smiles1}})


        c1 = c1.with_position_matrix(
            position_matrix=rdkit_c1.GetConformer().GetPositions(),
        )
        stk.XyzWriter().write(c1,f'{self.molecule_folder_path}//QM_GASPHASE//xtb//{self.molecule_name}.xyz')

        self.final_molecule = c1
        '''
        for a in self.final_molecule.get_building_blocks():
            for b in a.get_functional_groups():
                print(b)
        '''

        self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.BromoFactory()],)

        '''
        if 'Br' not in str(self.molecule_smile[self.molecule_name]):
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.IodoFactory()],)
            self.number += 1
            self.count.append(self.number)

        if 'Br' not in str(self.molecule_smile[self.molecule_name]) and 'I' not in str(self.molecule_smile[self.molecule_name]):
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.FluoroFactory()],)
            self.number += 1
            self.count.append(self.number)          
        
        if str(self.molecule_smile[self.molecule_name]).count('Br') >= 2:            
            print('Br:2')
            bromo_factory = stk.BromoFactory()        
            bromo_groups = tuple(bromo_factory.get_functional_groups(molecule=self.molecule_stk))
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[bromo_groups[0]])
            
                  
        if str(self.molecule_smile[self.molecule_name]).count('I') >= 2:            
            print('I:2')
            Iodo_factory = stk.IodoFactory()
            Iodo_groups = tuple(Iodo_factory.get_functional_groups(molecule=self.molecule_stk))
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[Iodo_groups[1]])                        
        
        if str(self.molecule_smile[self.molecule_name]).count('F') >= 2:            
            print('F:2')
            Fluoro_factory = stk.FluoroFactory()
            Fluoro_groups = tuple(Fluoro_factory.get_functional_groups(molecule=self.molecule_stk))
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[Fluoro_groups[0]])   
        '''
        for a in self.molecule_stk.get_functional_groups():
            deleter_ids = a.get_deleter_ids()
            for b in deleter_ids:                   
                deleter_id = b + 1
                for i in range(len(self.dihedral_number_list)):
                    self.dihedral_number_list[i] = tuple(x-1 if x > deleter_id else x for x in self.dihedral_number_list[i])
                else:
                    pass
        
        for a in self.molecule_stk.get_functional_groups():
            deleter_ids = a.get_deleter_ids()
            for b in deleter_ids:                   
                deleter_id = b + 1
                for i in range(len(self.dihedral_number_list)):
                    self.dihedral_number_list[i] = tuple(x-1 if x > deleter_id else x for x in self.dihedral_number_list[i])
                else:
                    pass

        def read_xyz(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            # 提取 atom 資訊 (忽略前兩行)
            atoms = [line.split() for line in lines[2:]]
            
            # 記錄 Br, I, F 的原子
            found_atoms = {'Br': [], 'I': [], 'F': []}
            
            for index, atom in enumerate(atoms, start=1):  # index 作為原子編號
                element = atom[0]  # 原子的元素符號
                if element in found_atoms:
                    found_atoms[element].append(index)
            
            # 找出最大原子編號的元素
            max_atom = None
            max_index = 0
            for element, indices in found_atoms.items():
                if indices:
                    max_for_element = max(indices)
                    if max_for_element > max_index:
                        max_index = max_for_element
                        max_atom = element
            
            return found_atoms, max_atom, max_index


        # 使用範例
        file_path = f'{self.molecule_folder_path}//QM_GASPHASE//xtb//{self.molecule_name}.xyz'
        found_atoms, max_atom, max_index = read_xyz(file_path)

        if found_atoms:
            print(f"Found atoms: {found_atoms}")
            if max_atom:
                print(f"Atom with the largest number is {max_atom} with index {max_index}")
        else:
            print("No Br, I, or F atoms found.")

        if max_atom is None:
            self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1, [stk.BromoFactory()])

        if  max_atom is not None:
            if len(found_atoms[max_atom]) == 1 and max_atom is not None:
                if max_atom == 'Br':
                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.BromoFactory()],)
                if max_atom == 'I':
                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.IodoFactory()],)
                if max_atom == 'F':
                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[stk.FluoroFactory()],)

        if  max_atom is not None:
            if len(found_atoms[max_atom]) > 1 :
                if max_atom == 'Br':
                    bromo_factory = stk.BromoFactory()        
                    bromo_groups = tuple(bromo_factory.get_functional_groups(molecule=self.molecule_stk))
                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[bromo_groups[-1]])
                if max_atom == 'I':
                    Iodo_factory = stk.IodoFactory()
                    Iodo_groups = tuple(Iodo_factory.get_functional_groups(molecule=self.molecule_stk))

                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[Iodo_groups[-1]]) 
                if max_atom == 'F':
                    Fluoro_factory = stk.FluoroFactory()
                    Fluoro_groups = tuple(Fluoro_factory.get_functional_groups(molecule=self.molecule_stk))
                    self.molecule_stk = stk.BuildingBlock.init_from_molecule(c1,[Fluoro_groups[-1]])


        #print(self.dihedral_number_list)
        #分子smiles寫入json檔
        with open(f"{self.molecule_folder_path}//STRUCTURE//{self.molecule_name}.json","w") as outfile:
            json.dump(self.molecule_smile, outfile, indent = 3)   

        return self.molecule_stk


    def build_molecule(self):
        for n in range(len(self.BB_stk_all)):
            if n == 0:
                self.molecule_stk = self.create(self.BB_stk_all[0], self.BB_stk_all[1])
            elif 0 < n < len(self.BB_stk_all) - 1:
                self.molecule_stk = self.create(self.molecule_stk, self.BB_stk_all[n + 1])
        
        return self.molecule_stk

    def show_molecule(self):
        rdkit_mol = self.molecule_stk.to_rdkit_mol()
        viewer = MoleculeViewer(rdkit_mol)
        viewer.show()
        return viewer
    


    def update_dihedral_number_pdb_create_SCAN(self):
        if not os.path.exists(f'molecule/{self.molecule_name}/QM_GASPHASE/{self.molecule_name}.log'):
            print(f"無法辨別找出對稱的dihedral請先完成gaussian計算")
            raise SystemExit


        n_atoms = self.molecule_stk.get_num_atoms()
        with open(f'{self.molecule_folder_path}/QM_GASPHASE/{self.molecule_name}.log', 'r') as file:
            for index, line in enumerate(file):
                if "Input orientation" in line:
                    start = index + 5
                    end = index + n_atoms + 5
                    break

            elements = []
            coords = []
            file.seek(0)  # 將檔案指標重設為開始位置
            for i, line in enumerate(file):
                if start <= i < end:
                    element = line.split()[1]
                    elements.append(int(element))
                    coord = line.split()[3:6]
                    coords.append([float(i) for i in coord])


        with open(f'{self.molecule_folder_path}/QM_GASPHASE/xtb/wbo','r') as file:

            b_orders = [[0] * n_atoms for _ in range(n_atoms)]  # 初始化为0
            file.seek(0)
            lines = file.readlines()
            for line in lines:
                if not line.split():
                    break
                atom_1 = int(line.split()[0])-1
                atom_2 = int(line.split()[1])-1
                bond_order = float(line.split()[2])

                
                b_orders[atom_1][atom_2] = bond_order

            elements = check_type_and_shape(elements, 'elements', int, (n_atoms,))
            coords = check_type_and_shape(coords, 'coords', float, (n_atoms, 3))
            b_orders = check_type_and_shape(b_orders, 'b_orders', float, (n_atoms, n_atoms))
            n_bonds = b_orders.sum(axis=1).round().astype(int)

            elements = np.array(elements)
            coords = np.asfarray(coords, float)
            coords = np.reshape(coords, (-1, 3))
            coords = np.array(coords)* Bohr


        graph = nx.Graph()
        for i_idx, i_elem in enumerate(elements):
            graph.add_node(i_idx,  elem=i_elem ,n_bonds=n_bonds[i_idx],
                                coords=coords[i_idx])
            # add bonds
            for j_idx, j_elem in enumerate(elements):
                b_order = b_orders[i_idx, j_idx]
                if b_order > 0.3:
                    id1, id2 = sorted([i_elem, j_elem])
                    b_order_half_rounded = np.round(b_order*2)/2
                    vec = coords[i_idx] - coords[j_idx]
                    dist = np.sqrt((vec**2).sum())
                    graph.add_edge(i_idx, j_idx, vector=vec, length=dist, order=b_order,
                                        type=f'{id1}({b_order_half_rounded}){id2}', n_rings=0)


        rings = nx.minimum_cycle_basis(graph) #循環
        rings3 = [r for r in rings if len(r) == 3]

        for i in range(n_atoms):
            graph.nodes[i]['n_ring'] = sum([i in ring for ring in rings]) #儲存原子的環數
        #
        for atoms in graph.edges: #全部的鍵結
            ring_members = [set(atoms).issubset(set(ring)) for ring in rings]
            graph.edges[atoms[0],atoms[1]]['n_rings'] = sum(ring_members) #鍵結再幾個環裡面
            graph.edges[atoms[0],atoms[1]]['in_ring'] = any(ring_members)
            graph.edges[atoms[0],atoms[1]]['in_ring3'] = any(set(atoms).issubset(set(ring))for ring in rings3)

        atom_ids = [[] for _ in range(n_atoms)]

        bonds, angles, dihedrals = [], [], []
        
        #辨認topology
        for i in range(n_atoms):
            neighbors = nx.bfs_tree(graph, source=i,  depth_limit=3).nodes
            for n in neighbors:
                paths = nx.all_simple_paths(graph, i, n, cutoff=3)
                for path in paths:

                    if len(path) == 2 and path[0] < path[1]: #後面的條件是保證不會重複辨認同樣的情況
                        bonds.append(path)
                    elif len(path) == 3 and path[0] < path[2]:
                        angles.append(path)
                    elif len(path) == 4 and path[1] < path[2]:
                        dihedrals.append(path)

        dihedrals.sort(key=lambda x: [x[1], x[2], x[0], x[3]]) #dihedral的第二項最小為優先順序，再來是第三項，以此類推

        #self.bonds = bonds
        #self.angles = angles
        #self.dihedrals = dihedrals
        #print(dihedrals)
        
        for i in range(n_atoms):
            neighbors = nx.bfs_tree(graph, i, depth_limit=4).nodes
            for n in neighbors:
                paths = nx.all_simple_paths(graph, i, n, cutoff=4)
                for path in map(nx.utils.pairwise, paths):
                    types = [graph.edges[edge[0],edge[1]]['type'] for edge in path]
                    atom_ids[i].append("-".join(types))
                atom_ids[i].sort()
        if 4 < 0:
            atom_ids = [i for i in range(n_atoms)]

        list = []

        for n in range(n_atoms):
            if n in [item for sub in list for item in sub]:
                continue
            eq = [i for i, a_id in enumerate(atom_ids) if atom_ids[n] == a_id]
            list.append(eq)  

        atomtype_dic = {}
        for n,atoms in enumerate(list):
            for atom in atoms:
                atomtype_dic[str(atom+1)] = n
                
        dihedral = []
        dihedral_atomtype = []
        dihedral_BB = []
        for d in self.dihedral_number_list:
            a1 = []
            a4 = []
            for i in self.molecule_stk.get_bonds():
                link1 = d[1]-1
                link2 = d[2]-1
                a2 = i.get_atom1().get_id()+1
                a3 = i.get_atom2().get_id()+1
                cond1 = (i.get_atom1().get_id() != link1 or i.get_atom2().get_id() != link2)
                cond2 = (i.get_atom1().get_id() != link2 or i.get_atom2().get_id() != link1)
                cond3 = i.get_atom1().get_id() == link1 or i.get_atom2().get_id() == link1
                cond4 = i.get_atom1().get_id() == link2 or i.get_atom2().get_id() == link2
                if cond1 and cond2:
                    if cond3 or cond4:
                        if i.get_atom1().get_id()+1 == link1+1:
                            a1.append(f'{i.get_atom2().get_id()+1}')
                        if i.get_atom2().get_id()+1 == link1+1:
                            a1.append(f'{i.get_atom1().get_id()+1}')
                        if i.get_atom1().get_id()+1 == link2+1:
                            a4.append(f'{i.get_atom2().get_id()+1}')
                        if i.get_atom2().get_id()+1 == link2+1:
                            a4.append(f'{i.get_atom1().get_id()+1}')

            for atom1 in a1:
                for atom4 in a4:
                    #print(f'{atom1} {d[1]} {d[2]} {atom4} ')
                    #print(f'atom type : {atomtype_dic[atom1]} {atomtype_dic[str(d[1])]} {atomtype_dic[str(d[2])]} {atomtype_dic[atom4]}')
                    #print(f'{self.atom_name_dic[str(atom1)]}  {self.atom_name_dic[str(d[1])]}  {self.atom_name_dic[str(d[2])]}  {self.atom_name_dic[str(atom4)]}')
                    #print(f'{self.BB_dic[str(atom1)]}  {self.BB_dic[str(d[1])]}  {self.BB_dic[str(d[2])]}  {self.BB_dic[str(atom4)]}')
                    dihedral.append(f'{atom1} {d[1]} {d[2]} {atom4}')
                    dihedral_atomtype.append(f'{atomtype_dic[atom1]} {atomtype_dic[str(d[1])]} {atomtype_dic[str(d[2])]} {atomtype_dic[atom4]}')
                    dihedral_BB.append(f'{self.BB_dic[str(atom1)]}  {self.BB_dic[str(d[1])]}  {self.BB_dic[str(d[2])]}  {self.BB_dic[str(atom4)]}')

        dihedral_dic = {}
        for a,b,c in zip(dihedral,dihedral_atomtype,dihedral_BB):
            dihedral_dic[a] = {'atomtype' : [b] , 'BB' : [c]}



        for atomnumber, value in dihedral_dic.items():
            atomtype = value['atomtype']

            atomtype_verse = [' '.join(item.split()[::-1]) for item in atomtype]
            
            BB = value.get('BB') 
            for atomnumber1,value1 in dihedral_dic.items():
                atomtype1 = value1['atomtype']
                BB1 = value1.get('BB')
                if (atomtype1 == atomtype_verse or atomtype1 == atomtype) and BB1 not in [v['BB'] for v in self.dihedral_dic.values()]:
                    self.dihedral_dic[atomnumber1] = {'atomtype': atomtype1, 'BB': BB1}   

        repeat = []
        for atom_number,value in self.dihedral_dic.items():
            atomtype = value['atomtype']

            atomtype_verse = [' '.join(item.split()[::-1]) for item in atomtype]
            if atomtype not in repeat and atomtype_verse not in repeat:
                repeat.append(atomtype)
                repeat.append(atomtype_verse)
                a1,a2,a3,a4 = atom_number.split()[0] , atom_number.split()[1] , atom_number.split()[2] , atom_number.split()[3]
                self.dihedral_number_pdb += f'REMARK SCAN  {a1:<3} {a2:<3} {a3:<3} {a4:<3}  0.0 360.0 10.0  3\n'
                

                SCANS_directory = f"{self.molecule_folder_path}/SCANS/{a1}_{a2}_{a3}_{a4}"

            if not os.path.exists(SCANS_directory):
                os.makedirs(SCANS_directory)
        
        #分子的dihedral編號寫入txt檔用來放到.pdb檔內
        with open(f"{self.molecule_folder_path}//STRUCTURE//dihedral_number.txt", 'w') as file:
            file.write(self.dihedral_number_pdb)

    def renameBB_copy(self):
        for name in self.BB_names:
            if name in self.name_counter:
                self.name_counter[name] += 1
                new_name = f"{name[:2]}{self.name_counter[name]}"
                while new_name in self.new_BBname:
                    self.name_counter[name] += 1
                    new_name = f"{name[:2]}{self.name_counter[name]}"
                self.new_BBname.append(new_name)
            else:
                self.name_counter[name] = 0
                self.new_BBname.append(name)
            for name, new_name in zip(self.BB_names, self.new_BBname):
                #source_file = f"..//RTP/{name}.rtp"
                source_file = f"qforce_rtp/{name}.rtp"
                destination_folder = f"{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff"
                destination_file = f"{destination_folder}//{new_name}.rtp"
                #print(f"Copying {source_file} to {destination_file}")
                os.makedirs(destination_folder, exist_ok=True)     
                shutil.copy(source_file, destination_file)

    def replace_BBname_in_rtp(self):
        for n in range(0, len(self.new_BBname)):
            with open(f"{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.new_BBname[n]}.rtp", 'r') as rtp_file:
                rtp_content = rtp_file.readlines()

                for i, line in enumerate(rtp_content):
                    if '[' in line and self.BB_names[n] in line:
                        rtp_content[i] = line.replace(self.BB_names[n], self.new_BBname[n])

            with open(f"{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{self.new_BBname[n]}.rtp", 'w') as rtp_file:
                rtp_file.writelines(rtp_content)   
 
    def pdb_file_information(self):          
        number = 0  
        for n, (a, b) in enumerate(zip(self.new_BBname, self.BB_names)):
            with open(f'{self.molecule_folder_path}//FORCEFIELD//{self.molecule_name}.ff//{a}.rtp', 'r') as rtp_file:
                rtp_content = rtp_file.read()

            # 有沒有讀到atoms
            in_atoms_section = False

            # 讀取rtp檔裡的資訊
            for line in rtp_content.split('\n'):
                if line.strip() == "[ atoms ]":  # 空行刪除是[ atoms ]
                    in_atoms_section = True
                    continue
                elif line.strip().startswith("["):  # 再次出現[
                    in_atoms_section = False

                if in_atoms_section and line.strip() != "":
                    number += 1
                    atom_name = line.strip().split()[0]
                    atom_type = line.strip().split()[1]
                    self.pdb_information.append(f"{atom_name:3} {a:8}{n+1}")
                    self.atom_name.append(f"{atom_name}")

                    self.atom_type.append(f"{atom_type}")
                    self.atom_name_dic[f'{number}'] = atom_name
                    self.BB_dic[f'{number}'] = a

        with open(f'{self.molecule_folder_path}//STRUCTURE//pdb.txt', 'w') as output_file:
            output_file.write('\n'.join(self.pdb_information))

    def establish_atp(self):
        unique_atom_types = list(OrderedDict.fromkeys(self.atom_type))
        with open(f"{self.molecule_folder_path}/FORCEFIELD/{self.molecule_name}.ff/atp.atp", 'w') as file:
            for atom_type in unique_atom_types:
                if '_C' in atom_type:
                    file.write(f"{atom_type:<8}  {12.011}\n")  
                if '_H' in atom_type:                     
                    file.write(f"{atom_type:<8}  {1.008}\n")
                if '_N' in atom_type:
                    file.write(f"{atom_type:<8}  {14.007}\n")
                if '_O' in atom_type:
                    file.write(f"{atom_type:<8}  {15.999}\n")
                if '_S' in atom_type:
                    file.write(f"{atom_type:<8}  {32.06}\n")

    def map_atomnames_to_dihedral(self):
        self.atomtype_dic = {str(i+1): value for i, value in enumerate(self.atom_type)}
        self.dihedral_rtp = []
        for index_tuple in self.dihedral_number_list:
            mapped_names = []
            for index in index_tuple:
                mapped_names.append(f'{self.atom_name[index - 1]}')
            self.dihedral_rtp.append(mapped_names)