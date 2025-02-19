import sys
import stk
import json
import re
import os
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, 
                             QLabel, QLineEdit, QTextEdit, QDialog, QFormLayout, QMessageBox, QHBoxLayout)
from PyQt6.QtCore import QSize
from cut_search import fragment_molecule, search_database
from link import link
from clean import clean_para, clean_all
from convert import info_to_dic
import py3Dmol
from rdkit import Chem


class EditBB(QDialog):
    def __init__(self, bb_name, smiles, is_missing_smiles, parent=None):
        super().__init__(parent)
        self.bb_name = bb_name
        self.is_missing_smiles = is_missing_smiles
        self.initUI(clean_para(smiles))

    def initUI(self, smiles):
        self.setWindowTitle(f"Edit Building Block: {self.bb_name}")
        layout = QVBoxLayout()

        form_layout = QFormLayout()

        self.name_input = QLineEdit(self.bb_name if self.bb_name != "Not Found" else "")
        form_layout.addRow("Molecule", self.name_input)

        self.smiles_input = QLineEdit(smiles)
        self.smiles_input.setEnabled(self.is_missing_smiles)
        form_layout.addRow("SMILES", self.smiles_input)

        self.bb_input = QLineEdit()
        form_layout.addRow("BB Name", self.bb_input)

        layout.addLayout(form_layout)

        # Connectors input
        connectors_layout = QHBoxLayout()
        connectors_layout.addWidget(QLabel("Connectors"))
        
        self.br_input = QLineEdit()
        self.i_input = QLineEdit()
        self.f_input = QLineEdit()
        
        connectors_layout.addWidget(QLabel("Br:"))
        connectors_layout.addWidget(self.br_input)
        connectors_layout.addWidget(QLabel("I:"))
        connectors_layout.addWidget(self.i_input)
        connectors_layout.addWidget(QLabel("F:"))
        connectors_layout.addWidget(self.f_input)
        
        layout.addLayout(connectors_layout)


        if self.is_missing_smiles:
            reminder_label = QLabel("IMPORTANT!! Remember to add the topologies file to the folders 'qforce_itp' and 'qforce_rtp'.")
            layout.addWidget(reminder_label)

        buttons = QHBoxLayout()
        save_button = QPushButton("Save")
        save_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        buttons.addWidget(save_button)
        buttons.addWidget(cancel_button)

        layout.addLayout(buttons)
        self.setLayout(layout)

    def get_data(self):
        connectors = {}
        for connector, input_field in [("Br", self.br_input), ("I", self.i_input), ("F", self.f_input)]:
            if input_field.text():
                try:
                    connectors[connector] = [int(x.strip()) for x in input_field.text().split(',') if x.strip()]
                except ValueError:
                    QMessageBox.warning(self, "Invalid Input", f"Invalid input for {connector}. Please enter comma-separated integers.")
                    return None

        return {
            "name": self.name_input.text(),
            "smiles": self.smiles_input.text(),
            "bb": self.bb_input.text(),
            "connectors": connectors
        }


class AutoFF(QWidget):
    def __init__(self):
        super().__init__()
        self.database = {}  # BB_info.json
        self.initUI()
        self.load_database()
        self.load_bb_dic()

    def load_database(self):
        try:
            with open('BB_info.json', 'r') as f:
                self.database = json.load(f)
        except FileNotFoundError:
            self.database = {}

    def save_database(self):
        with open('BB_info.json', 'w') as f:
            json.dump(self.database, f, indent=4)
        
        try:
            self.bb_dic = info_to_dic('BB_info.json', 'BB_dic.json')
            QMessageBox.information(self, "Information", "Database updated successfully. Check console for any warnings.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to update BB_dic.json: {str(e)}\n\nPlease check the console for more details.")
        
    def load_bb_dic(self):
        try:
            with open('BB_dic.json', 'r') as f:
                self.bb_dic = json.load(f)
        except FileNotFoundError:
            try:
                self.bb_dic = info_to_dic('BB_info.json', 'BB_dic.json')
                QMessageBox.information(self, "Information", "BB_dic.json created. Check console for any warnings.")
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to create BB_dic.json: {str(e)}\n\nPlease check the console for more details.")
                self.bb_dic = {}

    def initUI(self):
        self.resize(QSize(800, 600))  
        self.setWindowTitle('Molecule Builder')

        layout = QVBoxLayout()

        self.name_label = QLabel("Enter the name of your molecule:")
        layout.addWidget(self.name_label)

        self.name_input = QLineEdit(self)
        self.name_input.setFixedHeight(30)
        layout.addWidget(self.name_input)

        self.smiles_label = QLabel("Enter the canonical SMILES of your molecule:")
        layout.addWidget(self.smiles_label)

        self.smiles_input = QTextEdit(self)
        self.smiles_input.setFixedHeight(50)
        layout.addWidget(self.smiles_input)

        self.cut_button = QPushButton('Cut and Search', self)
        self.cut_button.clicked.connect(self.cut_and_search)
        layout.addWidget(self.cut_button)

        self.result_table = QTableWidget(self)
        self.result_table.setFixedHeight(250)  
        layout.addWidget(self.result_table)

        self.script_button = QPushButton('Generate Script', self)
        self.script_button.clicked.connect(self.generate_script)
        self.script_button.setEnabled(False)  
        layout.addWidget(self.script_button)

        self.script_output = QTextEdit(self)
        self.script_output.setFixedHeight(100)
        self.script_output.setReadOnly(False)
        layout.addWidget(self.script_output)

        self.confirm_button = QPushButton('Confirm Info', self)
        self.confirm_button.clicked.connect(self.confirm_action)
        self.confirm_button.setEnabled(False) 
        layout.addWidget(self.confirm_button)

        self.match_output = QTextEdit(self)
        self.match_output.setFixedHeight(50)
        layout.addWidget(self.match_output)

        self.setLayout(layout)

    def cut_and_search(self):
        molecule_name = self.name_input.text()
        smiles = self.smiles_input.toPlainText()

        if smiles and molecule_name:
            fragments = fragment_molecule(smiles)  
            if fragments:
                results = search_database(fragments, "BB_dic.json", "BB_info.json")
                self.show_results(results)
                self.script_button.setEnabled(True)
            else:
                self.result_table.setRowCount(0)
                self.result_table.setColumnCount(1)
                self.result_table.setHorizontalHeaderLabels(["Error"])
                self.result_table.setItem(0, 0, QTableWidgetItem("Invalid SMILES or unable to fragment."))
    
    def show_results(self, results):
        self.result_table.setRowCount(len(results))
        self.result_table.setColumnCount(4)
        self.result_table.setHorizontalHeaderLabels(["BB / Molecule", "Fragment SMILES", "Status", "Action"])
        self.result_table.setColumnWidth(0, 150)
        self.result_table.setColumnWidth(1, 300)
        self.result_table.setColumnWidth(2, 150)
        self.result_table.setColumnWidth(3, 100)

        for i, (bb_name, frag_smiles, status) in enumerate(results):
            self.result_table.setItem(i, 0, QTableWidgetItem(bb_name if bb_name else "Not Found"))
            self.result_table.setItem(i, 1, QTableWidgetItem(frag_smiles))  # Original SMILES with linkers and [*]
            
            if status == True:
                status_text = "Found"
            elif status == "Missing Connectors":
                status_text = "Missing Connectors"
            else:
                status_text = "Missing SMILES"
            
            self.result_table.setItem(i, 2, QTableWidgetItem(status_text))
            
            if status != True:
                edit_button = QPushButton("[+]")
                edit_button.clicked.connect(lambda _, row=i: self.edit_bb(row))
                self.result_table.setCellWidget(i, 3, edit_button)
            else:
                self.result_table.removeCellWidget(i, 3)

    def edit_bb(self, row):
        bb_name = self.result_table.item(row, 0).text()
        frag_smiles = self.result_table.item(row, 1).text()
        status = self.result_table.item(row, 2).text()
        
        is_missing_smiles = (status == "Missing SMILES")
        
        cleaned_smiles = clean_all(frag_smiles)
        #print(cleaned_smiles)
        dialog = EditBB(bb_name, cleaned_smiles, is_missing_smiles, self)
        
        if dialog.exec():
            data = dialog.get_data()
            
            if data is None: 
                return
            
            if not data["name"]:
                QMessageBox.warning(self, "Invalid Input", "Molecule name cannot be empty.")
                return
            
            if is_missing_smiles and not data["smiles"]:
                QMessageBox.warning(self, "Invalid Input", "SMILES cannot be empty for a new Molecule.")
                return
            
            if not data["connectors"]:
                QMessageBox.warning(self, "Invalid Input", "At least one connector must be specified.")
                return
            
            molecule_name = data["name"]
            smiles = clean_all(data["smiles"]) 
            bb_variant = data["bb"]
            connectors = data["connectors"]
            
            if molecule_name not in self.database:
                self.database[molecule_name] = {"smiles": smiles, "BB": {}}
            else:
                if is_missing_smiles:
                    self.database[molecule_name]["smiles"] = smiles
            
            self.database[molecule_name]["BB"][bb_variant] = connectors
            
            self.save_database()
            
            self.result_table.setItem(row, 0, QTableWidgetItem(molecule_name))
            self.result_table.setItem(row, 1, QTableWidgetItem(smiles))
            self.result_table.setItem(row, 2, QTableWidgetItem("Found"))
            self.result_table.removeCellWidget(row, 3)

        self.cut_and_search()

    def generate_script(self):
        molecule_name = self.name_input.text()
        smiles = self.smiles_input.toPlainText()

        fragments = fragment_molecule(smiles)
        if fragments:
            results = search_database(fragments)
            self.show_results(results)
        else:
            QMessageBox.warning(self, "Error", "Invalid SMILES or unable to fragment.")
            return

        bb_names = []
        for row in range(self.result_table.rowCount()):
            bb_name = self.result_table.item(row, 0).text()
            status = self.result_table.item(row, 2).text()
            if status == "Found":
                bb_names.append(bb_name)
            else:
                bb_names.append("Not Found")

        # Check itp & rtp files
        missing_files = []
        for bb in bb_names:
            if bb != "Not Found":
                if not os.path.exists(f"qforce_itp/{bb}.itp"):
                    missing_files.append(f"{bb}.itp")
                if not os.path.exists(f"qforce_rtp/{bb}.rtp"):
                    missing_files.append(f"{bb}.rtp")
        
        if missing_files:
            report = "The following files are missing:\n" + "\n".join(missing_files)
            report += "\nPlease add them to the appropriate folders before generating the script."
            self.script_output.setText(report)
            self.confirm_button.setEnabled(False)
            return

        # Generate script
        BB_connector = ["Br" if bb != "Not Found" else "Not Found" for bb in bb_names]

        script_text = "All required files ready, generating:\n\n"
        script_text += f'Molecule_Name = "{molecule_name}"\n'
        script_text += f'Building_Blocks = {bb_names}\n'
        script_text += f'Connectors = {BB_connector}\n'
        script_text += '\n#Note: We currently assume Br for all connectors. Please adjust to "I" or "F" if needed.'

        self.script_output.setText(script_text)
        self.confirm_button.setEnabled(True)


    def confirm_action(self):
        updated_script = self.script_output.toPlainText()
        script_lines = [line for line in updated_script.split('\n') if '=' in line]
        
        local_namespace = {}
        try:
            for line in script_lines:
                exec(line, globals(), local_namespace)
            
            molecule_name = local_namespace.get('Molecule_Name')
            bb_names = local_namespace.get('Building_Blocks')
            BB_connector_str = local_namespace.get('Connectors')

            if not all([molecule_name, bb_names, BB_connector_str]):
                raise ValueError("Missing required variables in the script")

            BB_connector = []
            for connector in BB_connector_str:
                if connector == "Br":
                    BB_connector.append(stk.BromoFactory())
                elif connector == "I":
                    BB_connector.append(stk.IodoFactory())
                elif connector == "F":
                    BB_connector.append(stk.FluoroFactory())
                else:
                    raise ValueError(f"Unsupported connector type: {connector}")

            FF = link(molecule_name, bb_names, BB_connector)
            self.visualization_window = FF.show_molecule()

            json_file_path = f"molecule/{molecule_name}/STRUCTURE/{molecule_name}.json"
            try:
                with open(json_file_path, 'r') as json_file:
                    data = json.load(json_file)
                    generated_smiles = data[molecule_name]['smiles']
            except FileNotFoundError:
                self.match_output.setText(f"Error: {json_file_path} not found.")
                return
            except KeyError:
                self.match_output.setText(f"Error: SMILES not found in {json_file_path}.")
                return

            input_smiles = self.smiles_input.toPlainText()
            input_mol = Chem.MolFromSmiles(input_smiles)
            generated_mol = Chem.MolFromSmiles(generated_smiles)

            if input_mol is None or generated_mol is None:
                self.match_output.setText("Error: Invalid SMILES string.")
                return

            input_inchi = Chem.MolToInchi(input_mol)
            generated_inchi = Chem.MolToInchi(generated_mol)

            if input_inchi == generated_inchi:
                self.match_output.setText("The structures match!")
            else:
                self.match_output.setText("The structures do not match.\nYou will need to adjust the connectors yourself.")

        except Exception as e:
            self.match_output.setText(f"Error: {str(e)}")

def main():
    app = QApplication(sys.argv)
    ex = AutoFF()
    ex.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()