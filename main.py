import sys
import stk
import json
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, QLabel, QLineEdit, QTextEdit
from PyQt6.QtCore import QSize
from cut_search import fragment_molecule, search_database
from link import link
import py3Dmol
from rdkit import Chem

class MoleculeApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

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
        self.smiles_input.setFixedHeight(60)
        layout.addWidget(self.smiles_input)

        self.cut_button = QPushButton('Cut and Search', self)
        self.cut_button.clicked.connect(self.cut_and_search)
        layout.addWidget(self.cut_button)

        self.result_table = QTableWidget(self)
        self.result_table.setFixedHeight(300)  
        layout.addWidget(self.result_table)

        self.script_button = QPushButton('Generate Script', self)
        self.script_button.clicked.connect(self.generate_script)
        self.script_button.setEnabled(False)  
        layout.addWidget(self.script_button)

        # Textbox for displaying generated script
        self.script_output = QTextEdit(self)
        self.script_output.setFixedHeight(100)
        self.script_output.setReadOnly(False)
        layout.addWidget(self.script_output)

        self.confirm_button = QPushButton('Confirm Info', self)
        self.confirm_button.clicked.connect(self.confirm_action)
        self.confirm_button.setEnabled(False) 
        layout.addWidget(self.confirm_button)

        # Show if the smiles matches
        self.match_output = QTextEdit(self)
        self.script_output.setFixedHeight(100)
        layout.addWidget(self.match_output)

        self.setLayout(layout)

    def cut_and_search(self):
        molecule_name = self.name_input.text()
        smiles = self.smiles_input.toPlainText()

        if smiles and molecule_name:
            # Step 1: Fragment the molecule
            fragments = fragment_molecule(smiles)
            if fragments:
                # Step 2: Search fragments in the database
                results = search_database(fragments)
                self.show_results(results)
                self.script_button.setEnabled(True)
            else:
                self.result_table.setRowCount(0)
                self.result_table.setColumnCount(1)
                self.result_table.setHorizontalHeaderLabels(["Error"])
                self.result_table.setItem(0, 0, QTableWidgetItem("Invalid SMILES or unable to fragment."))

    def show_results(self, results):
        self.result_table.setRowCount(len(results))
        self.result_table.setColumnCount(2)
        self.result_table.setHorizontalHeaderLabels(["BB Name", "Fragment SMILES"])
        self.result_table.setColumnWidth(0, 200)  
        self.result_table.setColumnWidth(1, 500) 

        for i, (bb_name, frag_smiles, _) in enumerate(results):
            self.result_table.setItem(i, 0, QTableWidgetItem(bb_name if bb_name else "Not Found"))
            self.result_table.setItem(i, 1, QTableWidgetItem(frag_smiles))

    def generate_script(self):
        molecule_name = self.name_input.text()
        bb_names = [self.result_table.item(row, 0).text() for row in range(self.result_table.rowCount())]

        # 假設所有接點都是 Br
        Br = "Br"
        BB_connector = [Br] * len(bb_names)

        # Generate script
        script_text = f"Molecule Name = \"{molecule_name}\"\n"
        script_text += f"LEGO/Building Blocks = {bb_names}\n"
        script_text += f"Connectors = {BB_connector}\n"
        script_text += f"\n #Note: We currently assume Br for all connectors. Please adjust to 'I' or 'F' if needed."

        self.script_output.setText(script_text)
        self.confirm_button.setEnabled(True)

    def confirm_action(self):
        updated_script = self.script_output.toPlainText()
        local_namespace = {}
        
        try:
            exec(updated_script, globals(), local_namespace)
            
            # Extract variables from the local namespace
            molecule_name = local_namespace.get('Molecule Name')
            bb_names = local_namespace.get('LEGO/Building Blocks')
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

            # Read the generated SMILES from the JSON file
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

            # Compare input SMILES with the generated one
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
                self.match_output.setText("The structures do not match.\n You will need to adjust based on the alert below.")

        except Exception as e:
            self.match_output.setText(f"Error: {str(e)}")
    

def main():
    app = QApplication(sys.argv)
    ex = MoleculeApp()
    ex.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
