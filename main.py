import sys
import stk
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, QLabel, QLineEdit, QTextEdit
from PyQt6.QtWidgets import QHBoxLayout
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
        # Set a longer window size for better display of fragmented SMILES
        self.resize(QSize(800, 600))  # Wider and taller window
        self.setWindowTitle('Molecule Builder')

        layout = QVBoxLayout()

        # Input for Molecule Name
        self.name_label = QLabel("Enter the name of your molecule:")
        layout.addWidget(self.name_label)

        # Larger input box for Molecule Name
        self.name_input = QLineEdit(self)
        self.name_input.setFixedHeight(30)
        layout.addWidget(self.name_input)

        # Input for SMILES
        self.smiles_label = QLabel("Enter the canonical SMILES of your molecule:")
        layout.addWidget(self.smiles_label)

        # Larger input box for SMILES with multiple lines
        self.smiles_input = QTextEdit(self)
        self.smiles_input.setFixedHeight(60)
        layout.addWidget(self.smiles_input)

        # Button to cut molecule and search database
        self.cut_button = QPushButton('Cut and Search', self)
        self.cut_button.clicked.connect(self.cut_and_search)
        layout.addWidget(self.cut_button)

        # Table for results, set larger for fragmented SMILES display
        self.result_table = QTableWidget(self)
        self.result_table.setFixedHeight(300)  # Larger height for displaying results
        layout.addWidget(self.result_table)

        # Button to generate the script
        self.script_button = QPushButton('Generate Script', self)
        self.script_button.clicked.connect(self.generate_script)
        self.script_button.setEnabled(False)  # Enable after showing results
        layout.addWidget(self.script_button)

        # Textbox for displaying the generated script
        self.script_output = QTextEdit(self)
        self.script_output.setFixedHeight(100)
        layout.addWidget(self.script_output)

        # Confirm button to generate files and visualize molecule
        self.confirm_button = QPushButton('Confirm', self)
        self.confirm_button.clicked.connect(self.confirm_action)
        self.confirm_button.setEnabled(False)  # Only enabled after generating script
        layout.addWidget(self.confirm_button)

        # Molecule visualization widget (py3Dmol integration)
        #self.molecule_viewer = QTextEdit(self)
        #self.molecule_viewer.setFixedHeight(200)  # Give space for the molecule visualization
        #layout.addWidget(self.molecule_viewer)

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
                self.script_button.setEnabled(True)  # Enable the script button
            else:
                self.result_table.setRowCount(0)
                self.result_table.setColumnCount(1)
                self.result_table.setHorizontalHeaderLabels(["Error"])
                self.result_table.setItem(0, 0, QTableWidgetItem("Invalid SMILES or unable to fragment."))

    def show_results(self, results):
        self.result_table.setRowCount(len(results))
        self.result_table.setColumnCount(2)
        self.result_table.setHorizontalHeaderLabels(["BB Name", "Fragment SMILES"])
        self.result_table.setColumnWidth(0, 200)  # Set wider column for BB Name
        self.result_table.setColumnWidth(1, 500)  # Set wider column for Fragment SMILES

        for i, (bb_name, frag_smiles, _) in enumerate(results):
            self.result_table.setItem(i, 0, QTableWidgetItem(bb_name if bb_name else "Not Found"))
            self.result_table.setItem(i, 1, QTableWidgetItem(frag_smiles))

    def generate_script(self):
        # Get the molecule name and fragments (BB names) from the table
        molecule_name = self.name_input.text()
        bb_names = [self.result_table.item(row, 0).text() for row in range(self.result_table.rowCount())]

        # Assuming all linkers are Bromine for now
        Br = "Br"
        BB_connector = [Br] * len(bb_names)

        # Generate the script
        script_text = f"molecule_name = \"{molecule_name}\"\n"
        script_text += f"BB_name = {bb_names}\n"
        script_text += f"BB_connector = {BB_connector}\n"

        # Show the generated script in the output box
        self.script_output.setText(script_text)

        # Enable the Confirm button to generate files
        self.confirm_button.setEnabled(True)

    def confirm_action(self):
        # Get the molecule name and fragments (BB names) from the table
        molecule_name = self.name_input.text()
        bb_names = [self.result_table.item(row, 0).text() for row in range(self.result_table.rowCount())]

        # Assuming all linkers are Bromine for now
        Br = stk.BromoFactory()
        BB_connector = [Br] * len(bb_names)

        # Call the link object to generate files
        FF = link(molecule_name, bb_names, BB_connector)

        # Visualize the molecule
        self.visualize_molecule(FF)

    '''
    def visualize_molecule(self, FF):
        # Get the final molecule in a format suitable for visualization (py3Dmol)
        mol = FF.molecule
        smiles = Chem.MolToSmiles(mol)

        # Use py3Dmol to visualize the molecule
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(smiles, "smi")
        viewer.setStyle({'stick': {}})
        viewer.zoomTo()
        viewer.render()
        
        # Display the molecule in the widget
        html_view = viewer.getHTML()
        self.molecule_viewer.setHtml(html_view)
    '''
    
def main():
    app = QApplication(sys.argv)
    ex = MoleculeApp()
    ex.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
