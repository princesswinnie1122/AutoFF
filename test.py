from PyQt6.QtWidgets import QWidget, QVBoxLayout, QApplication, QPushButton
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWebEngineCore import QWebEngineSettings
from PyQt6.QtCore import QUrl
import sys

class MoleculeViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('3Dmol.js Direct Molecule Viewer')
        self.setGeometry(100, 100, 800, 600)

        layout = QVBoxLayout()

        self.html_viewer = QWebEngineView(self)
        self.html_viewer.settings().setAttribute(QWebEngineSettings.WebAttribute.WebGLEnabled, True)
        layout.addWidget(self.html_viewer)

        debug_button = QPushButton('Print Debug Info', self)
        debug_button.clicked.connect(self.print_debug_info)
        layout.addWidget(debug_button)

        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <style>
                body { font-family: Arial, sans-serif; background-color: #f0f0f0; }
                h1 { color: blue; }
                #viewer { border: 1px solid black; width: 400px; height: 400px; position: relative; }
                #status { border: 1px solid #ccc; padding: 10px; margin-top: 20px; }
            </style>
        </head>
        <body>
            <h1>3Dmol.js Direct Molecule Viewer</h1>
            <div id="viewer"></div>
            <div id="status"></div>
            <script>
                function log(message) {
                    document.getElementById('status').innerHTML += '<p>' + message + '</p>';
                    console.log(message);
                }

                log('Script started');

                var viewer;
                var moleculeAdded = false;

                function initViewer() {
                    try {
                        log('Creating 3Dmol.js viewer');
                        viewer = $3Dmol.createViewer("viewer", {backgroundColor: 'lightgray'});
                        log('Viewer created');

                        // Create a simple molecule (water)
                        var atoms = [
                            {elem: 'O', x: 0, y: 0, z: 0},
                            {elem: 'H', x: 0.8, y: 0.6, z: 0},
                            {elem: 'H', x: -0.8, y: 0.6, z: 0}
                        ];
                        
                        viewer.addModel(atoms, "json");
                        viewer.setStyle({}, {sphere:{}, stick:{}});
                        viewer.zoomTo();
                        viewer.render();
                        log('Water molecule created and displayed');
                        moleculeAdded = true;
                    } catch (error) {
                        log('Error: ' + error.message);
                    }
                }

                $(document).ready(function() {
                    log('Document ready');
                    initViewer();
                });

                function getDebugInfo() {
                    return {
                        windowInnerSize: {width: window.innerWidth, height: window.innerHeight},
                        documentSize: {width: document.documentElement.clientWidth, height: document.documentElement.clientHeight},
                        viewerSize: {width: document.getElementById('viewer').clientWidth, height: document.getElementById('viewer').clientHeight},
                        userAgent: navigator.userAgent,
                        is3DmolDefined: (typeof $3Dmol !== 'undefined'),
                        moleculeAdded: moleculeAdded,
                        statusContent: document.getElementById('status').innerHTML
                    };
                }
            </script>
        </body>
        </html>
        """

        self.html_viewer.setHtml(html_content, QUrl("https://3dmol.csb.pitt.edu/"))
        self.html_viewer.loadFinished.connect(self.onLoadFinished)

        self.setLayout(layout)

    def onLoadFinished(self, ok):
        if ok:
            print("Page loaded successfully")
            self.html_viewer.page().runJavaScript(
                "document.getElementById('status').innerText",
                self.printResult
            )
        else:
            print("Page failed to load")

    def printResult(self, result):
        print("Status from webpage:", result)

    def print_debug_info(self):
        self.html_viewer.page().runJavaScript("getDebugInfo()", self.display_debug_info)

    def display_debug_info(self, info):
        print("Debug Info:")
        print(info)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    viewer = MoleculeViewer()
    viewer.show()
    sys.exit(app.exec())