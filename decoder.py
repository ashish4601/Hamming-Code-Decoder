import sys
import os
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets

# ========= GF(2) Helper Functions =========

def gf2_row_reduce(A):
    A = A.copy() % 2
    m, n = A.shape
    pivot_cols = []
    row = 0
    for col in range(n):
        pivot = None
        for r in range(row, m):
            if A[r, col] == 1:
                pivot = r
                break
        if pivot is None:
            continue
        A[[row, pivot]] = A[[pivot, row]]
        pivot_cols.append(col)
        for r in range(m):
            if r != row and A[r, col] == 1:
                A[r] = (A[r] + A[row]) % 2
        row += 1
        if row == m:
            break
    return A, pivot_cols

def gf2_nullspace(G):
    k, n = G.shape
    R, pivot_cols = gf2_row_reduce(G)
    pivot_cols = set(pivot_cols)
    free_cols = [j for j in range(n) if j not in pivot_cols]
    basis = []
    for free in free_cols:
        vec = np.zeros(n, dtype=int)
        vec[free] = 1
        for i, col in enumerate(sorted(pivot_cols)):
            vec[col] = (-R[i, free]) % 2
        basis.append(vec)
    if len(basis) == 0:
        return np.zeros((0, n), dtype=int)
    return np.array(basis, dtype=int)

def compute_pcm_from_generator(G):
    """
    For a generator matrix G in systematic form [I | P],
    compute the parity-check matrix (PCM) of the form:
       H = ([P^T | I])^T
    """
    k, n = G.shape
    r = n - k
    # Extract P from G = [I | P]
    P = G[:, k:]
    # Build [P^T | I] (dimensions: r x (k+r)) then take its transpose to get a (k+r) x r matrix.
    H = np.concatenate((P.T % 2, np.eye(r, dtype=int)), axis=1).T
    return H

# ========= Default Generator Matrix Generation =========

def generate_default_pgm(r):
    n = 2**r - 1
    k = n - r
    all_vecs = [np.array(list(map(int, format(i, '0{}b'.format(r)))), dtype=int)
                for i in range(1, 2**r)]
    basis = [np.eye(1, r, i, dtype=int).flatten() for i in range(r)]
    A_rows = [vec for vec in all_vecs if not any(np.array_equal(vec, b) for b in basis)]
    if len(A_rows) != k:
        raise ValueError("Dimension mismatch while constructing default A.")
    A = np.array(A_rows, dtype=int)
    I_k = np.eye(k, dtype=int)
    G = np.concatenate((I_k, A), axis=1)
    return G

# ========= PyQt5 GUI Classes =========

class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Hamming Decoder Main Menu")
        self.resize(400, 200)
        layout = QtWidgets.QVBoxLayout(self)

        # Input for r
        r_layout = QtWidgets.QHBoxLayout()
        r_label = QtWidgets.QLabel("Enter number of parity bits (r):")
        self.r_edit = QtWidgets.QLineEdit()
        r_layout.addWidget(r_label)
        r_layout.addWidget(self.r_edit)
        layout.addLayout(r_layout)

        self.compute_button = QtWidgets.QPushButton("Compute Parameters")
        self.compute_button.clicked.connect(self.compute_parameters)
        layout.addWidget(self.compute_button)

        self.param_label = QtWidgets.QLabel("")
        layout.addWidget(self.param_label)

        self.general_decoder_btn = QtWidgets.QPushButton("General Hamming Decoder")
        self.general_decoder_btn.clicked.connect(self.open_general_decoder)
        self.general_decoder_btn.setEnabled(False)

        self.ascii_decoder_btn = QtWidgets.QPushButton("ASCII Hamming Decoder (7,4)")
        self.ascii_decoder_btn.clicked.connect(self.open_ascii_decoder)
        layout.addWidget(self.general_decoder_btn)
        layout.addWidget(self.ascii_decoder_btn)

    def compute_parameters(self):
        try:
            r = int(self.r_edit.text().strip())
            if r < 2:
                raise ValueError("r must be at least 2.")
            n = 2**r - 1
            k = n - r
            self.param_label.setText(f"Computed: n = {n}, k = {k}")
            self.general_decoder_btn.setEnabled(True)
            self.r_value = r
            self.n_value = n
            self.k_value = k
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", str(e))

    def open_general_decoder(self):
        self.genDecoderWin = GeneralHammingDecoderWindow(self.r_value, self.n_value, self.k_value)
        self.genDecoderWin.show()

    def open_ascii_decoder(self):
        self.asciiDecoderWin = ASCIIHammingDecoderWindow()
        self.asciiDecoderWin.show()

class GeneralHammingDecoderWindow(QtWidgets.QWidget):
    def __init__(self, r, n, k):
        super().__init__()
        self.setWindowTitle("General Hamming Decoder (r-based)")
        self.resize(900, 700)
        self.r = r
        self.n = n
        self.k = k

        # Create default generator matrix (PGM) in systematic form [I | P] and its corresponding PCM.
        self.G = generate_default_pgm(self.r)  # shape: k x n
        self.PCM = compute_pcm_from_generator(self.G)  # shape: n x (n-k)

        layout = QtWidgets.QVBoxLayout(self)

        param_group = QtWidgets.QGroupBox("Parameters")
        param_layout = QtWidgets.QHBoxLayout()
        param_layout.addWidget(QtWidgets.QLabel(f"r = {self.r}, n = {self.n}, k = {self.k}"))
        param_group.setLayout(param_layout)
        layout.addWidget(param_group)

        matrices_group = QtWidgets.QGroupBox("Generator Matrix (PGM) and Parity-check Matrix (PCM)")
        matrices_layout = QtWidgets.QHBoxLayout()

        # PGM display using QPlainTextEdit.
        pgm_group = QtWidgets.QGroupBox("PGM")
        pgm_layout = QtWidgets.QVBoxLayout()
        self.pgm_display = QtWidgets.QPlainTextEdit()
        self.pgm_display.setReadOnly(True)
        self.pgm_display.setFont(QtGui.QFont("Courier", 10))
        pgm_layout.addWidget(self.pgm_display)
        # Controls for modifying PGM rows.
        shift_layout = QtWidgets.QHBoxLayout()
        shift_layout.addWidget(QtWidgets.QLabel("Row index (0-indexed):"))
        self.row_index_edit = QtWidgets.QLineEdit()
        self.row_index_edit.setFixedWidth(50)
        shift_layout.addWidget(self.row_index_edit)
        self.shift_up_btn = QtWidgets.QPushButton("Shift Row Up")
        self.shift_up_btn.clicked.connect(self.shift_row_up)
        shift_layout.addWidget(self.shift_up_btn)
        self.shift_down_btn = QtWidgets.QPushButton("Shift Row Down")
        self.shift_down_btn.clicked.connect(self.shift_row_down)
        shift_layout.addWidget(self.shift_down_btn)
        # Refresh PCM button.
        self.refresh_pcm_btn = QtWidgets.QPushButton("Refresh PCM")
        self.refresh_pcm_btn.clicked.connect(self.update_pcm_display)
        shift_layout.addWidget(self.refresh_pcm_btn)
        pgm_layout.addLayout(shift_layout)
        pgm_group.setLayout(pgm_layout)
        matrices_layout.addWidget(pgm_group)

        # PCM display using QPlainTextEdit.
        pcm_group = QtWidgets.QGroupBox("PCM")
        pcm_layout = QtWidgets.QVBoxLayout()
        self.pcm_display = QtWidgets.QPlainTextEdit()
        self.pcm_display.setReadOnly(True)
        self.pcm_display.setFont(QtGui.QFont("Courier", 10))
        pcm_layout.addWidget(self.pcm_display)
        pcm_group.setLayout(pcm_layout)
        matrices_layout.addWidget(pcm_group)

        matrices_group.setLayout(matrices_layout)
        layout.addWidget(matrices_group)

        # Codeword Input and Decoding section.
        codeword_group = QtWidgets.QGroupBox("Codeword Input and Decoding")
        cw_layout = QtWidgets.QVBoxLayout()
        file_btn = QtWidgets.QPushButton("Load Encoded Codeword from File (.txt)")
        file_btn.clicked.connect(self.load_codeword_file)
        cw_layout.addWidget(file_btn)
        cw_layout.addWidget(QtWidgets.QLabel("Or paste encoded binary text (only 0s and 1s):"))
        self.codeword_edit = QtWidgets.QTextEdit()
        cw_layout.addWidget(self.codeword_edit)
        self.decode_btn = QtWidgets.QPushButton("Decode Codeword")
        self.decode_btn.clicked.connect(self.decode_codeword)
        cw_layout.addWidget(self.decode_btn)
        cw_layout.addWidget(QtWidgets.QLabel("Decoded Bit String (concatenated message bits):"))
        self.decoded_out = QtWidgets.QTextEdit()
        self.decoded_out.setReadOnly(True)
        cw_layout.addWidget(self.decoded_out)
        cw_layout.addWidget(QtWidgets.QLabel("Syndromes (one per codeword block):"))
        self.syndrome_out = QtWidgets.QTextEdit()
        self.syndrome_out.setReadOnly(True)
        cw_layout.addWidget(self.syndrome_out)
        codeword_group.setLayout(cw_layout)
        layout.addWidget(codeword_group)

        self.update_pgm_display()
        self.update_pcm_display()

    def matrix_to_text(self, M):
        lines = []
        for row in M:
            lines.append(" ".join(str(int(x)) for x in row))
        return "\n".join(lines)

    def update_pgm_display(self):
        text = self.matrix_to_text(self.G)
        self.pgm_display.setPlainText(text)

    def update_pcm_display(self):
        # Simply update the PCM display using the current self.PCM.
        text = self.matrix_to_text(self.PCM)
        self.pcm_display.setPlainText(text)

    def shift_row_up(self):
        try:
            idx = int(self.row_index_edit.text().strip())
            if idx <= 0 or idx >= self.G.shape[0]:
                QtWidgets.QMessageBox.critical(self, "Error", "Invalid row index for shift up.")
                return
            # Permute rows in PGM.
            self.G[[idx-1, idx]] = self.G[[idx, idx-1]]
            # Permute the corresponding rows in the top block of PCM.
            pcm_top = self.PCM[:self.k, :].copy()
            pcm_top[[idx-1, idx]] = pcm_top[[idx, idx-1]]
            self.PCM[:self.k, :] = pcm_top
            self.update_pgm_display()
            self.update_pcm_display()
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", str(e))

    def shift_row_down(self):
        try:
            idx = int(self.row_index_edit.text().strip())
            if idx < 0 or idx >= self.G.shape[0]-1:
                QtWidgets.QMessageBox.critical(self, "Error", "Invalid row index for shift down.")
                return
            # Permute rows in PGM.
            self.G[[idx, idx+1]] = self.G[[idx+1, idx]]
            # Permute the corresponding rows in the top block of PCM.
            pcm_top = self.PCM[:self.k, :].copy()
            pcm_top[[idx, idx+1]] = pcm_top[[idx+1, idx]]
            self.PCM[:self.k, :] = pcm_top
            self.update_pgm_display()
            self.update_pcm_display()
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Error", str(e))

    def load_codeword_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open Codeword File", "", "Text Files (*.txt)")
        if path:
            with open(path, 'r') as f:
                data = f.read().strip()
            self.codeword_edit.setPlainText(data)

    def decode_codeword(self):
        cw_text = self.codeword_edit.toPlainText().strip().replace(" ", "").replace("\n", "")
        if set(cw_text) - set("01"):
            QtWidgets.QMessageBox.critical(self, "Error", "Codeword must contain only 0s and 1s.")
            return
        if len(cw_text) % self.n != 0:
            QtWidgets.QMessageBox.critical(self, "Error", f"Length of codeword must be a multiple of n = {self.n}.")
            return
        
        decoded_bits = ""
        syndrome_list = []
        num_blocks = len(cw_text) // self.n
        # H is computed from PCM; note that here we use PCM.T to form H.
        H = self.PCM.T  # H has shape (n-k) x n.
        for b in range(num_blocks):
            block_str = cw_text[b*self.n:(b+1)*self.n]
            block = np.array(list(map(int, list(block_str))), dtype=int)
            syndrome = (H @ block) % 2
            syndrome_str = ''.join(map(str, syndrome.tolist()))
            syndrome_list.append(syndrome_str)
            if np.any(syndrome):
                error_index = None
                for j in range(self.n):
                    if np.array_equal(H[:, j], syndrome):
                        error_index = j
                        break
                if error_index is not None:
                    block[error_index] ^= 1  # correct error by flipping bit
            decoded_bits += ''.join(map(str, block[:self.k].tolist()))
        self.decoded_out.setPlainText(decoded_bits)
        self.syndrome_out.setPlainText("\n".join(syndrome_list))

class ASCIIHammingDecoderWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ASCII Hamming (7,4) Decoder with Systematic Format")
        self.resize(700, 600)
        layout = QtWidgets.QVBoxLayout(self)

        formula_group = QtWidgets.QGroupBox("Parity Bit XOR Formulas")
        formula_layout = QtWidgets.QVBoxLayout()
        formula_text = ("For each 4-bit nibble (m₁ m₂ m₃ m₄):\n"
                        "p₁ = m₁ ⊕ m₂ ⊕ m₄\n"
                        "p₂ = m₁ ⊕ m₃ ⊕ m₄\n"
                        "p₃ = m₂ ⊕ m₃ ⊕ m₄\n"
                        "Encoded codeword format: [m₁, m₂, m₃, m₄, p₁, p₂, p₃]")
        formula_layout.addWidget(QtWidgets.QLabel(formula_text))
        formula_group.setLayout(formula_layout)
        layout.addWidget(formula_group)

        sample_group = QtWidgets.QGroupBox("Sample Codeword for 'Hello'")
        sample_layout = QtWidgets.QVBoxLayout()
        sample_codeword = self.compute_sample_codeword("Hello")
        sample_layout.addWidget(QtWidgets.QLabel(sample_codeword))
        sample_group.setLayout(sample_layout)
        layout.addWidget(sample_group)

        input_group = QtWidgets.QGroupBox("Enter/Load Encoded Codeword (length must be a multiple of 7)")
        input_layout = QtWidgets.QVBoxLayout()
        file_btn = QtWidgets.QPushButton("Load Encoded Codeword from File (.txt)")
        file_btn.clicked.connect(self.load_codeword_file)
        input_layout.addWidget(file_btn)
        input_layout.addWidget(QtWidgets.QLabel("Or paste encoded binary text:"))
        self.codeword_edit = QtWidgets.QTextEdit()
        input_layout.addWidget(self.codeword_edit)
        self.decode_btn = QtWidgets.QPushButton("DECODE MESSAGE")
        self.decode_btn.clicked.connect(self.decode_message)
        input_layout.addWidget(self.decode_btn)
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)

        output_group = QtWidgets.QGroupBox("Decoded Message")
        output_layout = QtWidgets.QVBoxLayout()
        self.decoded_out = QtWidgets.QTextEdit()
        self.decoded_out.setReadOnly(True)
        output_layout.addWidget(self.decoded_out)
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)

    def load_codeword_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open Codeword File", "", "Text Files (*.txt)")
        if path:
            with open(path, 'r') as f:
                data = f.read().strip()
            self.codeword_edit.setPlainText(data)

    def encode_hamming74(self, nibble):
        if len(nibble) != 4 or any(c not in "01" for c in nibble):
            raise ValueError("Nibble must be a 4-bit binary string.")
        m1, m2, m3, m4 = map(int, nibble)
        p1 = m1 ^ m2 ^ m4
        p2 = m1 ^ m3 ^ m4
        p3 = m2 ^ m3 ^ m4
        return f"{m1}{m2}{m3}{m4}{p1}{p2}{p3}"

    def compute_sample_codeword(self, message):
        sample_code = ""
        for ch in message:
            char_bits = format(ord(ch), '08b')
            nibble1 = char_bits[:4]
            nibble2 = char_bits[4:]
            code1 = self.encode_hamming74(nibble1)
            code2 = self.encode_hamming74(nibble2)
            sample_code += code1 + code2
        return sample_code

    def decode_hamming74(self, codeword):
        if len(codeword) != 7 or any(c not in "01" for c in codeword):
            raise ValueError("Each codeword must be a 7-bit binary string.")
        bits = np.array(list(map(int, list(codeword))), dtype=int)
        Q = np.array([[1, 1, 0],
                      [1, 0, 1],
                      [0, 1, 1],
                      [1, 1, 1]], dtype=int)
        H = np.concatenate((Q.T, np.eye(3, dtype=int)), axis=1)
        syndrome = (H @ bits) % 2
        if np.any(syndrome):
            for j in range(7):
                if np.array_equal(H[:, j], syndrome):
                    bits[j] ^= 1
                    break
        return ''.join(map(str, bits[:4]))

    def decode_message(self):
        cw_text = self.codeword_edit.toPlainText().strip().replace(" ", "").replace("\n", "")
        if set(cw_text) - set("01"):
            QtWidgets.QMessageBox.critical(self, "Error", "Codeword must contain only 0s and 1s.")
            return
        if len(cw_text) % 7 != 0:
            QtWidgets.QMessageBox.critical(self, "Error", "Length of codeword must be a multiple of 7.")
            return
        num_blocks = len(cw_text) // 7
        nibbles = []
        for b in range(num_blocks):
            block = cw_text[b*7:(b+1)*7]
            nibble = self.decode_hamming74(block)
            nibbles.append(nibble)
        if len(nibbles) % 2 != 0:
            QtWidgets.QMessageBox.critical(self, "Error", "Decoded nibble count is odd; cannot form complete ASCII characters.")
            return
        message = ""
        for i in range(0, len(nibbles), 2):
            byte = nibbles[i] + nibbles[i+1]
            message += chr(int(byte, 2))
        self.decoded_out.setPlainText(message)

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main_win = MainWindow()
    main_win.show()
    sys.exit(app.exec_())
