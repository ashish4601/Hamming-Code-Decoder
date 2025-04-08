# Hamming Code Decoder with ASCII Message Recovery and PGM Row Interchange

A complete Hamming Code Decoder application built using **PyQt5** and **NumPy**.  
This tool corrects single-bit errors in Hamming-encoded binary data, decodes ASCII messages, and allows row interchange in grayscale PGM images.

---

## 🚀 Tech Stack

- 🐍 **Python 3**
- 📦 **NumPy** – for bit-level operations and matrix handling
- 🖥️ **PyQt5** – for building the user interface (if GUI is included)

---

## 🔧 Features

### ✅ Hamming Decoder
- Accepts a list of binary strings (codewords).
- Uses `r` parity bits to detect and correct one-bit errors.
- Returns the original message bits.

### 🔡 ASCII Message Decoder (1-bit error per 7 bits)
- Each ASCII character (8 bits) is split into two 4-bit parts.
- Encoded using Hamming (7,4), with max 1 error per 7 bits.
- Corrects and reconstructs the original ASCII message.

### 🖼️ PGM Row Interchange (Optional)
- Allows custom row-swapping (e.g., for error simulation or recovery).
- Outputs a modified PGM with the new row order.

---

## 📥 Inputs

1. **r**: Number of parity bits used in encoding (typically 3 for (7,4)).
2. **codewords**: List of binary Hamming code strings.
3. **ascii_codewords**: List of 7-bit Hamming blocks (1 error per 7 bits).

---

## 📤 Outputs

- Original message bits from corrected codewords.
- Corrected ASCII message from encoded 7-bit blocks.
- Updated `.pgm` file with modified rows (if used).

---
