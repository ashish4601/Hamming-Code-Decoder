# Hamming Code Decoder with ASCII Message Recovery and PGM Row Interchange

This project implements a **Hamming Code Decoder** that:
- Corrects **single-bit errors** in standard Hamming codewords.
- Recovers **ASCII messages** encoded with two (7,4) Hamming blocks per character.
- Optionally supports **row interchange** in grayscale **PGM image** files.

---

## 🔧 Features

### ✅ Hamming Decoder
- Accepts a list of binary strings (codewords).
- Uses `r` parity bits to detect and fix one-bit errors.
- Returns the original message bits after correction.

### 🔡 ASCII Message Decoder (1-bit error per 7 bits)
- Each 8-bit ASCII character is split into two 4-bit segments.
- Each 4-bit segment is Hamming (7,4) encoded and may contain at most one error.
- The decoder corrects errors and reconstructs the ASCII string.

### 🖼️ PGM Row Interchange (Optional)
- Supports loading a `.pgm` (grayscale) image file.
- Allows users to interchange rows in the image.
- Saves the modified image back in PGM format.

---

## 📥 Inputs

1. **r**: Number of parity bits (usually 3 for Hamming(7,4)).
2. **codewords**: List of received Hamming-encoded binary strings.
3. **ascii_codewords**: List of 7-bit strings (each encoded with 1-bit error).
4. **pgm_image (optional)**: Path to a PGM file for row interchange functionality.

---

## 📤 Outputs

- Decoded message bits from Hamming code.
- Corrected ASCII string from 7-bit encoded blocks.
- Modified PGM file (if row swapping is used).

---

## 🚀 How to Run

1. Clone the repo:
```bash
git clone https://github.com/yourusername/hamming-code-decoder.git
cd hamming-code-decoder
