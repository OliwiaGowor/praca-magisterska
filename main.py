 
import numpy as np
import cv2
import time
import math
import os
from scipy.stats import entropy as scipy_entropy

# ===============================================
# Helper functions for security and performance metrics
# ===============================================

def calculate_entropy(img):
    """Calculate image entropy"""
    hist, _ = np.histogram(img.flatten(), bins=256, range=(0, 256))
    hist = hist / np.sum(hist)
    return scipy_entropy(hist, base=2)

def calculate_correlation(img1, img2):
    """Calculate correlation coefficient between two images"""
    img1 = img1.flatten()
    img2 = img2.flatten()
    return np.corrcoef(img1, img2)[0, 1]

def calculate_NPCR(img1, img2):
    """Calculate NPCR (Number of Pixels Change Rate)"""
    diff = img1 != img2
    return np.sum(diff) / diff.size * 100

def calculate_noise(img1, img2):
    """Simple noise measurement (mean pixel difference)"""
    return np.mean(np.abs(img1.astype(np.int16) - img2.astype(np.int16)))

# ===============================================
# Encryption algorithm placeholders
# ===============================================

def classical_encryption(img):
    """TODO: Implement classical encryption (e.g., AES, DES, XOR, etc.)"""
    # Placeholder – simulate processing time
    time.sleep(0.5)
    return np.flipud(img)  # temporary modification

def chaotic_encryption(img):
    """TODO: Implement chaotic encryption (e.g., Logistic Map, Chen, Lorenz, etc.)"""
    # Placeholder – simulate processing time
    time.sleep(0.7)
    return np.fliplr(img)  # temporary modification

# ===============================================
# Main comparison logic
# ===============================================

def compare_methods(image_path, results_file="results.txt"):
    if not os.path.exists(image_path):
        print(f"Error: file not found - {image_path}")
        return

    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        print("Error: failed to read image.")
        return

    results = []

    # ---- Classical encryption ----
    start = time.time()
    encrypted_classical = classical_encryption(img)
    classical_time = time.time() - start

    results.append({
        "method": "classical",
        "entropy": calculate_entropy(encrypted_classical),
        "correlation": calculate_correlation(img, encrypted_classical),
        "NPCR": calculate_NPCR(img, encrypted_classical),
        "noise": calculate_noise(img, encrypted_classical),
        "time": classical_time,
        "performance": img.size / classical_time
    })

    # ---- Chaotic encryption ----
    start = time.time()
    encrypted_chaotic = chaotic_encryption(img)
    chaotic_time = time.time() - start

    results.append({
        "method": "chaotic",
        "entropy": calculate_entropy(encrypted_chaotic),
        "correlation": calculate_correlation(img, encrypted_chaotic),
        "NPCR": calculate_NPCR(img, encrypted_chaotic),
        "noise": calculate_noise(img, encrypted_chaotic),
        "time": chaotic_time,
        "performance": img.size / chaotic_time
    })

    # ===============================================
    # Save results to file
    # ===============================================
    with open(results_file, "w", encoding="utf-8") as f:
        f.write("Comparison of Classical and Chaotic Image Encryption Methods\n")
        f.write("=" * 70 + "\n\n")

        for r in results:
            f.write(f"Method: {r['method']}\n")
            f.write(f"  Entropy: {r['entropy']:.4f}\n")
            f.write(f"  Correlation: {r['correlation']:.4f}\n")
            f.write(f"  NPCR: {r['NPCR']:.2f}%\n")
            f.write(f"  Noise: {r['noise']:.4f}\n")
            f.write(f"  Time: {r['time']:.4f} s\n")
            f.write(f"  Performance: {r['performance']:.2f} pixels/s\n")
            f.write("-" * 70 + "\n")

    print(f"Results saved to file: {results_file}")

# ===============================================
# Console entry point
# ===============================================
if __name__ == "__main__":
    path = input("Enter image path (e.g., lena.png): ").strip()
    compare_methods(path)
