# Analiza porównawcza algorytmów klasycznych i chaotycznych w szyfrowaniu obrazów

Repozytorium zawiera kod źródłowy oraz środowisko testowe przygotowane w ramach pracy dyplomowej. Projekt realizuje kompleksową analizę porównawczą efektywności i bezpieczeństwa klasycznych algorytmów kryptograficznych oraz nowoczesnych metod bazujących na teorii chaosu deterministycznego.

## O projekcie
Celem projektu jest zbadanie przydatności algorytmów chaotycznych w szyfrowaniu multimediów (obrazów cyfrowych). Środowisko umożliwia:
1. Szyfrowanie i deszyfrowanie obrazów różnymi metodami.
2. Przeprowadzenie kryptoanalizy statystycznej (histogramy, korelacja, entropia).
3. Analizę wrażliwości (NPCR, UACI, wykładniki Lapunowa).
4. Testy odporności na ataki destrukcyjne (szum, wycinanie danych).
5. Porównanie wydajności czasowej implementacji (MATLAB vs C/MEX).

## Zaimplementowane algorytmy

W projekcie zaimplementowano 10 algorytmów podzielonych na dwie grupy:

### 1. Algorytmy Klasyczne (Standardy)
Implementacje wariantowe (czysty MATLAB oraz zoptymalizowane C-MEX):
* **AES-CBC** (Advanced Encryption Standard)
* **DES-CBC** (Data Encryption Standard)
* **Blowfish-CBC**
* **ChaCha20**

### 2. Algorytmy Chaotyczne
Zgodnie z najnowszą literaturą przedmiotu:
* **CCS** – Circular Chaotic Shifts [Moysis et al., 2025]
* **BL-CS** – Bit-Level Cubic Shuffling [Moysis et al., 2022]
* **HC-AD** – Hyperchaotic Adaptive Diffusion [Tang, 2025]
* **3D-LC** – 3D Logistic-Chirikov System [Ramakrishna et al., 2022]
* **2S-LM** – Two-Stage Logistic Map [Hu & Tian, 2020]
* **HC-CS** – Hyperchaos & 2D Compressed Sensing [Liu et al., 2024]

## Wymagania
* **MATLAB** (R2020b lub nowszy)
* **Image Processing Toolbox**
* Komputator C/C++ (np. MinGW) – wymagany tylko do rekompilacji plików MEX (opcjonalnie).
