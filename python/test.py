import numpy as np
import hashlib
from PIL import Image
import matplotlib.pyplot as plt

# =========================================================================
# Tłumaczenie pliku: smht.m
# =========================================================================
def smht(x, a, b, c, d):
    """
    Port 1:1 funkcji smht.m
    """
    numerator = np.exp(a * x) - np.exp(-b * x)
    denominator = np.exp(c * x) + np.exp(-d * x)
    return numerator / denominator

# =========================================================================
# Tłumaczenie algorytmu: encryption_decryption_algorithm.m
# =========================================================================

def encrypt_chaotic(plaintext_image):
    """
    Tłumaczenie sekcji "Encryption Algorithm"
    """
    print("Rozpoczęcie szyfrowania chaotycznego...")

    # Używamy float64, aby odpowiadało 'double' z MATLABa
    plaintext = np.array(plaintext_image).astype(np.float64)
    rows, cols = plaintext.shape

    # Parametry (takie same jak w .m)
    DD = 10**10
    ten16 = 10**16

    # --- STEP 1 - Key initialization ---
    # Odpowiednik DataHash(plaintext,'SHA-256')
    # Musimy użyć .tobytes('C') aby dopasować kolejność bajtów MATLABa
    hash_obj = hashlib.sha256(plaintext.tobytes('C'))
    hash_hex = hash_obj.hexdigest()

    # Odpowiednik hex2dec
    h1 = int(hash_hex[0:13], 16)
    h2 = int(hash_hex[13:26], 16)
    h3 = int(hash_hex[26:39], 16)
    h4 = int(hash_hex[39:52], 16)
    h5 = int(hash_hex[52:64], 16) # Do końca hasha (64 znaki)

    x = (h1 + h2 + h5) / ten16
    y = (h3 + h4 + h5) / ten16

    a = 5.0
    b = 5.0
    N = 1.0
    A = 0.84
    B = 0.75
    C = 1.0
    D = 1.0

    # Zapisz klucze. Zwrócimy je do użycia w deszyfrowaniu.
    keys = (x, y, a, b, N, A, B, C, D)
    
    # "Rozgrzewanie" mapy chaotycznej
    for _ in range(49):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)

    xold = x
    yold = y
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    x = np.fmod(x + yold, 1)
    y = np.fmod(xold + y, 1)

    # --- STEP 2 - Arrange... bit planes ---
    x1 = np.zeros(8)
    for i in range(8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        x1[i] = x

    # Odpowiednik [~,indexX]=sort(x1)
    # +1 aby dopasować 1-bazowe indeksowanie bitów (jak w bitget)
    indexX = np.argsort(x1) + 1 

    # Odpowiednik bitget i łączenia
    bit_planes = []
    # Konwertujemy do uint8 TYLKO na czas operacji bitowych
    plaintext_uint8 = plaintext.astype(np.uint8) 
    for bit_num in indexX:
        # (plaintext >> (bit_num - 1)) & 1 to odpowiednik bitget(plaintext, bit_num)
        plane = (plaintext_uint8 >> (bit_num - 1)) & 1
        bit_planes.append(plane)
    
    # Odpowiednik logical([plane1, plane2, ...])
    plainbits = np.hstack(bit_planes).astype(bool)

    # --- STEP 3 & 4 - circshift... ---
    # Oblicz posrow
    posrow = np.zeros(rows)
    for j in range(rows):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        posrow[j] = np.floor(8 * cols * x)
    
    # circshift na posrow
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    pos = int(np.floor(rows * x))
    posrow = np.roll(posrow, shift=pos) # Odpowiednik circshift(posrow, pos)

    # Oblicz poscol
    poscol = np.zeros(cols * 8)
    for j in range(cols * 8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        poscol[j] = np.floor(rows * x)

    # circshift na poscol
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    pos = int(np.floor(cols * 8 * x))
    poscol = np.roll(poscol, shift=pos)

    # Zastosuj circshift wiersz po wierszu
    # (Nie da się tego zwektoryzować, bo każdy wiersz ma inny shift)
    for j in range(rows):
        plainbits[j, :] = np.roll(plainbits[j, :], shift=int(posrow[j]))

    # Zastosuj circshift kolumna po kolumnie
    for j in range(cols * 8):
        plainbits[:, j] = np.roll(plainbits[:, j], shift=int(poscol[j]))

    # --- STEP 5 - Append... bit planes ---
    y1 = np.zeros(8)
    for i in range(8):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        y1[i] = y
    
    # +1 aby dopasować 1-bazowe indeksowanie MATLABa
    indexY = np.argsort(y1) + 1 

    # Rekonstrukcja obrazu (składanie płaszczyzn bitowych)
    plainbitsshuffled = np.zeros((rows, cols), dtype=np.float64)
    for i in range(8):
        # bit_num_to_get to numer płaszczyzny (1-8) z posortowanej listy
        bit_num_to_get = indexY[i] 
        
        # Wyciągnij odpowiednią płaszczyznę z plainbits
        start_col = (bit_num_to_get - 1) * cols
        end_col = bit_num_to_get * cols
        plane = plainbits[:, start_col:end_col]
        
        # Dodaj ją z odpowiednią wagą (2**i)
        plainbitsshuffled += plane.astype(np.float64) * (2**i)

    # --- STEP 6 - Compute the substitution matrix (XOR) ---
    prbgmatvec = np.zeros(rows * cols)
    # Ta pętla nie może być zwektoryzowana z powodu zależności 'y'
    for i in range(rows * cols):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        prbgmatvec[i] = np.fmod(np.floor(DD * y), 256)

    # circshift na wektorze
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    shiftval = int(np.floor(rows * cols * y))
    prbgmatvec = np.roll(prbgmatvec, shift=shiftval)

    # Reshape do macierzy
    prbgmat = prbgmatvec.reshape((rows, cols))

    # Końcowy XOR
    # Musimy konwertować na uint8 na czas operacji bitowej
    encrypted = np.bitwise_xor(plainbitsshuffled.astype(np.uint8), prbgmat.astype(np.uint8))
    
    print("Szyfrowanie zakończone.")
    
    # Zwracamy zaszyfrowany obraz oraz klucze potrzebne do deszyfrowania
    return encrypted, keys, (rows, cols), DD


def decrypt_chaotic(encrypted_image, keys, shape, DD):
    """
    Tłumaczenie sekcji "Decryption Algorithm"
    """
    print("Rozpoczęcie deszyfrowania chaotycznego...")
    
    encrypted = np.array(encrypted_image).astype(np.float64)
    rows, cols = shape
    
    # Rozpakuj klucze
    x, y, a, b, N, A, B, C, D = keys

    # --- "Rozgrzewanie" mapy chaotycznej (identyczne jak w szyfrowaniu) ---
    for _ in range(49):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    
    xold = x
    yold = y
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    x = np.fmod(x + yold, 1)
    y = np.fmod(xold + y, 1)

    # --- Odwrócenie STEP 6 (XOR) ---
    # Najpierw musimy wygenerować ten sam 'prbgmat'
    
    # Oblicz indexY (potrzebny do kroku 5)
    y1 = np.zeros(8)
    for i in range(8):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        y1[i] = y
    indexY = np.argsort(y1) + 1 # 1-based index

    # Generuj prbgmatvec (identyczny kod jak w szyfrowaniu)
    prbgmatvec = np.zeros(rows * cols)
    for i in range(rows * cols):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        prbgmatvec[i] = np.fmod(np.floor(DD * y), 256)

    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    shiftval = int(np.floor(rows * cols * y))
    prbgmatvec = np.roll(prbgmatvec, shift=shiftval)

    prbgmat = prbgmatvec.reshape((rows, cols))

    # Odwróć XOR
    plainbitsshuffled = np.bitwise_xor(encrypted.astype(np.uint8), prbgmat.astype(np.uint8))
    plainbitsshuffled = plainbitsshuffled.astype(np.float64) # Powrót do float dla operacji

    # --- Odwrócenie STEP 5 (Append planes) ---
    # Musimy rozłożyć 'plainbitsshuffled' z powrotem na 'plainbits'
    
    # Wyciągnij 8 płaszczyzn bitowych
    planes_extracted = []
    for i in range(8):
        planes_extracted.append((plainbitsshuffled.astype(np.uint8) >> i) & 1)
        
    plainbits = np.zeros((rows, cols * 8), dtype=bool)
    
    # Użyj indexY, aby umieścić płaszczyzny z powrotem na ich
    # "oryginalnych" pozycjach w 'plainbits'
    for i in range(8):
        original_plane_num = indexY[i] # Np. 2
        start_col = (original_plane_num - 1) * cols
        end_col = original_plane_num * cols
        
        # Płaszczyzna 'i' (np. 0-LSB) idzie na miejsce płaszczyzny '2'
        plainbits[:, start_col:end_col] = planes_extracted[i]
        
    # --- Odwrócenie STEP 2, 3, 4 ---
    
    # Oblicz indexX (identycznie jak w szyfrowaniu)
    x1 = np.zeros(8)
    for i in range(8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        x1[i] = x
    indexX = np.argsort(x1) + 1 # 1-based index

    # Oblicz posrow (identycznie jak w szyfrowaniu)
    posrow = np.zeros(rows)
    for j in range(rows):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        posrow[j] = np.floor(8 * cols * x)
    
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    pos = int(np.floor(rows * x))
    posrow = np.roll(posrow, shift=pos)

    # Oblicz poscol (identycznie jak w szyfrowaniu)
    poscol = np.zeros(cols * 8)
    for j in range(cols * 8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        poscol[j] = np.floor(rows * x)

    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    pos = int(np.floor(cols * 8 * x))
    poscol = np.roll(poscol, shift=pos)

    # Odwróć circshift kolumn (użyj ujemnego shift)
    for j in range(cols * 8):
        plainbits[:, j] = np.roll(plainbits[:, j], shift=-int(poscol[j]))

    # Odwróć circshift wierszy (użyj ujemnego shift)
    for j in range(rows):
        plainbits[j, :] = np.roll(plainbits[j, :], shift=-int(posrow[j]))

    # --- Odwrócenie STEP 2 (Assemble) ---
    # Złóż 'plainbits' z powrotem w ostateczny obraz
    decrypted = np.zeros((rows, cols), dtype=np.uint8)
    
    # Wyciągnij 8 "oryginalnych" płaszczyzn (już posortowanych)
    planes_to_assemble = []
    for i in range(8):
        start_col = i * cols
        end_col = (i + 1) * cols
        planes_to_assemble.append(plainbits[:, start_col:end_col])
        
    # Użyj indexX, aby złożyć je w odpowiedniej kolejności bitowej
    for i in range(8): # i = 0 (LSB), i = 1 (drugi bit), ...
        original_plane_num = i + 1 # Chcemy znaleźć płaszczyznę 1, potem 2, ...
        
        # Znajdź, na którą pozycję 'indexX' trafiła 'original_plane_num'
        # To jest odpowiednik MATLABowego: find(indexX == original_plane_num)
        find_where_it_went = np.where(indexX == original_plane_num)[0][0]
        
        # Weź płaszczyznę z tej pozycji
        plane_to_add = planes_to_assemble[find_where_it_went]
        
        # Dodaj ją z odpowiednią wagą
        decrypted += plane_to_add.astype(np.uint8) * (2**i)
        
    print("Deszyfrowanie zakończone.")
    return decrypted


# =========================================================================
# === GŁÓWNY BLOK WYKONAWCZY (Odpowiednik "Run Section") ===
# =========================================================================
if __name__ == "__main__":
    
    try:
        # Wczytaj obraz, upewnij się, że jest w skali szarości ('L')
        img = Image.open('cameraman.tif').convert('L')
        print(f"Wczytano obraz: {img.size[0]}x{img.size[1]}")
    except FileNotFoundError:
        print("Błąd: Plik 'cameraman.tif' nie został znaleziony.")
        print("Proszę umieścić obraz w tym samym katalogu co skrypt.")
        exit()

    # --- Szyfrowanie ---
    encrypted_img, keys_bundle, shape, DD_val = encrypt_chaotic(img)
    
    # --- Deszyfrowanie ---
    decrypted_img = decrypt_chaotic(encrypted_img, keys_bundle, shape, DD_val)

    # --- Weryfikacja ---
    original_array = np.array(img)
    is_correct = np.array_equal(original_array, decrypted_img)
    
    if is_correct:
        print("\nWeryfikacja: SUKCES! Obraz oryginalny i odszyfrowany są identyczne.")
    else:
        print("\nWeryfikacja: BŁĄD! Obraz oryginalny i odszyfrowany różnią się.")

    # --- Wyświetlanie wyników ---
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    ax1.imshow(original_array, cmap='gray')
    ax1.set_title("Oryginał (Plaintext)")
    
    ax2.imshow(encrypted_img, cmap='gray')
    ax2.set_title("Szyfrogram (Ciphertext)")
    
    ax3.imshow(decrypted_img, cmap='gray')
    ax3.set_title("Odszyfrowany (Decrypted)")
    
    plt.show()
