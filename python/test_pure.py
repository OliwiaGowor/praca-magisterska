import numpy as np
import hashlib
from PIL import Image
import matplotlib.pyplot as plt
import time
from scipy.io import savemat
from skimage.measure import shannon_entropy
from skimage.metrics import peak_signal_noise_ratio as psnr, structural_similarity as ssim

# === IMPORTY PURE PYTHON ===
import pyaes
# (Usunięto Blowfish i ChaCha20)

# =========================================================================
# === TŁUMACZENIE FUNKCJI POMOCNICZYCH ===
# =========================================================================

# === USUNIĘTO NUMBA JIT ===
def smht(x, a, b, c, d):
    """ Port 1:1 funkcji smht.m (TERAZ CZYSTY PYTHON) """
    numerator = np.exp(a * x) - np.exp(-b * x)
    denominator = np.exp(c * x) + np.exp(-d * x)
    if denominator == 0:
        return 0.0
    return numerator / denominator

def calculate_entropy(img):
    return shannon_entropy(np.uint8(img))

def calculate_correlation(img, num_samples=5000):
    img_d = img.astype(np.float64)
    rows, cols = img_d.shape
    try:
        x_h = np.random.randint(0, rows, num_samples)
        y_h = np.random.randint(0, cols - 1, num_samples)
        v1_h = img_d[x_h, y_h]
        v2_h = img_d[x_h, y_h + 1]
        corr_h = np.corrcoef(v1_h, v2_h)[0, 1]
        x_v = np.random.randint(0, rows - 1, num_samples)
        y_v = np.random.randint(0, cols, num_samples)
        v1_v = img_d[x_v, y_v]
        v2_v = img_d[x_v + 1, y_v]
        corr_v = np.corrcoef(v1_v, v2_v)[0, 1]
        x_d = np.random.randint(0, rows - 1, num_samples)
        y_d = np.random.randint(0, cols - 1, num_samples)
        v1_d = img_d[x_d, y_d]
        v2_d = img_d[x_d + 1, y_d + 1]
        corr_d = np.corrcoef(v1_d, v2_d)[0, 1]
    except (np.linalg.LinAlgError, ValueError):
        corr_h, corr_v, corr_d = 0.0, 0.0, 0.0
    if np.isnan(corr_h): corr_h = 0.0
    if np.isnan(corr_v): corr_v = 0.0
    if np.isnan(corr_d): corr_d = 0.0
    return corr_h, corr_v, corr_d


def calculate_npcr_uaci(c1, c2):
    c1_d = c1.astype(np.float64)
    c2_d = c2.astype(np.float64)
    rows, cols = c1.shape
    num_pixels = rows * cols
    diff = (c1_d != c2_d)
    npcr = np.sum(diff) / num_pixels * 100
    uaci = np.sum(np.abs(c1_d - c2_d)) / (num_pixels * 255) * 100
    return npcr, uaci

# === FUNKCJE POMOCNICZE: Padding i Chunker ===
def pad(data_to_pad, block_size):
    """Dodaje dopełnienie PKCS#7"""
    padding_len = block_size - (len(data_to_pad) % block_size)
    padding = bytes([padding_len] * padding_len)
    return data_to_pad + padding

def unpad(padded_data, block_size):
    """Usuwa dopełnienie PKCS#7"""
    padding_len = padded_data[-1]
    if padding_len > block_size or padding_len == 0:
        print("BŁĄD: Niepoprawny padding podczas unpad!")
        return padded_data 
    return padded_data[:-padding_len]

def chunker(data, block_size):
    """Generator dzielący dane na bloki (używany tylko przez pyaes)"""
    for i in range(0, len(data), block_size):
        yield data[i:i + block_size]

# =========================================================================
# === TŁUMACZENIE ALGORYTMU CHAOTYCZNEGO (WERSJA BEZ JIT) ===
# =========================================================================

def encrypt_circular_chaotic(plaintext_image_pil):
    """ Tłumaczenie 'encrypt_circular_chaotic.m' (Czysty Python, bez JIT) """
    
    plaintext = np.array(plaintext_image_pil).astype(np.float64)
    rows, cols = plaintext.shape
    DD = 10**10
    ten16 = 10**16

    # --- STEP 1: Key initialization ---
    hash_obj = hashlib.sha256(plaintext.tobytes('C'))
    hash_hex = hash_obj.hexdigest()
    h1 = int(hash_hex[0:13], 16); h2 = int(hash_hex[13:26], 16)
    h3 = int(hash_hex[26:39], 16); h4 = int(hash_hex[39:52], 16)
    h5 = int(hash_hex[52:64], 16) 
    x_init = (h1 + h2 + h5) / ten16; y_init = (h3 + h4 + h5) / ten16
    a=5.0; b=5.0; N=1.0; A=0.84; B=0.75; C=1.0; D=1.0
    keys_tuple = (x_init, y_init, a, b, N, A, B, C, D) 
    
    x, y = x_init, y_init
    for _ in range(49):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    xold=x; yold=y
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    x=np.fmod(x+yold,1); y=np.fmod(xold+y,1)

    # --- STEP 2: Arrange bit planes ---
    x1 = np.zeros(8)
    for i in range(8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        x1[i] = x
    indexX = np.argsort(x1) + 1 
    plaintext_uint8 = plaintext.astype(np.uint8)
    bit_planes = []
    for bit_num in indexX:
        plane = (plaintext_uint8 >> (bit_num - 1)) & 1
        bit_planes.append(plane)
    plainbits = np.hstack(bit_planes).astype(bool)

    # --- STEP 3 & 4: circshift ---
    posrow = np.zeros(rows)
    for j in range(rows):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        posrow[j] = np.floor(8 * cols * x)
    
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    posrow = np.roll(posrow, shift=int(np.floor(rows * x)))

    poscol = np.zeros(cols * 8)
    for j in range(cols * 8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        poscol[j] = np.floor(rows * x)

    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    poscol = np.roll(poscol, shift=int(np.floor(cols * 8 * x)))

    for j in range(rows):
        plainbits[j, :] = np.roll(plainbits[j, :], shift=int(posrow[j]))
    for j in range(cols * 8):
        plainbits[:, j] = np.roll(plainbits[:, j], shift=int(poscol[j]))

    # --- STEP 5: Append bit planes ---
    y1 = np.zeros(8)
    for i in range(8):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        y1[i] = y
    indexY = np.argsort(y1) + 1 

    plainbitsshuffled = np.zeros((rows, cols), dtype=np.float64)
    for i in range(8):
        bit_num_to_get = indexY[i] 
        start_col = (bit_num_to_get - 1) * cols
        end_col = bit_num_to_get * cols
        plane = plainbits[:, start_col:end_col]
        plainbitsshuffled += plane.astype(np.float64) * (2**i)

    # --- STEP 6: Substitution (XOR) ---
    prbgmatvec = np.zeros(rows * cols)
    for i in range(rows * cols):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        prbgmatvec[i] = np.fmod(np.floor(DD * y), 256)

    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    prbgmatvec = np.roll(prbgmatvec, shift=int(np.floor(rows * cols * y)))
    prbgmat = prbgmatvec.reshape((rows, cols))

    encrypted = np.bitwise_xor(plainbitsshuffled.astype(np.uint8), prbgmat.astype(np.uint8))
    
    keys_struct = {
        'keys_tuple': keys_tuple,
        'shape': (rows, cols),
        'DD': DD
    }
    return encrypted, keys_struct


def decrypt_circular_chaotic(encrypted_image, keys_struct):
    """ Tłumaczenie 'decrypt_circular_chaotic.m' (Czysty Python, bez JIT) """
    
    encrypted = np.array(encrypted_image).astype(np.float64)
    keys_tuple = keys_struct['keys_tuple']
    rows, cols = keys_struct['shape']
    DD = keys_struct['DD']
    x_init, y_init, a, b, N, A, B, C, D = keys_tuple

    # "Rozgrzewanie"
    x, y = x_init, y_init
    for _ in range(49):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    xold=x; yold=y
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    x=np.fmod(x+yold,1); y=np.fmod(xold+y,1)
    
    # --- Odwrócenie Kroków (Generowanie kluczy) ---
    x1 = np.zeros(8)
    for i in range(8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        x1[i] = x
    indexX = np.argsort(x1) + 1 

    posrow = np.zeros(rows)
    for j in range(rows):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        posrow[j] = np.floor(8 * cols * x)
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    posrow = np.roll(posrow, shift=int(np.floor(rows * x)))

    poscol = np.zeros(cols * 8)
    for j in range(cols * 8):
        x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
        poscol[j] = np.floor(rows * x)
    x = np.fmod(x + a + b * smht(x / N, A, B, C, D), N)
    poscol = np.roll(poscol, shift=int(np.floor(cols * 8 * x)))
    
    y1 = np.zeros(8)
    for i in range(8):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        y1[i] = y
    indexY = np.argsort(y1) + 1 
        
    prbgmatvec = np.zeros(rows * cols)
    for i in range(rows * cols):
        y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
        prbgmatvec[i] = np.fmod(np.floor(DD * y), 256)
    y = np.fmod(y + a + b * smht(y / N, A, B, C, D), N)
    prbgmat_shift_val = int(np.floor(rows * cols * y))

    # --- Odwrócenie Krok 6 (XOR) ---
    prbgmatvec = np.roll(prbgmatvec, shift=prbgmat_shift_val)
    prbgmat = prbgmatvec.reshape((rows, cols))
    plainbitsshuffled = np.bitwise_xor(encrypted.astype(np.uint8), prbgmat.astype(np.uint8))
    plainbitsshuffled = plainbitsshuffled.astype(np.float64)

    # --- Odwrócenie Krok 5 (Append planes) ---
    planes_extracted = []
    for i in range(8):
        planes_extracted.append((plainbitsshuffled.astype(np.uint8) >> i) & 1)
    plainbits = np.zeros((rows, cols * 8), dtype=bool)
    for i in range(8):
        original_plane_num = indexY[i] 
        start_col = (original_plane_num - 1) * cols
        end_col = original_plane_num * cols
        plainbits[:, start_col:end_col] = planes_extracted[i]
        
    # --- Odwrócenie Krok 3/4 (circshift) ---
    for j in range(cols * 8):
        plainbits[:, j] = np.roll(plainbits[:, j], shift=-int(poscol[j]))
    for j in range(rows):
        plainbits[j, :] = np.roll(plainbits[j, :], shift=-int(posrow[j]))

    # --- Odwrócenie Krok 2 (Assemble) ---
    decrypted = np.zeros((rows, cols), dtype=np.uint8)
    planes_to_assemble = []
    for i in range(8):
        start_col = i * cols
        end_col = (i + 1) * cols
        planes_to_assemble.append(plainbits[:, start_col:end_col])
    for i in range(8): 
        original_plane_num = i + 1
        find_where_it_went = np.where(indexX == original_plane_num)[0][0]
        plane_to_add = planes_to_assemble[find_where_it_went]
        decrypted += plane_to_add.astype(np.uint8) * (2**i)
        
    return decrypted

# =========================================================================
# === GŁÓWNY SKRYPT PORÓWNAWCZY (Pure Python vs Pure Python) ===
# =========================================================================
def main():
    
    N_RUNS = 3 # UWAGA: Ustawiono na 3, ponieważ testy są BARDZO wolne
    
    print(f'Rozpoczynanie {N_RUNS} przebiegów testowych (Pure Python vs Pure Python)...')
    print('*** TO MOŻE POTRWAĆ KILKA MINUT ***')
    
    try:
        img_orig_pil = Image.open('cameraman.tif').convert('L')
        img_orig = np.array(img_orig_pil)
    except FileNotFoundError:
        print("BŁĄD: Nie znaleziono 'cameraman.tif'.")
        img_orig = np.random.randint(0, 256, (256, 256), dtype=np.uint8)
        img_orig_pil = Image.fromarray(img_orig)
        
    img_size = img_orig.shape
    rows, cols = img_size

    entropy_orig = calculate_entropy(img_orig)
    corr_h_orig, corr_v_orig, corr_d_orig = calculate_correlation(img_orig)

    img_mod = img_orig.copy()
    img_mod[0, 0] = np.bitwise_xor(img_mod[0, 0], 1)
    img_mod_pil = Image.fromarray(img_mod) 

    titles = ['AES-CBC (Pure Py)', 'Chaotyczny (Pure Py)']
    num_algos = len(titles)

    wall_times_enc = np.full((num_algos, N_RUNS), np.nan)
    wall_times_dec = np.full((num_algos, N_RUNS), np.nan)
    cpu_times_enc = np.full((num_algos, N_RUNS), np.nan)
    cpu_times_dec = np.full((num_algos, N_RUNS), np.nan)
    
    npcr_values = np.full((num_algos, N_RUNS), np.nan)
    uaci_values = np.full((num_algos, N_RUNS), np.nan)
    entropy_values = np.full((num_algos, N_RUNS), np.nan)
    corr_h_values = np.full((num_algos, N_RUNS), np.nan)
    corr_v_values = np.full((num_algos, N_RUNS), np.nan)
    corr_d_values = np.full((num_algos, N_RUNS), np.nan)

    images_enc = [None] * num_algos
    images_dec = [None] * num_algos

    img_bytes_orig = img_orig.tobytes('C')
    img_bytes_mod = img_mod.tobytes('C')
    
    for run_idx in range(N_RUNS):
        print(f'--- Rozpoczynanie przebiegu {run_idx + 1} / {N_RUNS} ---')
        
        # --- TEST 1: AES-256-CBC (pyaes) ---
        print("Testowanie AES (pyaes)...")
        key_aes = np.random.bytes(32) # 256 bitów
        iv_aes = np.random.bytes(16)
        
        img_bytes_aes_pad = pad(img_bytes_orig, 16)
        img_bytes_aes_mod_pad = pad(img_bytes_mod, 16)
        
        start_wall = time.perf_counter(); start_cpu = time.process_time()
        cipher_enc_aes = pyaes.AESModeOfOperationCBC(key_aes, iv_aes)
        ct_chunks_aes = []
        for chunk in chunker(img_bytes_aes_pad, 16):
            ct_chunks_aes.append(cipher_enc_aes.encrypt(chunk))
        ct_aes_C1_padded = b"".join(ct_chunks_aes)
        wall_times_enc[0, run_idx] = time.perf_counter() - start_wall
        cpu_times_enc[0, run_idx] = time.process_time() - start_cpu
        
        cipher_enc_aes_2 = pyaes.AESModeOfOperationCBC(key_aes, iv_aes)
        ct_chunks_aes_2 = []
        for chunk in chunker(img_bytes_aes_mod_pad, 16):
            ct_chunks_aes_2.append(cipher_enc_aes_2.encrypt(chunk))
        ct_aes_C2_padded = b"".join(ct_chunks_aes_2)
        
        start_wall = time.perf_counter(); start_cpu = time.process_time()
        cipher_dec_aes = pyaes.AESModeOfOperationCBC(key_aes, iv_aes)
        pt_chunks_aes = []
        for chunk in chunker(ct_aes_C1_padded, 16):
            pt_chunks_aes.append(cipher_dec_aes.decrypt(chunk))
        pt_aes_padded = b"".join(pt_chunks_aes)
        pt_aes_bytes = unpad(pt_aes_padded, 16)
        wall_times_dec[0, run_idx] = time.perf_counter() - start_wall
        cpu_times_dec[0, run_idx] = time.process_time() - start_cpu
        
        C1 = np.frombuffer(ct_aes_C1_padded, dtype=np.uint8)[:rows*cols].reshape(img_size)
        C2 = np.frombuffer(ct_aes_C2_padded, dtype=np.uint8)[:rows*cols].reshape(img_size)
        npcr_values[0, run_idx], uaci_values[0, run_idx] = calculate_npcr_uaci(C1, C2)
        entropy_values[0, run_idx] = calculate_entropy(C1)
        corr_h_values[0, run_idx], corr_v_values[0, run_idx], corr_d_values[0, run_idx] = calculate_correlation(C1)
        images_enc[0] = C1
        images_dec[0] = np.frombuffer(pt_aes_bytes, dtype=np.uint8).reshape(img_size)
        

        # --- TEST 3: ALGORYTM CHAOTYCZNY (Pure Python) ---
        print("Testowanie Chaotyczny (Pure Python)...")
        try:
            start_wall = time.perf_counter(); start_cpu = time.process_time()
            C1_chaos, keys_chaos_C1 = encrypt_circular_chaotic(img_orig_pil)
            wall_times_enc[2, run_idx] = time.perf_counter() - start_wall
            cpu_times_enc[2, run_idx] = time.process_time() - start_cpu
            
            C2_chaos, keys_chaos_C2 = encrypt_circular_chaotic(img_mod_pil)
            
            start_wall = time.perf_counter(); start_cpu = time.process_time()
            pt_chaos = decrypt_circular_chaotic(C1_chaos, keys_chaos_C1)
            wall_times_dec[2, run_idx] = time.perf_counter() - start_wall
            cpu_times_dec[2, run_idx] = time.process_time() - start_cpu
            
            if not np.array_equal(img_orig, pt_chaos):
                print(f"BŁĄD KRYTYCZNY (Przebieg {run_idx+1}): Deszyfrowanie chaotyczne nie powiodło się!")
            
            npcr_values[2, run_idx], uaci_values[2, run_idx] = calculate_npcr_uaci(C1_chaos, C2_chaos)
            entropy_values[2, run_idx] = calculate_entropy(C1_chaos)
            corr_h_values[2, run_idx], corr_v_values[2, run_idx], corr_d_values[2, run_idx] = calculate_correlation(C1_chaos)
            images_enc[2] = C1_chaos
            images_dec[2] = pt_chaos
        except Exception as e:
            print(f"BŁĄD podczas testu algorytmu chaotycznego: {e}")

    print('--- Wszystkie przebiegi zakończone! ---\n')

    # --- OBLICZANIE ŚREDNICH I ZAPIS DO PLIKU ---
    avg_wall_times_enc = np.nanmean(wall_times_enc, axis=1)
    avg_wall_times_dec = np.nanmean(wall_times_dec, axis=1)
    avg_cpu_times_enc = np.nanmean(cpu_times_enc, axis=1)
    avg_cpu_times_dec = np.nanmean(cpu_times_dec, axis=1)
    
    avg_npcr = np.nanmean(npcr_values, axis=1)
    avg_uaci = np.nanmean(uaci_values, axis=1)
    avg_entropy = np.nanmean(entropy_values, axis=1)
    avg_corr_h = np.nanmean(corr_h_values, axis=1)
    avg_corr_v = np.nanmean(corr_v_values, axis=1)
    avg_corr_d = np.nanmean(corr_d_values, axis=1)
    
    filename = 'encryption_results_PURE_PYTHON.mat'
    results_dict = {
        'titles': titles, 'N_RUNS': N_RUNS,
        'wall_times_enc': wall_times_enc, 'wall_times_dec': wall_times_dec,
        'cpu_times_enc': cpu_times_enc, 'cpu_times_dec': cpu_times_dec,
        'npcr_values': npcr_values, 'uaci_values': uaci_values,
        'entropy_values': entropy_values,
        'corr_h_values': corr_h_values, 'corr_v_values': corr_v_values, 'corr_d_values': corr_d_values,
        'entropy_orig': entropy_orig, 'corr_h_orig': corr_h_orig, 'corr_v_orig': corr_v_orig, 'corr_d_orig': corr_d_orig,
        'avg_wall_times_enc': avg_wall_times_enc, 'avg_wall_times_dec': avg_wall_times_dec,
        'avg_cpu_times_enc': avg_cpu_times_enc, 'avg_cpu_times_dec': avg_cpu_times_dec,
        'avg_npcr': avg_npcr, 'avg_uaci': avg_uaci,
        'avg_entropy': avg_entropy,
        'avg_corr_h': avg_corr_h, 'avg_corr_v': avg_corr_v, 'avg_corr_d': avg_corr_d
    }
    savemat(filename, results_dict)
    print(f'Wszystkie surowe wyniki ({N_RUNS} przebiegów) zostały zapisane w pliku: {filename}')

    # --- WIZUALIZACJA ŚREDNICH WYNIKÓW ---
    print('Generowanie wizualizacji na podstawie średnich wyników...')
    subplot_height = num_algos + 2 
    fig = plt.figure(figsize=(18, 12 + num_algos * 1.5))
    
    ax_bar = plt.subplot2grid((subplot_height, 3), (0, 1), colspan=2)
    bar_width = 0.35
    index = np.arange(num_algos)
    
    ax_bar.bar(index - bar_width/2, avg_cpu_times_enc, bar_width, label='Szyfrowanie (CPU Time)')
    ax_bar.bar(index + bar_width/2, avg_cpu_times_dec, bar_width, label='Deszyfrowanie (CPU Time)')
    
    ax_bar.set_ylabel('Czas CPU (sekundy)')
    ax_bar.set_title(f'Porównanie OBCIĄŻENIA CPU (Pure Python, Średnia z {N_RUNS} przebiegów, Obraz: {rows}x{cols})')
    ax_bar.set_xticks(index)
    ax_bar.set_xticklabels(titles, rotation=10)
    ax_bar.legend(loc='upper left')
    if (not np.isnan(avg_cpu_times_enc).all()) and (np.nanmax(avg_cpu_times_enc) / np.nanmin(avg_cpu_times_enc[avg_cpu_times_enc>0])) > 50:
         ax_bar.set_yscale('log')
         ax_bar.set_ylabel('Czas CPU (sekundy) - Skala Log')
    ax_bar.grid(axis='y', linestyle='--', alpha=0.7)
    
    ax_orig = plt.subplot2grid((subplot_height, 3), (1, 0))
    ax_orig.imshow(img_orig, cmap='gray', vmin=0, vmax=255)
    title_str = f'Oryginał\nEntropia: {entropy_orig:.4f}\nKorelacja H/V/D: {corr_h_orig:.3f} / {corr_v_orig:.3f} / {corr_d_orig:.3f}'
    ax_orig.set_title(title_str)
    ax_orig.set_ylabel('WEJŚCIE', fontweight='bold')
    
    ax_mod = plt.subplot2grid((subplot_height, 3), (2, 0))
    ax_mod.imshow(img_mod, cmap='gray', vmin=0, vmax=255)
    ax_mod.set_title('Oryginał (zmieniony 1px)')
    ax_mod.set_ylabel('WEJŚCIE MOD', fontweight='bold')

    for i in range(num_algos):
        row_idx = i + 2 
        ax_enc = plt.subplot2grid((subplot_height, 3), (row_idx, 1))
        if images_enc[i] is not None:
            ax_enc.imshow(images_enc[i], cmap='gray', vmin=0, vmax=255)
        
        title_str = (
            f'Szyfrogram (Średni Czas [Wall]: {avg_wall_times_enc[i]:.4f} s)\n'
            f'Śr. Ent: {avg_entropy[i]:.4f}, Śr. NPCR: {avg_npcr[i]:.2f}%, Śr. UACI: {avg_uaci[i]:.2f}%\n'
            f'Śr. Korelacja H/V/D: {avg_corr_h[i]:.3f} / {avg_corr_v[i]:.3f} / {avg_corr_d[i]:.3f}'
        )
        ax_enc.set_title(title_str, fontsize=9)
        ax_enc.set_ylabel(titles[i], fontweight='bold')
        
        ax_dec = plt.subplot2grid((subplot_height, 3), (row_idx, 2))
        if images_dec[i] is not None:
            ax_dec.imshow(images_dec[i], cmap='gray', vmin=0, vmax=255)
        ax_dec.set_title(f'Deszyfrogram\n(Średni Czas [Wall]: {avg_wall_times_dec[i]:.4f} s)')

    plt.subplot2grid((subplot_height, 3), (0, 0)).axis('off')
    plt.subplot2grid((subplot_height, 3), (1, 1)).axis('off')
    plt.subplot2grid((subplot_height, 3), (1, 2)).axis('off')
    if num_algos > 1: plt.subplot2grid((subplot_height, 3), (3, 0)).axis('off')
    if num_algos > 2: plt.subplot2grid((subplot_height, 3), (4, 0)).axis('off')
    if num_algos > 3: plt.subplot2grid((subplot_height, 3), (5, 0)).axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle('Kompletne Porównanie Kryptograficzne (Pure Python)', fontsize=16, fontweight='bold')
    plt.savefig('encryption_benchmark_PURE_PYTHON.png')
    print('Wykresy zostały zapisane do pliku: encryption_benchmark_PURE_PYTHON.png')
    plt.show()

if __name__ == "__main__":
    main()
