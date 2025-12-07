#include "mex.h"
#include <openssl/evp.h>
#include <openssl/provider.h> // <-- NOWY INCLUDE DLA OPENSSL 3.0+
#include <string.h> 

/*
 * Uniwersalny wrapper kryptograficzny dla MATLABa
 * Wersja 3 (Robust): Ładuje "legacy provider" dla OpenSSL 3.0+
 */

// --- Globalna flaga, aby załadować dostawców tylko raz ---
static int openssl_initialized = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // --- KROK 0: Jednorazowa inicjalizacja OpenSSL (dla Blowfish, RC4, etc.) ---
    // To zostanie wykonane tylko raz, przy pierwszym wywołaniu tej funkcji MEX
    if (!openssl_initialized) {
        // Jawnie ładujemy "legacy" (dla Blowfish, RC4, DES)
        // i "default" (dla AES, ChaCha20).
        // W OpenSSL < 3.0 te funkcje nie istnieją, ale kod i tak
        // będzie działał (po prostu się nie załadują).
        // W OpenSSL 3.0+ to jest WYMAGANE.
        OSSL_PROVIDER_load(NULL, "legacy");
        OSSL_PROVIDER_load(NULL, "default");
        
        openssl_initialized = 1; 
        // Ustaw flagę, nawet jeśli ładowanie się nie powiodło,
        // aby nie próbować w kółko.
    }

    // --- 1. Walidacja wejść ---
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CRYPTO:nrhs", "Wymagane 5 argumentów: (plaintext, key, iv, cipher_name, do_encrypt)");
    }
    // ... (reszta walidacji bez zmian) ...
    if (!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]) || !mxIsUint8(prhs[2])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "plaintext, key i iv muszą być typu uint8.");
    }
    if (!mxIsChar(prhs[3])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "Nazwa szyfru (cipher_name) musi być łańcuchem znaków (string).");
    }
    if (!mxIsLogicalScalar(prhs[4])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "do_encrypt musi być wartością logiczną (true/false).");
    }

    // --- 2. Pobranie danych z MATLAB-a ---
    const unsigned char* plaintext = (const unsigned char*)mxGetData(prhs[0]);
    size_t plaintext_len = mxGetNumberOfElements(prhs[0]);
    const unsigned char* key = (const unsigned char*)mxGetData(prhs[1]);
    size_t key_len = mxGetNumberOfElements(prhs[1]);
    const unsigned char* iv = (const unsigned char*)mxGetData(prhs[2]);
    size_t iv_len = mxGetNumberOfElements(prhs[2]);
    char cipher_name[64];
    mxGetString(prhs[3], cipher_name, sizeof(cipher_name));
    int do_encrypt = (int)mxIsLogicalScalarTrue(prhs[4]);
    
    // --- 3. Wybór szyfru OpenSSL ---
    const EVP_CIPHER* cipher = EVP_get_cipherbyname(cipher_name);
    
    if (cipher == NULL) {
        // Jeśli algorytm nie jest znaleziony, to jest teraz główny podejrzany
        mexErrMsgIdAndTxt("CRYPTO:cipher", "Nieznany szyfr: %s. (Jeśli to Blowfish/RC4, upewnij się, że masz OpenSSL z 'legacy provider').", cipher_name);
    }

    // --- 4. Walidacja klucza i IV (bez zmian) ---
    int required_iv_len = EVP_CIPHER_iv_length(cipher);
    if (required_iv_len > 0 && required_iv_len != iv_len) {
         mexErrMsgIdAndTxt("CRYPTO:iv", "Zła długość IV. Szyfr %s wymaga %d bajtów, podano %d.", 
                           cipher_name, required_iv_len, iv_len);
    }
    int cipher_flags = EVP_CIPHER_flags(cipher);
    if (!(cipher_flags & EVP_CIPH_VARIABLE_LENGTH)) {
        if (EVP_CIPHER_key_length(cipher) != key_len) {
            mexErrMsgIdAndTxt("CRYPTO:key", "Zła długość klucza. Szyfr %s wymaga DOKŁADNIE %d bajtów, podano %d.", 
                              cipher_name, EVP_CIPHER_key_length(cipher), key_len);
        }
    }

    // --- 5. Alokacja pamięci wyjściowej ---
    size_t max_output_len = plaintext_len + EVP_CIPHER_block_size(cipher);
    plhs[0] = mxCreateNumericMatrix(1, max_output_len, mxUINT8_CLASS, mxREAL);
    unsigned char* output = (unsigned char*)mxGetData(plhs[0]);
    int output_len = 0;
    int final_len = 0;
    
    // --- 6. Wykonanie operacji (ta 3-etapowa inicjalizacja jest poprawna) ---
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (ctx == NULL) mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CIPHER_CTX_new failed");

    // Krok 1: Inicjalizuj typ szyfru ORAZ tryb (enc/dec)
    if (1 != EVP_CipherInit_ex(ctx, cipher, NULL, NULL, NULL, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        // Ten błąd (Krok 1 failed) był spowodowany brakiem providera, co naprawiliśmy w Kroku 0
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherInit_ex (Krok 1: set cipher) failed. Prawdopodobnie błąd providera OpenSSL.");
    }

    // Krok 2: JAWNIE ustaw długość klucza
    if (1 != EVP_CIPHER_CTX_set_key_length(ctx, (int)key_len)) {
         EVP_CIPHER_CTX_free(ctx);
         mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CIPHER_CTX_set_key_length failed");
    }
    
    // Krok 3: Ustaw klucz i IV
    if (1 != EVP_CipherInit_ex(ctx, NULL, NULL, key, iv, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherInit_ex (Krok 3: set key/iv) failed");
    }

    // Ustawienie paddingu (bez zmian)
    if (EVP_CIPHER_mode(cipher) == EVP_CIPH_STREAM_CIPHER || 
        EVP_CIPHER_mode(cipher) == EVP_CIPH_CTR_MODE) {
        EVP_CIPHER_CTX_set_padding(ctx, 0); 
    }

    // Przetwarzanie (bez zmian)
    if (1 != EVP_CipherUpdate(ctx, output, &output_len, plaintext, (int)plaintext_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherUpdate failed");
    }
    if (1 != EVP_CipherFinal_ex(ctx, output + output_len, &final_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherFinal_ex failed (sprawdź klucz/padding)");
    }
    EVP_CIPHER_CTX_free(ctx);
    
    // --- 7. Ustawienie finalnego rozmiaru wyjścia ---
    size_t total_len = (size_t)(output_len + final_len);
    mxSetN(plhs[0], total_len);
}