 
#include "mex.h"
#include <openssl/evp.h>
#include <string.h> // Dla strcmp

/*
 * Funkcja pomocnicza do wyboru szyfru na podstawie nazwy trybu i długości klucza
 */
const EVP_CIPHER* get_evp_cipher(size_t key_len, const char* mode) {
    if (key_len == 16) { // AES-128
        if (strcmp(mode, "CBC") == 0) return EVP_aes_128_cbc();
        if (strcmp(mode, "ECB") == 0) return EVP_aes_128_ecb();
        if (strcmp(mode, "CTR") == 0) return EVP_aes_128_ctr();
    } else if (key_len == 24) { // AES-192
        if (strcmp(mode, "CBC") == 0) return EVP_aes_192_cbc();
        if (strcmp(mode, "ECB") == 0) return EVP_aes_192_ecb();
        if (strcmp(mode, "CTR") == 0) return EVP_aes_192_ctr();
    } else if (key_len == 32) { // AES-256
        if (strcmp(mode, "CBC") == 0) return EVP_aes_256_cbc();
        if (strcmp(mode, "ECB") == 0) return EVP_aes_256_ecb();
        if (strcmp(mode, "CTR") == 0) return EVP_aes_256_ctr();
    }
    return NULL; // Nieobsługiwany tryb lub długość klucza
}

/*
 * Główna funkcja MEX
 * Wejścia: (plaintext, key, iv, mode_str, do_encrypt)
 * Wyjścia: [ciphertext]
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --- 1. Walidacja wejść ---
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("AES:nrhs", "Wymagane 5 argumentów: (plaintext, key, iv, mode, do_encrypt)");
    }
    if (!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]) || !mxIsUint8(prhs[2])) {
        mexErrMsgIdAndTxt("AES:type", "plaintext, key i iv muszą być typu uint8.");
    }
    if (!mxIsChar(prhs[3])) {
        mexErrMsgIdAndTxt("AES:type", "Tryb (mode) musi być łańcuchem znaków (string).");
    }
    if (!mxIsLogicalScalar(prhs[4])) {
        mexErrMsgIdAndTxt("AES:type", "do_encrypt musi być wartością logiczną (true/false).");
    }

    // --- 2. Pobranie danych z MATLAB-a ---
    const unsigned char* plaintext = (const unsigned char*)mxGetData(prhs[0]);
    size_t plaintext_len = mxGetNumberOfElements(prhs[0]);

    const unsigned char* key = (const unsigned char*)mxGetData(prhs[1]);
    size_t key_len = mxGetNumberOfElements(prhs[1]);

    const unsigned char* iv = (const unsigned char*)mxGetData(prhs[2]);
    size_t iv_len = mxGetNumberOfElements(prhs[2]);

    char mode_str[8]; // Bufor na nazwę trybu
    mxGetString(prhs[3], mode_str, sizeof(mode_str));

    int do_encrypt = (int)mxIsLogicalScalarTrue(prhs[4]);

    if (iv_len != 16 && strcmp(mode_str, "ECB") != 0) {
         mexErrMsgIdAndTxt("AES:iv", "IV musi mieć 16 bajtów dla trybów innych niż ECB.");
    }

    // --- 3. Wybór szyfru OpenSSL ---
    const EVP_CIPHER* cipher = get_evp_cipher(key_len, mode_str);
    if (cipher == NULL) {
        mexErrMsgIdAndTxt("AES:mode", "Nieobsługiwana kombinacja trybu lub długości klucza.");
    }

    // --- 4. Alokacja pamięci wyjściowej ---
    size_t max_output_len = plaintext_len + EVP_CIPHER_block_size(cipher);
    plhs[0] = mxCreateNumericMatrix(1, max_output_len, mxUINT8_CLASS, mxREAL);
    unsigned char* output = (unsigned char*)mxGetData(plhs[0]);

    int output_len = 0;
    int final_len = 0;

    // --- 5. Wykonanie operacji AES ---
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (ctx == NULL) mexErrMsgIdAndTxt("AES:openssl", "EVP_CIPHER_CTX_new failed");

    if (1 != EVP_CipherInit_ex(ctx, cipher, NULL, key, iv, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("AES:openssl", "EVP_CipherInit_ex failed");
    }

    if (strcmp(mode_str, "CTR") == 0) {
        EVP_CIPHER_CTX_set_padding(ctx, 0); 
    }

    if (1 != EVP_CipherUpdate(ctx, output, &output_len, plaintext, (int)plaintext_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("AES:openssl", "EVP_CipherUpdate failed");
    }

    if (1 != EVP_CipherFinal_ex(ctx, output + output_len, &final_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("AES:openssl", "EVP_CipherFinal_ex failed (sprawdź klucz/padding)");
    }

    EVP_CIPHER_CTX_free(ctx);

    // --- 6. Ustawienie finalnego rozmiaru wyjścia w MATLAB-ie ---
    size_t total_len = (size_t)(output_len + final_len);
    mxSetN(plhs[0], total_len); // Zmniejszamy bufor do rzeczywistego rozmiaru
}
