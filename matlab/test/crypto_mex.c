#include "mex.h"
#include <openssl/evp.h>
#include <openssl/provider.h> 
#include <string.h> 
#include <openssl/err.h> // Dodano do czyszczenia kolejki błędów

/*
 * Universal crypto wrapper for MATLAB
 * Version 4 (Attack-Resilient): Ignores padding errors during decryption
 */

static int openssl_initialized = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // --- STEP 0: Initialization ---
    if (!openssl_initialized) {
        OSSL_PROVIDER_load(NULL, "legacy");
        OSSL_PROVIDER_load(NULL, "default");
        openssl_initialized = 1; 
    }

    // --- 1. Input validation ---
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CRYPTO:nrhs", "Wymagane 5 argumentów.");
    }
    
    if (!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]) || !mxIsUint8(prhs[2])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "plaintext, key i iv muszą być typu uint8.");
    }
    
    // --- 2. Get data ---
    const unsigned char* plaintext = (const unsigned char*)mxGetData(prhs[0]);
    size_t plaintext_len = mxGetNumberOfElements(prhs[0]);
    const unsigned char* key = (const unsigned char*)mxGetData(prhs[1]);
    size_t key_len = mxGetNumberOfElements(prhs[1]);
    const unsigned char* iv = (const unsigned char*)mxGetData(prhs[2]);
    size_t iv_len = mxGetNumberOfElements(prhs[2]);
    
    char cipher_name[64];
    mxGetString(prhs[3], cipher_name, sizeof(cipher_name));
    int do_encrypt = (int)mxIsLogicalScalarTrue(prhs[4]);
    
    // --- 3. Cipher Setup ---
    const EVP_CIPHER* cipher = EVP_get_cipherbyname(cipher_name);
    if (cipher == NULL) {
        mexErrMsgIdAndTxt("CRYPTO:cipher", "Nieznany szyfr: %s", cipher_name);
    }

    // --- 4. Allocate output ---
    size_t max_output_len = plaintext_len + EVP_CIPHER_block_size(cipher);
    plhs[0] = mxCreateNumericMatrix(1, max_output_len, mxUINT8_CLASS, mxREAL);
    unsigned char* output = (unsigned char*)mxGetData(plhs[0]);
    int output_len = 0;
    int final_len = 0;
    
    // --- 5. Context Setup ---
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CIPHER_CTX_new failed");

    if (1 != EVP_CipherInit_ex(ctx, cipher, NULL, NULL, NULL, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "CipherInit failed");
    }
    
    EVP_CIPHER_CTX_set_key_length(ctx, (int)key_len);
    
    if (1 != EVP_CipherInit_ex(ctx, NULL, NULL, key, iv, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "Key/IV init failed");
    }

    // Padding settings
    // Dla strumieniowych (ChaCha, RC4) i CTR padding jest zawsze 0
    if (EVP_CIPHER_mode(cipher) == EVP_CIPH_STREAM_CIPHER || 
        EVP_CIPHER_mode(cipher) == EVP_CIPH_CTR_MODE) {
        EVP_CIPHER_CTX_set_padding(ctx, 0); 
    }

    // --- 6. Execute Update ---
    if (1 != EVP_CipherUpdate(ctx, output, &output_len, plaintext, (int)plaintext_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherUpdate failed");
    }

    // --- 7. Execute Final (Modified for Attacks) ---
    // Podczas ataków (do_encrypt == 0) padding często jest uszkodzony.
    // Zamiast wyrzucać błąd, ignorujemy końcówkę.
    
    if (1 != EVP_CipherFinal_ex(ctx, output + output_len, &final_len)) {
        if (do_encrypt) {
            // Przy SZYFROWANIU błąd finalizacji to poważny problem.
            EVP_CIPHER_CTX_free(ctx);
            mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherFinal_ex failed during encryption");
        } else {
            // Przy DESZYFROWANIU błąd oznacza zazwyczaj zły padding (skutek ataku).
            // Ignorujemy błąd, czyścimy kolejkę błędów OpenSSL i akceptujemy to, co mamy.
            ERR_clear_error(); 
            final_len = 0; // Odrzucamy ostatni, uszkodzony blok/padding
        }
    }
    
    EVP_CIPHER_CTX_free(ctx);
    
    // --- 8. Resize output ---
    size_t total_len = (size_t)(output_len + final_len);
    mxSetN(plhs[0], total_len);
}