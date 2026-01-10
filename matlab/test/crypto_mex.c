#include "mex.h"
#include <openssl/evp.h>
#include <openssl/provider.h> 
#include <string.h> 

/*
 * Universal crypto wrapper for MATLAB
 * Version 3 (Robust): Loads "legacy provider" for OpenSSL 3.0+
 */

// --- Global flag to load providers only once ---
static int openssl_initialized = 0;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // --- STEP 0: One-time OpenSSL initialization (for Blowfish, RC4, etc.) ---
    // This executes only once, on the first call to this MEX function
    if (!openssl_initialized) {
        // Explicitly load "legacy" (for Blowfish, RC4, DES)
        // and "default" (for AES, ChaCha20).
        // In OpenSSL < 3.0 these functions don't exist, but the code 
        // will still run (they simply won't load anything).
        // In OpenSSL 3.0+ this is REQUIRED.
        OSSL_PROVIDER_load(NULL, "legacy");
        OSSL_PROVIDER_load(NULL, "default");
        
        openssl_initialized = 1; 
        // Set flag even if loading failed to avoid retrying endlessly.
    }

    // --- 1. Input validation ---
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CRYPTO:nrhs", "Wymagane 5 argumentów: (plaintext, key, iv, cipher_name, do_encrypt)");
    }
    
    if (!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]) || !mxIsUint8(prhs[2])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "plaintext, key i iv muszą być typu uint8.");
    }
    if (!mxIsChar(prhs[3])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "Nazwa szyfru (cipher_name) musi być łańcuchem znaków (string).");
    }
    if (!mxIsLogicalScalar(prhs[4])) {
        mexErrMsgIdAndTxt("CRYPTO:type", "do_encrypt musi być wartością logiczną (true/false).");
    }

    // --- 2. Get data from MATLAB ---
    const unsigned char* plaintext = (const unsigned char*)mxGetData(prhs[0]);
    size_t plaintext_len = mxGetNumberOfElements(prhs[0]);
    const unsigned char* key = (const unsigned char*)mxGetData(prhs[1]);
    size_t key_len = mxGetNumberOfElements(prhs[1]);
    const unsigned char* iv = (const unsigned char*)mxGetData(prhs[2]);
    size_t iv_len = mxGetNumberOfElements(prhs[2]);
    char cipher_name[64];
    mxGetString(prhs[3], cipher_name, sizeof(cipher_name));
    int do_encrypt = (int)mxIsLogicalScalarTrue(prhs[4]);
    
    // --- 3. Select OpenSSL cipher ---
    const EVP_CIPHER* cipher = EVP_get_cipherbyname(cipher_name);
    
    if (cipher == NULL) {
        // If algorithm is not found, this is the main suspect
        mexErrMsgIdAndTxt("CRYPTO:cipher", "Nieznany szyfr: %s. (Jeśli to Blowfish/RC4, upewnij się, że masz OpenSSL z 'legacy provider').", cipher_name);
    }

    // --- 4. Validate Key and IV ---
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

    // --- 5. Allocate output memory ---
    size_t max_output_len = plaintext_len + EVP_CIPHER_block_size(cipher);
    plhs[0] = mxCreateNumericMatrix(1, max_output_len, mxUINT8_CLASS, mxREAL);
    unsigned char* output = (unsigned char*)mxGetData(plhs[0]);
    int output_len = 0;
    int final_len = 0;
    
    // --- 6. Execute operation ---
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (ctx == NULL) mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CIPHER_CTX_new failed");

    // Step 1: Initialize cipher type AND mode (enc/dec)
    if (1 != EVP_CipherInit_ex(ctx, cipher, NULL, NULL, NULL, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        // Error here usually means provider issue (fixed in Step 0)
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherInit_ex (Krok 1: set cipher) failed. Prawdopodobnie błąd providera OpenSSL.");
    }

    // Step 2: EXPLICITLY set key length
    if (1 != EVP_CIPHER_CTX_set_key_length(ctx, (int)key_len)) {
         EVP_CIPHER_CTX_free(ctx);
         mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CIPHER_CTX_set_key_length failed");
    }
    
    // Step 3: Set Key and IV
    if (1 != EVP_CipherInit_ex(ctx, NULL, NULL, key, iv, do_encrypt)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherInit_ex (Krok 3: set key/iv) failed");
    }

    // Padding setup
    if (EVP_CIPHER_mode(cipher) == EVP_CIPH_STREAM_CIPHER || 
        EVP_CIPHER_mode(cipher) == EVP_CIPH_CTR_MODE) {
        EVP_CIPHER_CTX_set_padding(ctx, 0); 
    }

    // Processing
    if (1 != EVP_CipherUpdate(ctx, output, &output_len, plaintext, (int)plaintext_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherUpdate failed");
    }
    if (1 != EVP_CipherFinal_ex(ctx, output + output_len, &final_len)) {
        EVP_CIPHER_CTX_free(ctx);
        mexErrMsgIdAndTxt("CRYPTO:openssl", "EVP_CipherFinal_ex failed (sprawdź klucz/padding)");
    }
    EVP_CIPHER_CTX_free(ctx);
    
    // --- 7. Set final output size ---
    size_t total_len = (size_t)(output_len + final_len);
    mxSetN(plhs[0], total_len);
}
