import matplotlib.pyplot as plt

def read_double_file(filename):
    """Reads a file containing double values and handles INF and NaN values properly."""
    with open(filename, 'r') as f:
        content = f.read().replace('\n', ' ')  # Read entire file and remove newlines

    words = content.split()  # Split into individual words
    data = []
    
    for word in words:
        try:
            # Handle special cases (INF, -INF, NaN)
            if "INF" in word or "inf" in word:
                value = float("inf") if "-" not in word else float("-inf")
            elif "NaN" in word or "#J" in word or "#IND" in word:  # Catching NaN representations
                value = float("0")  # Convert to Python's NaN
                #print(f"Warning: Found NaN in {filename}, replacing with NaN.")
            else:
                value = float(word)  # Normal conversion
            
            data.append(value)

        except ValueError:
            print(f"Skipping invalid entry: {word}")  # Debugging output

    print(f"Extracted {len(data)} valid numbers from {filename}.")
    return data

snr = read_double_file("Output_SNR.txt")
evm_dB = read_double_file("Output_EVM_AGC.txt")
evm_agc_dB = read_double_file("Output_EVM_AGC_DB.txt")
ber = read_double_file("Output_BER.txt")

plt.figure(figsize=(10, 5))
plt.plot(snr, evm_dB, label="EVM_dB", marker="o", linestyle="-")
plt.plot(snr, evm_agc_dB, label="EVM_AGC_dB", marker="s", linestyle="--")
plt.xlabel("SNR")
plt.ylabel("EVM (dB)")
plt.title("EVM_dB and EVM_AGC_dB vs SNR")
plt.grid(True)
plt.legend()
plt.show()

# Plot BER vs SNR
plt.figure(figsize=(10, 5))
plt.semilogy(snr, ber, label="BER", marker="^", linestyle="-", color="red")
plt.xlabel("SNR")
plt.ylabel("BER (log scale)")
plt.title("BER vs SNR")
plt.grid(True)
plt.legend()
plt.show()

print(snr)
print(evm_dB)
print(evm_agc_dB)
print(ber)