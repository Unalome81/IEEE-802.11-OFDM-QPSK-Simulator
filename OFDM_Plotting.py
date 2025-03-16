import matplotlib.pyplot as plt

def read_values_from_file(filename):
    with open(filename, "r") as file:
        values = [float(line.strip()) for line in file if line.strip()]
    return values

snr = read_values_from_file("Output_SNR.txt")
evm_dB = read_values_from_file("Output_EVM_AGC.txt")
evm_agc_dB = read_values_from_file("Output_EVM_AGC_DB.txt")
ber = read_values_from_file("Output_BER.txt")

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