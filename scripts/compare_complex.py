import math

def read_complex_file(filename):
    """Reads a file containing complex numbers and returns a list of (real, imag) tuples."""
    with open(filename, 'r') as f:
        content = f.read().replace('\n', ' ')  # Read entire file and remove newlines

    words = content.split()  # Split into individual words (real, sign, imaginary)

    if len(words) % 3 != 0:
        print("Error: File does not contain a multiple of 3 words!")
        return []

    data = []
    for i in range(0, len(words), 3):  # Process in chunks of 3
        real_part = float(words[i])  # First word is real part
        imag_sign = 1 if words[i+1] == '+' else -1  # Second word is sign
        imag_part = float(words[i+2][:-1]) * imag_sign  # Remove 'i' from third word & convert

        data.append((real_part, imag_part))  # Store as tuple (real, imag)

    print(f"Extracted {len(data)} complex numbers.")  # Debugging output
    return data

def complex_abs(real, imag):
    """Computes the absolute value (magnitude) of a complex number."""
    return math.sqrt(real ** 2 + imag ** 2)

def compare_complex_files(file1, file2, tolerance=1e-6):
    """Compares two files containing complex numbers and prints differences."""
    data1 = read_complex_file(file1)
    data2 = read_complex_file(file2)

    if len(data1) != len(data2):
        print(f"Files have different lengths: {len(data1)} vs {len(data2)}")
        return

    max_error = 0.0
    total_error = 0.0

    for i, ((r1, i1), (r2, i2)) in enumerate(zip(data1, data2)):
        error = complex_abs(r1 - r2, i1 - i2)  # Compute absolute error
        total_error += error
        max_error = max(max_error, error)

    avg_error = total_error / len(data1)

    print(f"Max Error: {max_error:.6e}")
    print(f"Average Error: {avg_error:.6e}")

    chk = 0

    # Print mismatched values if they exceed tolerance
    for i, ((r1, i1), (r2, i2)) in enumerate(zip(data1, data2)):
        error = complex_abs(r1 - r2, i1 - i2)
        if error > tolerance:
            print(f"Mismatch at index {i}: ({r1:.15e} + {i1:.15e}i) vs ({r2:.15e} + {i2:.15e}i) | Error: {error:.6e}")
            chk = 1

    if chk == 0:
        print("\nNo Mismatch Found!")
if __name__ == "__main__":
    compare_complex_files("Code_Output.txt", "Matlab_Output.txt", tolerance= 1e-6)
