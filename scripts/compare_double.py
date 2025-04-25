import math

def read_double_file(filename):
    """Reads a file containing double values and returns a list of floats."""
    with open(filename, 'r') as f:
        content = f.read().replace('\n', ' ')  # Read entire file and remove newlines

    words = content.split()  # Split into individual words
    data = [float(word) for word in words]  # Convert to float
    print(f"Extracted {len(data)} numbers.")  # Debugging output
    return data

def compare_double_files(file1, file2, tolerance=1e-6):
    """Compares two files containing double values and prints differences."""
    data1 = read_double_file(file1)
    data2 = read_double_file(file2)

    if len(data1) != len(data2):
        print(f"Files have different lengths: {len(data1)} vs {len(data2)}")
        return

    max_error = 0.0
    total_error = 0.0

    for i, (d1, d2) in enumerate(zip(data1, data2)):
        error = abs(d1 - d2)  # Compute absolute difference
        total_error += error
        max_error = max(max_error, error)

    avg_error = total_error / len(data1)

    print(f"Max Error: {max_error:.6e}")
    print(f"Average Error: {avg_error:.6e}")

    # Print mismatched values if they exceed tolerance
    for i, (d1, d2) in enumerate(zip(data1, data2)):
        error = abs(d1 - d2)
        if error > tolerance:
            print(f"Mismatch at index {i}: {d1:.15e} vs {d2:.15e} | Error: {error:.6e}")

if __name__ == "__main__":
    compare_double_files("Code_Output.txt", "Matlab_Output.txt", tolerance=1e-6)
