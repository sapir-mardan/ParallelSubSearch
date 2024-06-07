# ParallelSubSearch

ParallelSubSearch is a Python project designed to illustrate and explore the algorithmic design and parallelization of functions to find the longest contiguous sub-sequence within a sequence. This project includes a performance comparison between serial and parallel execution approaches.

## Description

The project consists of two main parts:
1. **Algorithmic Design**: Implementing functions `longest_contiguous` and `all_longest` to find specific types of sub-sequences within provided data.
2. **Parallel Execution**: Enhancing performance by implementing a parallel approach to process large datasets using Python's `concurrent.futures` module.

## Getting Started

### Dependencies

- Python 3.7 or higher
- pytest for running unit tests

### Installing

- Clone the repository to your local machine:
  ```bash
  git clone https://github.com/sapirmardan/ParallelSubSearch.git

- Install the required packages:
  ```bash
  pip install numpy biopython pytest

### Executing Program
- To run the tests and validate functionality:
  ```bash
  python3 -m pytest --doctest-modules -v
- To execute the program and compare the serial vs parallel performance:
  ```bash
  python3 longest_sub.py

## Implementation Details

### Functions:

- `longest_contiguous_version1(seq, stop_codon="", gap="", case_sensitive=True)`: Calculates the longest contiguous sub-sequence for each unique letter in the sequence. It is case-sensitive by default.
- `longest_contiguous_version2(seq, stop_codon="", gap="", case_sensitive=True)`: Similar to version 1 but optimized for performance using NumPy for large sequences.
- `all_longest`: This function finds all sub-sequences which are longest by length and not necessarily contiguous.
- `run_large_test_serial`: Executes the data processing in a serial manner.
- `run_test_parallel`: Executes the data processing in parallel using `concurrent.futures` to improve performance.

## Files

- `longest_sub.py`: Contains the main logic for sub-sequence identification.
- `test_longest_sub.py`: Contains unit tests for the functions implemented in `longest_sub.py`.

### Example

Here is a quick example of how to use `longest_contiguous`:

```python
from longest_sub import longest_contiguous_version1
sequence = Seq("AACCGGTTAACCGGTT")
result = longest_contiguous_version1(sequence)
print(result)
