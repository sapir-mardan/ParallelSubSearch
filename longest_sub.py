import random
import time
from collections import defaultdict
import numpy as np
import os

from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor

#############################################
# explore one function structure, version 1
#############################################

def longest_contiguous_version1(seq, stop_codon="", gap="", case_sensitive=True): 
    """
    Count the length of the longest contiguous subsequences.

    This function finds each contiguous subsequence in ``seq``,
    e.g. ``aaaa`` or ``bb`` and counts its length and returns the
    letter with the longest contiguous subsequence and its length.

    It returns a dictionary with the key being the letter and the
    value being the length of the longest subsequence. If two letters
    have longest subsequences of the same length, both are included
    separately in the dictionary.

    Each unique letter in ``seq`` is counted separately so the
    sequence "aabbcc" will find the longest subsequences of each of
    "a", "b" and "c".

    If ``seq`` is empty then an empty dicitonary is returned.

    Args:
        seq (Bio.Seq.Seq): A sequence of characters.
        stop_codon (str): A single character symbol for stop codon in sequence, usually "*". Default is sequences contains no stop codons, empty str, e.g. stop_codon="". 
        gap (str):  A single character symbol for gap in sequence, usually "-". Default is sequences contains no gaps, empty str, e.g. gap="". 
        case_sensitive (bool): Lowercase and uppercase letters are treated as distinct characters if True (default in BioPython's Seq objects are case sensitive).
                                ***   Explanation: ***
                                Lowercase and uppercase letters are treated as distinct characters. This means that "A" and "a" are considered different characters when comparing sequences.
                                While the letters may differ in case, their representation in the sequence indicates the same nucleotide or amino acid. 
                                For example, both "A" and "a" represent adenine.


    Returns:
        dict: The count of the longest subsequences, keyed by letter.

    Examples:
        >>> longest_contiguous("aacbbb")
        {'b': 3}
        >>> longest_contiguous("aabbbaabbc")
        {'b': 3}
        >>> longest_contiguous("aaabbbaabb")
        {'a': 3, 'b': 3}
        >>> longest_contiguous("aAaAbb", case_sensitive=True)
        {'b': 2}
        >>> longest_contiguous("aAaAbb", case_sensitive=False)
        {'a': 4}
    """
    
    # check if seq is from Bio.Seq.Seq type
    if not isinstance(seq, Seq):
        raise ValueError("seq is not an instance of Seq (Bio.Seq.Seq)")
    
    # check for empty seq
    if len(seq) == 0:
        return {}
    
    # change to str and all lower case
    seq = str(seq.lower())

    # check stop_codon or gap are empty or a single character 
    if len(stop_codon) > 1 or len(gap) > 1:
        raise ValueError("stop_codon and gap must be either empty or single characters")
    
    # check for non alphabetic characters
    # check if no symbols other than stop codon,or gap.
    if not seq.replace(stop_codon, "").replace(gap, "").isalpha():
        raise ValueError(f"Sequence must contain only letters, '{stop_codon}', or '{gap}' ")
    

    
    # dictionary for counting lengths
    longest_counts = defaultdict(int)
    current_letter = None
    current_count = 0
    
    # Loop through each letter in the sequence
    for letter in seq:
        # add counts
        if letter == current_letter:
            current_count += 1
        else:
            # # If a new letter sequence starts, update longest count if current is greater, and reset tracking.
            if current_count > longest_counts[current_letter]:
                longest_counts[current_letter] = current_count
            current_letter = letter
            current_count = 1
    
    # update the last sequence
    if current_count > longest_counts[current_letter]:
        longest_counts[current_letter] = current_count

    longest_counts.pop(None, None)

    # Find the maximum count to filter the dictionary
    if longest_counts:  # Ensure the dictionary is not empty
        max_count = max(longest_counts.values())
        longest_counts = {k: v for k, v in longest_counts.items() if v == max_count}

    return longest_counts


################################################
# explore another function structure, version 2
################################################

def longest_contiguous_version2(seq, stop_codon="", gap="", case_sensitive=True):
    """
    Count the length of the longest contiguous subsequences.

    This function finds each contiguous subsequence in ``seq``,
    e.g. ``aaaa`` or ``bb`` and counts its length and returns the
    letter with the longest contiguous subsequence and its length.

    It returns a dictionary with the key being the letter and the
    value being the length of the longest subsequence. If two letters
    have longest subsequences of the same length, both are included
    separately in the dictionary.

    Each unique letter in ``seq`` is counted separately so the
    sequence "aabbcc" will find the longest subsequences of each of
    "a", "b" and "c".

    If ``seq`` is empty then an empty dicitonary is returned.

    Args:
        seq (Bio.Seq.Seq): A sequence of characters.
        stop_codon (str): A single character symbol for stop codon in sequence, usually "*". Default is sequences contains no stop codons, empty str, e.g. stop_codon="". 
        gap (str):  A single character symbol for gap in sequence, usually "-". Default is sequences contains no gaps, empty str, e.g. gap="". 
        case_sensitive (bool): Lowercase and uppercase letters are treated as distinct characters if True (default in BioPython's Seq objects are case sensitive).
                                ***   Explanation: ***
                                Lowercase and uppercase letters are treated as distinct characters. This means that "A" and "a" are considered different characters when comparing sequences.
                                While the letters may differ in case, their representation in the sequence indicates the same nucleotide or amino acid. 
                                For example, both "A" and "a" represent adenine.


    Returns:
        dict: The count of the longest subsequences, keyed by letter.

    Examples:
        >>> longest_contiguous("aacbbb")
        {'b': 3}
        >>> longest_contiguous("aabbbaabbc")
        {'b': 3}
        >>> longest_contiguous("aaabbbaabb")
        {'a': 3, 'b': 3}
        >>> longest_contiguous("aAaAbb", case_sensitive=True)
        {'b': 2}
        >>> longest_contiguous("aAaAbb", case_sensitive=False)
        {'a': 4}
    """
    
    # check if seq is from Bio.Seq.Seq type
    if not isinstance(seq, Seq):
        raise ValueError(f"Sequence must contain only letters, '{stop_codon}', or '{gap}' ")
    
    # check if empty
    if len(seq) == 0:
        return {}

    # change to str and if not case sensitive to all lower case
    seq = str(seq.lower()) if not case_sensitive else str(seq)

    # check stop_codon or gap are empty or a single character 
    if len(stop_codon) > 1 or len(gap) > 1:
        raise ValueError("stop_codon and gap must be either empty or single characters")
    
    # check for non alphabetic characters
    # check if no symbols other than stop codon,or gap.
    if not seq.replace(stop_codon, "").replace(gap, "").isalpha():
        raise ValueError(f"Sequence must contain only letters, '{stop_codon}', or '{gap}")
    

    # convert sequence into a numpy array of ASCII values for efficient computation.
    # The dtype is exactly of size byte matching ascii, ensuring that no conversion or scaling is necessary.
    # np.frombuffer is used to avoid copying, and creating instead almost instantaneously an array view to minimize memory allocation
    seq_arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    # create a boolean array marking where consecutive letters change
    different_letter_arr = np.r_[True, seq_arr[:-1] != seq_arr[1:]]
    # identify indices in the sequence where a new character run starts.
    seq_start_indices_arr = np.flatnonzero(different_letter_arr)
    # calculate lengths of each contiguous sequence.
    seq_lengths_arr = np.diff(np.r_[seq_start_indices_arr, len(seq_arr)])
    # find the longest sequences
    max_length = np.max(seq_lengths_arr)
    # Extract the ASCII values of letters with the longest contiguous sequences.
    letters = seq_arr[seq_start_indices_arr[seq_lengths_arr == max_length]]

    return {chr(letter): max_length for letter in set(letters)}


#########################################
#   define the final longest_contiguous()
#########################################

def select_longest_contiguous(longest_contiguous_version1, longest_contiguous_version2):
    """
    Returns the faster version out of two versions of longest_contiguous().
    Faster execution times is optimized for efficiency, particularly when dealing with large datasets, such as numerous sequences.
    """
    
    # lambda function to create a single sequence
    create_seq = lambda: Seq("".join(random.choices(["g", "a", "t", "c"], k=10000000)))

    seq = create_seq()

    # time the duration of first function
    start = time.perf_counter()  # start the performance timer
    answer1 = longest_contiguous_version1(seq)  # Run the analysis serially
    duration_first_function = time.perf_counter() - start   # stop the performance timer

    # time the duration of second function
    start = time.perf_counter()  # start the performance timer
    answer2 = longest_contiguous_version2(seq)  # Run the analysis serially
    duration_second_function = time.perf_counter() - start   # stop the performance timer

    # check which function is more efficient 
    if duration_first_function < duration_second_function:
        version = "version 1"
    else:
        version = "version 2"
        
    # print the timing and chosen function
    print("Timing two versions of function 'longest_contiguous()'")
    print(f"  First function: {duration_first_function:.2f} seconds - {answer1}")
    print(f"  Second function: {duration_second_function:.2f} seconds - {answer2}")
    print(f"{version} of 'longest_contiguous' functions is optimized for efficiency, resulting in faster execution times particularly when dealing with large datasets, such as numerous sequences.\n")

# assigning  longest_contiguous as version 2
longest_contiguous = longest_contiguous_version2

def all_longest(each_longest):
    """
    Find the longest subsequences across multiple sequences.

    This function combines the outputs from ``longest_contiguous``
    to a single dictionary. For each letter in each of the input
    dictionaries, the output will contain that letter as a key and
    the largest of all associated values as its value.

    Args:
        each_longest (list): A list of longest subsequence counts.

    Returns:
        dict: The largest values for each letter from the inputs.

    Examples:
        >>> all_longest([{"a": 4}, {"b": 6}])
        {'a': 4, 'b': 6}
        >>> all_longest([{"a": 4}, {"b": 6, "a": 2}])
        {'a': 4, 'b': 6}
        >>> all_longest([{"a": 4}, {"b": 10, "a": 2}])
        {'a': 4, 'b': 10}
        >>> all_longest([{"a": 4}, {"b": 6, "a": 2}, {"b": 10}])
        {'a': 4, 'b': 10}
    """
    #create empty dict
    all_longest_dict = defaultdict(int)
    for d in each_longest:
        for k, v in d.items():
            # Update the stored value for each character, if the current value is greater than the stored value.
            all_longest_dict[k] = max(all_longest_dict[k], v)
    return dict(all_longest_dict)


def generate_random_test_data(num_seq=100_000, seq_length=1000):
    """
    Generate a large amount of example data to test the above functions.

    Args:
        num_seq (int): The number of sequnces to return.
        seq_length (int): The length of each individual sequence.
    """
    sequences = []
    for i in range(num_seq):
        sequences.append(Seq("".join(random.choices(["g", "a", "t", "c"], k=seq_length))))
    return sequences

def run_large_test_serial(sequences):
    """
    Run through a list of provided sequences in series and return the result.
    """
    longest_each_sequence = []
    for seq in sequences:
        longest_each_sequence.append(longest_contiguous(seq))

    answer = all_longest(longest_each_sequence)

    return answer


def run_large_test_parallel(sequences):
    """
    Run through a list of provided sequences in parallel and return the result.
    """
    with ProcessPoolExecutor() as pool:
        # Calculate optimal chunk size for dividing work among processes based on available CPU cores
        CHUNKSIZE = len(sequences) // os.cpu_count()
        # Parallel execution
        each_longest = pool.map(longest_contiguous, sequences, chunksize=CHUNKSIZE)
    return all_longest(each_longest)


if __name__ == "__main__":
    
    #  show the final longest_contiguous() as the fasted function from 2 versions:
    select_longest_contiguous(longest_contiguous_version1, longest_contiguous_version2)

    print("Generating test squences")
    sequences = generate_random_test_data()

    print("Running serial test")
    start = time.perf_counter()  # start the performance timer
    answer1 = run_large_test_serial(sequences)  # Run the analysis serially
    duration_seconds = time.perf_counter() - start  # stop the performance timer
    print(f"  Serial: {duration_seconds:.2f} seconds - {answer1}")

    print("Running parallel test")
    start = time.perf_counter()
    answer2 = run_large_test_parallel(sequences)  # Run the parallel analysis
    duration_seconds = time.perf_counter() - start
    print(f"Parallel: {duration_seconds:.2f} seconds - {answer2}")
    
    print("\nGenerating larger number of test sequences")
    #generate more sequences to show run_large_test_parallel is faster when dealing with bigger data
    sequences = generate_random_test_data(num_seq=300_000, seq_length=1000)
    start = time.perf_counter()  # start the performance timer
    answer1 = run_large_test_serial(sequences)  # Run the analysis serially
    duration_seconds = time.perf_counter() - start  # stop the performance timer
    print(f"  Serial: {duration_seconds:.2f} seconds - {answer1}")

    print("Running parallel test")
    start = time.perf_counter()
    answer2 = run_large_test_parallel(sequences)  # Run the parallel analysis
    duration_seconds = time.perf_counter() - start
    print(f"Parallel: {duration_seconds:.2f} seconds - {answer2}")
