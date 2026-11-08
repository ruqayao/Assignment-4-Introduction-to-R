# Helper Function: Load FASTA Sequence
load_fasta_sequence <- function(filename) {
  # Read all lines from the file
  lines <- readLines(filename)
  
  # Filter out header lines starting with '>'
  seq_lines <- lines[!grepl("^>", lines)]
  
  # Combine into a single string and convert to uppercase
  full_seq <- paste(seq_lines, collapse = "")
  return(toupper(full_seq))
}

# Load the sequence from Part 1
file_path <- "chr1_GL383518v1_alt.fa"
dna_seq <- load_fasta_sequence(file_path)

# Part 1: Sequence Indexing
cat(" PART 1: Sequence Indexing \n")

# The 10th letter is at index 10
tenth_letter <- substr(dna_seq, 10, 10)
cat("10th letter of the sequence:", tenth_letter, "\n")

# The 758th letter is at index 758
letter_758 <- substr(dna_seq, 758, 758)
cat("758th letter of the sequence:", letter_758, "\n\n")

# Part 2: Reverse Complement
cat(" PART 2: Reverse Complement \n")

# Split sequence into individual characters
seq_chars <- strsplit(dna_seq, NULL)[[1]]

# Reverse the sequence order
rev_chars <- rev(seq_chars)

# Translate bases using Watson-Crick-Franklin pairing rules
complement_base <- function(base) {
  if (base == "A") {
    return("T")
  } else if (base == "T") {
    return("A")
  } else if (base == "C") {
    return("G")
  } else if (base == "G") {
    return("C")
  } else {
    return(base)
  }
}

rev_comp_chars <- sapply(rev_chars, complement_base)
rev_comp_seq <- paste(rev_comp_chars, collapse = "")

# 79th letter
letter_79 <- substr(rev_comp_seq, 79, 79)
cat("79th letter of reverse complement:", letter_79, "\n")

# 500th through 800th letters
slice_500_800 <- substr(rev_comp_seq, 500, 800)
cat("500th through 800th letters:", slice_500_800, "\n\n")

# Part 3: Kilobase Nucleotide Frequency List
cat(" PART 3: Kilobase List Creation \n")

seq_len <- nchar(dna_seq)
kb_list <- list()

# Loop through the sequence in steps of 1000 base pairs
start_pos <- seq(1, seq_len, by = 1000)

for (i in seq_along(start_pos)) {
  s_start <- start_pos[i]
  s_end <- min(s_start + 999, seq_len)
  segment <- substr(dna_seq, s_start, s_end)
  
  # Count individual nucleotides
  seg_chars <- strsplit(segment, NULL)[[1]]
  count_a <- sum(seg_chars == "A")
  count_c <- sum(seg_chars == "C")
  count_g <- sum(seg_chars == "G")
  count_t <- sum(seg_chars == "T")
  
  counts_vector <- c(A = count_a, C = count_c, G = count_g, T = count_t)
  
  # Use starting position as the list key name
  key_name <- paste("kb", s_start, sep = "_")
  kb_list[[key_name]] <- counts_vector
}

cat("Finished creating kilobase frequency list.\n\n")

# Part 4: Data Frame Construction and Row Sum Analysis
cat(" PART 4: Data Frame & Analysis \n")

# Part 4a & 4b: Combine all list vectors into a single data frame
all_kb_df <- do.call(rbind, kb_list)

# Part 4c: Calculate the sum of each row
row_sums <- rowSums(all_kb_df)

# Display the first few row sums to verify
print(head(row_sums))
print(tail(row_sums))

# ANSWERS TO CONCEPTUAL QUESTIONS:
# Question 1: What is the expected sum for each list?
# Answer: The expected sum for each row (list) is 1000, because each kilobase 
#         block is defined as containing 1000 base pairs.

# Question 2: Are there any lists whose sums are not equal to the expected value?
# Answer: Yes. The very last row (representing the trailing sequence segment) 
#         has a sum less than 1000 because the total sequence length is not an 
#         exact multiple of 1000.

# Question 3: Provide a general explanation for the differences in your expected
#             results and your observed results.
# Answer: Genomic sequences rarely divide evenly into exact 1000 base pair blocks. 
#         After partitioning full 1000-bp windows, any remaining nucleotides 
#         form a shorter final block, causing its nucleotide sum to fall below 1000.