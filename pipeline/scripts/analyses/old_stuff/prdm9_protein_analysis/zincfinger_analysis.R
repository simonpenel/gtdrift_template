#options(repos = c(CRAN = "https://cloud.r-project.org/"))
#.libPaths("~/R/libs")
#install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(Biostrings)
library(stringr)

# zincfinger_analysis.R
args <- commandArgs(trailingOnly = TRUE)
protein_file <- args[1]
output_file <- args[2]


if (!file.exists(protein_file)) {
  stop(paste("Error: The input file does not exist:", protein_file))
}


print(paste("Path to the input:", protein_file))


# This code contains three functions. In the overall, the input is a FASTA file and the output is a data frame
# which contains, for each aminoacidic sequence, its ID, its total number of ZF (*Num_ZF*), the diversity index
# of the positions in contact with the DNA from the ZFs wich are in the longest tandem array (*ZFDa*), 
# the diversity index of the positions in contact with the DNA from all ZF in the sequence (*ZFDg*),
# the number of ZFs of 28 AA (*Num_ZF_28*), the number of ZFs of 29 AA (*Num_ZF_29*)
# and the number ZF that are in the longest tandem array (*Max_Tandem_Length*)

calc.ZFD = function(aa.seq, allele.name, equation = 2, tandem.only = TRUE, min.ZFs = 4, return.prop = TRUE) {
  # equation = 1 --> true heterozygosity (sampling without replacement)
  # equation = 2 --> approximation for heterozygosity in large populations (pi)
  # tandem.only --> when TRUE will only use the longest set of tandem ZFs, when FALSE it considers all ZFs
  # min.ZFs --> function will return NA if there are not enough ZFs (or enough tandem ZFs if tandem.only is true)
  
  
  # Patterns to identify ZFs with 28 and 29 AAs
  c2h2.pattern.28 = paste(c(rep("[a-z]", 2), "c", rep("[a-z]", 2), "c", rep("[a-z]", 12), "h", rep("[a-z]", 3), "h", rep("[a-z]", 5)), collapse = "")
  c2h2.pattern.29 = paste(c(rep("[a-z]", 2), "c", rep("[a-z]", 2), "c", rep("[a-z]", 12), "h", rep("[a-z]", 4), "h", rep("[a-z]", 5)), collapse = "")
  # For X2-CXXC-X12-HXXXH-X5 ZnF syntax, the DNA binding residue positions on the alpha helix are 12,14,15,18.
  # For X7-CXXC-X12-HXXXH ZnF syntax, the DNA binding residue positions on the alpha helix are 17,19,20,23.
  
  
  aa.seq = tolower(aa.seq)
  
  # Extract ZFs of length 28 AA and 29 AA
  ZFs.28 = str_extract_all(aa.seq, c2h2.pattern.28)[[1]]
  ZFs.29 = str_extract_all(aa.seq, c2h2.pattern.29)[[1]]
  
  # Remove the 20th position in any 29 AA ZF to convert it to 28 AA
  if (length(ZFs.29) > 0) {
    ZFs.29 = sapply(ZFs.29, function(zf) {
      paste0(substr(zf, 1, 19), substr(zf, 21, 29))  # Remove the 20th AA
    })
  }
  
  # Count ZFs
  total.ZFs.28 = length(ZFs.28)
  total.ZFs.29 = length(ZFs.29)
  total.ZFs = total.ZFs.28 + total.ZFs.29
  
  # If the total ZFs are less than min.ZFs, return a message and stop processing
  if (total.ZFs < min.ZFs) {
    message("There are fewer than four copies of ZF in the sequence.")
    return(list(result.ZFDa = NA, result.ZFDg = NA, ZF.sequences = NULL, 
                num.ZFs = total.ZFs, num.ZFs.28 = total.ZFs.28, num.ZFs.29 = total.ZFs.29, 
                max.tandem.length = NA))
  }
  
  # Combine both sets of ZFs and handle tandem ZFs
  max.tandem.length = 0  # Initialize the longest tandem length
  if (tandem.only) {
    ZF.pos.28 = str_locate_all(aa.seq, c2h2.pattern.28)[[1]]
    ZF.pos.29 = str_locate_all(aa.seq, c2h2.pattern.29)[[1]]
    
    # Ensure that the ZF.pos arrays have values before processing
    if (nrow(ZF.pos.28) > 1) {
      is.tandem.28 = sapply(ZF.pos.28[, 1], function(x) {
        (x + 28) %in% ZF.pos.28[, 1] || (x + 27) %in% ZF.pos.28[, 1] || (x + 29) %in% ZF.pos.28[, 1]
      })
      if (sum(is.tandem.28) > 0) {
        tandem_data = process_tandem_ZFs(ZF.pos.28, 28)
        max.tandem.length = max(tandem_data)  # Update the maximum tandem length
      }
    }
    
    if (nrow(ZF.pos.29) > 1) {
      is.tandem.29 = sapply(ZF.pos.29[, 1], function(x) {
        (x + 29) %in% ZF.pos.29[, 1] || (x + 28) %in% ZF.pos.29[, 1] || (x + 30) %in% ZF.pos.29[, 1]
      })
      if (sum(is.tandem.29) > 0) {
        tandem_data = process_tandem_ZFs(ZF.pos.29, 29)
        max.tandem.length = max(max.tandem.length, max(tandem_data))  # Update the maximum tandem length
      }
    }
  }
  
  # Output the count of original ZFs (28 and 29)
  message(paste("Original ZFs: ", total.ZFs.28, " 28 AA ZFs and ", total.ZFs.29, " 29 AA ZFs (converted to 28 AA)"))
  
  # Calculate ZFDa
  ZFs = c(ZFs.28, ZFs.29)
  ZF.length = 28 # Start with the default value of 28 AA
  
  if (length(ZFs.28) >= min.ZFs & length(ZFs.29) >= min.ZFs) {
    message("Both 28 AA and 29 AA ZFs were found. Handling insertion as a variation.")
    ZF.length = 29
  } else if (length(ZFs.28) >= min.ZFs) {
    ZFs = ZFs.28
    ZF.length = 28
  } else if (length(ZFs.29) >= min.ZFs) {
    ZFs = ZFs.29
    ZF.length = 28  # After the modification, ZFs.29 are treated as 28 AA
  }
  
  # Calculate ZFDa (formerly pi)
  ZFs.ZFDa = sapply(1:ZF.length, FUN = function(x) {
    1 - sum((table(substring(ZFs, x, x)) / length(ZFs))^2)
  })
  
  # Calculate ZFDg using all ZFs (both 28 and 29 AA)
  all_ZFs = c(ZFs.28, ZFs.29)  # Use all ZFs for ZFDg
  ZFs.ZFDg = sapply(1:28, FUN = function(x) {
    1 - sum((table(substring(all_ZFs, x, x)) / length(all_ZFs))^2)
  })
  
  # Calculate the value of ZFDg
  ZFDg = sum(ZFs.ZFDg[c(12, 14, 15, 18)]) / sum(ZFs.ZFDg)
  
  # Return the number of ZFs, their sequences, and the calculated ZFDs
  return(list(result.ZFDa = sum(ZFs.ZFDa[c(12, 14, 15, 18)]) / sum(ZFs.ZFDa), 
              result.ZFDg = ZFDg, 
              ZF.sequences = ZFs, 
              num.ZFs = total.ZFs, 
              num.ZFs.28 = total.ZFs.28, 
              num.ZFs.29 = total.ZFs.29,
              max.tandem.length = max.tandem.length)) # Add maximum tandem length
}

# Auxiliary function to process tandem ZFs
process_tandem_ZFs = function(ZF.pos, ZF.length) {
  num.ZFs = numeric(nrow(ZF.pos))
  group = 1
  num.ZFs[group] = 1  # Start the count
  
  for (i in 2:nrow(ZF.pos)) {
    if ((ZF.pos[i, 1] - ZF.length) == ZF.pos[i - 1, 1] ||  # ZF in tandem
        (ZF.pos[i, 1] - 1) == ZF.pos[i - 1, 1] ||  # Starts one amino acid before
        (ZF.pos[i, 1] + 1) == ZF.pos[i - 1, 1]) {  # Starts one amino acid after
      num.ZFs[group] = num.ZFs[group] + 1
    } else {  # If the ZF is not in tandem with the previous one
      group = group + 1
      num.ZFs[group] = 1  # Reset the count for the new group
    }
  }
  
  return(num.ZFs)
}

# Function to manage a FASTA as an imput and to generate a data.frame as an output
calc.ZFD.FASTA = function(fasta_file, return.prop = TRUE) {
  # Read the FASTA file
  fasta_seqs <- readAAStringSet(fasta_file) 
  
  # Verifi if the FASTA file is empty
  if (length(fasta_seqs) == 0) {
    return(list(resultados = data.frame(ID = "# 0 candidate PRDM9 proteins"), num_sequences = 0))
  }
  
  # Obtain the number of sequences of the FASTA
  num_sequences <- length(fasta_seqs)
  
  # Generate an empty data frame to store the results
  resultados <- data.frame(
    ID = character(),
    Num_ZF = integer(),
    ZFDa = numeric(),
    ZFDg = numeric(),
    Num_ZF_28 = integer(),
    Num_ZF_29 = integer(),
    Max_Tandem_Length = integer(),
    stringsAsFactors = FALSE
  )
  
  # Iteration over each sequence of the FASTA file
  for (i in 1:length(fasta_seqs)) {
    seq_id <- names(fasta_seqs)[i]  # Obtener el ID de la secuencia
    aa.seq <- tolower(as.character(fasta_seqs[[i]]))
    
    # Call to the function calc.ZFD
    calc_results <- calc.ZFD(aa.seq, seq_id, return.prop = return.prop)
    
    # Auxiliary code to verify if calc.ZFD has returned valid results
    if (is.null(calc_results) || !all(c("num.ZFs", "result.ZFDa", "result.ZFDg") %in% names(calc_results))) {
      warning(paste("Results could not be obtained for the sequence:", seq_id))
      
      resultados <- rbind(resultados, data.frame(
        ID = seq_id,
        Num_ZF = 0,
        ZFDa = NA,
        ZFDg = NA,
        Num_ZF_28 = 0,
        Num_ZF_29 = 0,
        Max_Tandem_Length = NA,
        stringsAsFactors = FALSE
      ))
      next  
    }
    
    # Extract the results of the function calc.ZFD
    total_zf <- calc_results$num.ZFs
    result_zfda <- calc_results$result.ZFDa
    result_zfdg <- calc_results$result.ZFDg
    max_tandem_length <- calc_results$max.tandem.length 
    
    # Add the results to the data frame
    resultados <- rbind(resultados, data.frame(
      ID = seq_id,
      Num_ZF = total_zf,
      ZFDa = result_zfda,
      ZFDg = result_zfdg,
      Num_ZF_28 = calc_results$num.ZFs.28,
      Num_ZF_29 = calc_results$num.ZFs.29,
      Max_Tandem_Length = max_tandem_length,
      stringsAsFactors = FALSE
    ))
  }
  
  # Return the results and the number of sequences
  return(list(resultados = resultados, num_sequences = num_sequences))
}

# Input FASTA file
fasta_file <- protein_file

# Run the analisis and save the results in a .csv file
ZFD_analysis <- calc.ZFD.FASTA(fasta_file)

output_conn <- file(output_file, open = "w")

# Add the header line with the number of sequences
if (ZFD_analysis$num_sequences > 0) {
  cat(paste("#", ZFD_analysis$num_sequences, "candidate PRDM9 proteins\n"), file = output_conn)
} else {
  # If the FASTA file is empty, write "0 candidate PRDM9 proteins"
  cat("# 0 candidate PRDM9 proteins\n", file = output_conn)
}


if (ZFD_analysis$num_sequences > 0) {
  write.csv(ZFD_analysis$resultados, file = output_conn, row.names = FALSE)
}

close(output_conn)

