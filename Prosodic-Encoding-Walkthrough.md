# Introduction

In this file, we discuss the relevant software used to apply the
prosodic encoding algorithm to texts in Section 4.1 of *Context Tree
Classification and Clustering*.

## Required libraries and sourced files

One needs to install the spacyr package and run `spacy_install()` and
`spacy_initialize()`. A fair warning – spacyR may give some problems
downloading, as there is often difficulty linking the package to a
active python environment.

    if (!requireNamespace("spacyr", quietly = TRUE)) install.packages("spacyr"); library(spacyr)
    #Sys.setenv(SPACY_PYTHON = "~/myenvs/myenv/bin") # this is for my own virtual python environment for spacyr
    spacy_install()

    ## Warning in spacy_install(): Skipping installation. Use `force` to force
    ## installation or update. Or use `spacy_download_langmodel()` if you just want to
    ## install a model.

    spacy_initialize()

    ## successfully initialized (spaCy Version: 3.7.6, language model: en_core_web_sm)

    if (!requireNamespace("tokenizers", quietly = TRUE)) install.packages("tokenizers"); library(tokenizers)
    if (!requireNamespace("english", quietly = TRUE)) install.packages("english"); library(english)

    source('DictionaryExtractionV2.R')
    source('findNumSyl.R')

## Web Dictionary web-scraping

A prerequisite to applying the Prosodic encoding algorithm is to know
both the number of syllables and the location of a word’s stressed
syllable for all words in a text. To do this, we implemented a
web-scraping algorithm provided in the file *DictionaryExtractionV2.R*,
where the functions `sAsENG` and `sAsBRIT` (sAs = stress AND syllables)
are used to web-scrape the Merriam Webster and Cambridge dictionaries
respectively. The two functions `sAsENG` and `sAsBRIT` require that the
user specify the part of speech a particular word pertains to. In order
to find the part of speech of a particular word in a sentence, we
leverage the `spaCyr` package that uses NLP models to identify parts of
speeches of tokenized words within sentences.

For illustration, given the word “garage” and its part of speech “noun”,
the “stress AND syllables” can be found with the following call:

    Word = "garage"
    POS = "NOUN"

    eng_SAS = sAsENG(Word = Word, type = POS)

    brit_SAS = sAsBRIT(Word = Word, type = POS)

    #OUTPUT: 

    # [1] "Word: garage |  Entry: garage"
    # [1] "numSyl: 2 |  Stress: 2"
    # [1] "Word: garage |  Entry: garage"
    # [1] "numSyl: 2 |  Stress: 1"

Notice that the word “garage” differs in the pronunciation between
American and UK English.

## Phonological Encoding Function

We now provide the function to perform the phonological encoding. The
function is long; however, we will illustrate with an example in the
chunk after the function initialization. The phonological encoding
process is described with detail in the manuscript *Context Tree
Classification and Clustering*.

    Phonological_Encoding = function(input_file, output_file = "", dictionary_type = c("Webster", "Cambridge"), elements = c("the", "a", "at")) {
      
      # Parse_Token function: Determines the prosodic encoding of a word based on its syllables, stress, and whether it begins a prosodic word
      Parse_Token <- function(Word, Num_Sil, Stressed_syllable, Output_file, Phonetic_word) {       
        prosodyString = ""
        
        # If the word has syllables
        if (Num_Sil > 0) {
          
          # Special case for one-syllable words
          if (Num_Sil == 1) {
            
            # Check if the word is one of the unstressed pronouns or determiners in the elements list
            if (is.element(toupper(Word), elements)) {
              
              # If it is the beginning of a prosodic word, assign 2
              if (Phonetic_word == 1) { 
                prosodyString = "2 "
                Phonetic_word <- 0
              } else { 
                prosodyString = "0 " # Middle of prosodic word, assign 0
              }
              
            } else { 
              
              # Not in the elements list; assign a different encoding
              if (Phonetic_word == 1) { 
                prosodyString = "3 " # Beginning of prosodic word, assign 3
              } else { 
                prosodyString = "1 " # Middle of prosodic word, assign 1
                Phonetic_word = 1
              }
            }
            
          } else {
            # Case for words with more than one syllable
            prosodyString = ""
            
            for(i in 1: Num_Sil) {
              # If the syllable is stressed and starts a prosodic word, assign 3
              if ((Stressed_syllable == i) && (Phonetic_word == 1) && i == 1) {
                prosodyString = paste(prosodyString, "3 ", sep="")
                Phonetic_word = 0
                
                # If the syllable is stressed but not starting a prosodic word, assign 1
              } else if ((Stressed_syllable == i) && (Phonetic_word == 0)) {
                prosodyString = paste(prosodyString, "1 ", sep="")
                
                # If the syllable is not stressed but starts a prosodic word, assign 2
              } else if ((Stressed_syllable != i) && (Phonetic_word == 1)) {
                prosodyString = paste(prosodyString, "2 ", sep="")
                Phonetic_word = 0
                
                # Unstressed syllables within a prosodic word are assigned 0
              } else if ((Stressed_syllable != i) && (Phonetic_word == 0)) {
                prosodyString = paste(prosodyString, "0 ", sep="")
              }
            }
            Phonetic_word = 1
          }    
        }
        
        # Write prosodic encoding to the output file and return values
        cat(prosodyString, file = Output_file, append = TRUE) 
        return(c(Phonetic_word, prosodyString))
      }
      
      ##################################################################################################
      ##################################################################################################
      ######## START HERE: Main body of the Phonological_Encoding function
      ##################################################################################################
      ##################################################################################################
      
      file <- input_file
      Output_file = 0
      
      # If no output file is provided, create a default output file name based on the input file name
      if(output_file == "") {
        Output_file <- paste(substr(file, 1, nchar(file) - 4), "_output.txt", sep = "")
        close(file(Output_file, open = "w")) # Clear the output file
      } else {
        Output_file = output_file
      }
      
      # Create a log file for debugging and tracking processing
      Output_log <- paste(substr(file, 1, nchar(file) - 4), "_log.txt", sep = "")
      
      # Read the entire input file
      text <- readChar(file, file.info(file)$size)
      
      # Tokenize text into sentences
      sentences = tokenize_sentences(text) 
      sink(Output_log)
      
      # Convert numbers to words in the text
      for (i in 1: length(sentences[[1]])) {
        sentence = unlist(tokenize_words(sentences[[1]][i]))
        sentence = gsub(",", "", sentence)
        suppressWarnings(numbers <- sapply(lapply(sentence, as.numeric), is.na))
        indecies = which(!numbers)
        
        # Replace numeric values with their text representation
        if (length(indecies) > 0) {
          for (j in indecies) {
            sentence[j] = toString(as.english(as.numeric(sentence[j])))
          }
        }
        
        # Rebuild sentences with converted numbers
        sentences[[1]][i] = toString(paste(paste(sentence, collapse = ' '), ".", sep = ""))
      }
      
      ################################################################################   
      
      # Process each sentence to compute prosodic encodings
      for (i in 1:length(sentences[[1]])) {
        
        Phonetic_word = 1 # Initialize flag for tracking prosodic words
        tokenize = spacy_parse(sentences[[1]][i], pos = TRUE) # Tokenize sentence into words and POS tags
        
        for (j in 1: dim(tokenize)[1]) { # Loop through each word in the sentence
          word = tokenize$token[j]
          pos = tokenize$pos[j]
          
          # Skip punctuation and certain contractions
          if (pos != "PUNCT" && pos != "SYM" && !is.element(tolower(word), c("'s", "’s", "n't", "n’t", "’t", "'t", "'ll", "’ll", "’d", "'d", "’n", "'n", "’ve", "'ve", "’re", "'re"))) {
            
            sink(Output_log, type = "output", split = TRUE, append = TRUE)
            
            # Default syllable and stress information
            sAs = 0
            
            # Retrieve syllable and stress information based on dictionary type
            if (dictionary_type == "Webster") {
              sAs = sAsENG(word, pos) # U.S. dictionary
            } else {
              sAs = sAsBRIT(word, pos) # British dictionary (if used)
            }
            
            sink()
            
            # If the word is not found in the dictionary, attempt to handle suffixes and prefixes
            if (is.na(sAs[1])) {
              result = WordAppendeges(word) # Analyze suffixes and prefixes
              Num_Sil <- as.numeric(result[3])
              Word <- result[2]
              fullWord = result[1]
              
              # Handle silent 'e' at the end of words
              Word = separateSilentE(Word, eDictionary)
              
              # Count syllables in the word
              for (k in 1: length(Word))
                Num_Sil = Num_Sil + findNumSil(Word[k])
              
              if (Num_Sil < 1)
                Num_Sil = 1
              
              sAs[1] = Num_Sil
              sAs[2] = 1
              
              sink(Output_log, append = TRUE, type = "output", split = TRUE)
              print("Word was not found in dictionary.  Revert to original function")
              print(paste("Word:", result[1]))
              print(paste("numSyl:", Num_Sil,  "|  Stress:", 1))
              sink()
            }
            
          }
          
          # Special case for abbreviations
          if (!is.na(sAs[1]) && sAs[1] == "abb") {
            sink(Output_log, append = TRUE, type = "output", split = TRUE)
            print(paste("Word:", word, "entered as abbreviation"))
            sink()
            
            # Each letter in the abbreviation gets a prosodic encoding of 3
            for (k in 1: nchar(word)) {
              cat("3 ", file = Output_file, append = TRUE)
            }
            
            Phonetic_word = 1 # Reset for next prosodic word
            
            sink(Output_log, type = "output", split = TRUE, append = TRUE)
            print(paste("Prosodic encoding for '", word, "' is: ", gsub(",", "", toString(rep("3 ", nchar(word)))), ".  Phonetic Word= ", Phonetic_word, sep = ""))
            print("--------------------------------------------------------------------")
            sink()
            
          } else if (pos != "PUNCT" && pos != "SYM" && !is.element(tolower(word), c("'s", "’s", "n't", "n’t", "’t", "'t", "'ll", "’s"))) {
            
            # Parse the token and assign prosodic encoding
            prosody = Parse_Token(word, sAs[1], sAs[2], Output_file, Phonetic_word)
            Phonetic_word = prosody[1]
            encoding = prosody[2]
            
            sink(Output_log, type = "output", split = TRUE, append = TRUE)
            print(paste("Prosodic encoding for '", word, "' is: ", encoding, ".  Phonetic Word= ", Phonetic_word, sep = ""))
            print("--------------------------------------------------------------------")
            sink()
            
          } else if (is.element(word, c(".", "?", "!"))) {
            
            # Sentence-ending punctuation gets a special encoding of 4
            cat("4 ", file = Output_file, append = TRUE)
            Phonetic_word = 1 # Reset for next prosodic word
            
            sink(Output_log, type = "output", split = TRUE, append = TRUE)
            print(paste("Prosodic encoding for '", word, "' is: 4.  Phonetic Word= ", Phonetic_word, sep = ""))
            print("--------------------------------------------------------------------")
            sink()
            
          }
        }
      }
    }    

To showcase the function, let us first create a txt file containing two
sentences.

    text = "The man went to the store.  He later went for a walk."
    write(text, file = "sample_text.txt")

Now, let us apply the prosodic encoding function

    suppressWarnings(Phonological_Encoding(input_file = "sample_text.txt", 
                          output_file = "sample_text_prosodic_output.txt", 
                          dictionary_type = "Webster"))

After running the Phonological encoding process, we will see that two
files will be generated – the Output file with the phonological text
encoding, and the output log that includes details from the encoding
process.

Let us look at the output from the two files.

The `sample_text_prosodic_output.txt` file will contain the phonological
encoding of the text.

    suppressWarnings(output_file_text <- readLines("sample_text_prosodic_output.txt"))

    output_file_text

    ## [1] "3 3 3 3 3 3 0 0 4 3 3 0 3 3 3 3 4 "

The `sample_text_log.txt` file will contain the log of the number of
syllables and stress of each word, and the corresponding phonological
encoding.

    output_log_file_text = readLines("sample_text_log.txt")

    output_log_file_text

    ##  [1] "[1] \"Word: the |  Entry: the\""                                             
    ##  [2] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ##  [3] "[1] \"Prosodic encoding for 'the' is: 3 .  Phonetic Word= 1\""               
    ##  [4] "[1] \"--------------------------------------------------------------------\""
    ##  [5] "[1] \"Word: man |  Entry: man\""                                             
    ##  [6] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ##  [7] "[1] \"Prosodic encoding for 'man' is: 3 .  Phonetic Word= 1\""               
    ##  [8] "[1] \"--------------------------------------------------------------------\""
    ##  [9] "[1] \"Calling Cambridge dictionary because phonetic spelling was not given\""
    ## [10] "[1] \"Word: went |  Entry: went\""                                           
    ## [11] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [12] "[1] \"Prosodic encoding for 'went' is: 3 .  Phonetic Word= 1\""              
    ## [13] "[1] \"--------------------------------------------------------------------\""
    ## [14] "[1] \"Word: to |  Entry: to\""                                               
    ## [15] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [16] "[1] \"Prosodic encoding for 'to' is: 3 .  Phonetic Word= 1\""                
    ## [17] "[1] \"--------------------------------------------------------------------\""
    ## [18] "[1] \"Word: the |  Entry: the\""                                             
    ## [19] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [20] "[1] \"Prosodic encoding for 'the' is: 3 .  Phonetic Word= 1\""               
    ## [21] "[1] \"--------------------------------------------------------------------\""
    ## [22] "[1] \"Word: store |  Entry: store\""                                         
    ## [23] "[1] \"numSyl: 3 |  Stress: 1\""                                              
    ## [24] "[1] \"Prosodic encoding for 'store' is: 3 0 0 .  Phonetic Word= 1\""         
    ## [25] "[1] \"--------------------------------------------------------------------\""
    ## [26] "[1] \"Prosodic encoding for '.' is: 4.  Phonetic Word= 1\""                  
    ## [27] "[1] \"--------------------------------------------------------------------\""
    ## [28] "[1] \"Word: he |  Entry: he\""                                               
    ## [29] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [30] "[1] \"Prosodic encoding for 'he' is: 3 .  Phonetic Word= 1\""                
    ## [31] "[1] \"--------------------------------------------------------------------\""
    ## [32] "[1] \"Word: later |  Entry: later\""                                         
    ## [33] "[1] \"numSyl: 2 |  Stress: 1\""                                              
    ## [34] "[1] \"Prosodic encoding for 'later' is: 3 0 .  Phonetic Word= 1\""           
    ## [35] "[1] \"--------------------------------------------------------------------\""
    ## [36] "[1] \"Calling Cambridge dictionary because phonetic spelling was not given\""
    ## [37] "[1] \"Word: went |  Entry: went\""                                           
    ## [38] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [39] "[1] \"Prosodic encoding for 'went' is: 3 .  Phonetic Word= 1\""              
    ## [40] "[1] \"--------------------------------------------------------------------\""
    ## [41] "[1] \"Word: for |  Entry: for\""                                             
    ## [42] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [43] "[1] \"Prosodic encoding for 'for' is: 3 .  Phonetic Word= 1\""               
    ## [44] "[1] \"--------------------------------------------------------------------\""
    ## [45] "[1] \"Word: a |  Entry: a\""                                                 
    ## [46] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [47] "[1] \"Prosodic encoding for 'a' is: 3 .  Phonetic Word= 1\""                 
    ## [48] "[1] \"--------------------------------------------------------------------\""
    ## [49] "[1] \"Word: walk |  Entry: walk\""                                           
    ## [50] "[1] \"numSyl: 1 |  Stress: 1\""                                              
    ## [51] "[1] \"Prosodic encoding for 'walk' is: 3 .  Phonetic Word= 1\""              
    ## [52] "[1] \"--------------------------------------------------------------------\""
    ## [53] "[1] \"Prosodic encoding for '.' is: 4.  Phonetic Word= 1\""                  
    ## [54] "[1] \"--------------------------------------------------------------------\""

## Full disclaimer

The original algorithms developed to extract the phonetic spellings of
words from Merriam-Webster and Cambridge dictionary webpages have been
updated in response to recent modifications to these sites. Initially
implemented in 2018, the original code has since become outdated as the
structure and format of the webpages have evolved over time.
Consequently, the revised algorithms **have not** undergone extensive
validation for web scraping and may be prone to errors. In contrast, the
original version successfully processed hundreds of text documents with
minimal issues.

## Proposal for future prosodic encoding

For future implementations of a prosodic encoding algorithm in text
analysis, it may be beneficial to first explore existing literature for
recurrent neural network architectures capable of predicting both
syllable count and stress patterns, contingent on the language’s country
of origin. Given that we previously utilized spacyR for part-of-speech
tagging, it is plausible that a comparable architectural approach could
be leveraged to extract the desired phonological information, which we
had initially attempted to acquire via web scraping