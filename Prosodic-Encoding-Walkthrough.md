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
syllable for all words in a text. To do this, we originally implemented
a web-scraping algorithm provided in the file
*DictionaryExtractionV2.R*, where the functions `sAsENG` and `sAsBRIT`
(sAs = stress AND syllables) were used to web-scrape the Merriam Webster
and Cambridge dictionaries respectively. However, due to recent changes
in the web format of the Merriam Webster dictionary, I have augmented
the function `sAsENG` to instead extract the United States English
pronunciation from the Cambridge dictionary (which is provided). The two
functions `sAsENG` and `sAsBRIT` require that the user specify the part
of speech a particular word pertains to. In order to find the part of
speech of a particular word in a sentence, we leverage the `spaCyr`
package that uses NLP models to identify parts of speeches of tokenized
words within sentences.

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
     Parse_Token <- function(Word, Num_Sil, Stressed_syllable, Output_file, Phonetic_word)
    {       
        # Initialize an empty string to hold the prosodic encoding for the word
        prosodyString = ""
        
        # Proceed only if the word has one or more syllables
        if (Num_Sil > 0) # Ensure Num_Sil is positive
        {
            # Special case for words with exactly one syllable
            if (Num_Sil == 1)
            {   
                #### Check if the word is a non-stressed pronoun or determiner ####
                if (is.element(toupper(Word), c("ME", "US", "THEM", "HIM", "HER", "YOU", "IT", "A", "AN", "OF", "THE", "FOR"))) 
                {
                    # Check if this is the beginning of a prosodic word
                    if (Phonetic_word == 1) 
                    { 
                        # Assign prosody code 2 (indicates start of a prosodic word)
                        prosodyString = "2 "
                        Phonetic_word <- 0 # Reset Phonetic_word flag
                    } 
                    else 
                    {
                        # Assign prosody code 0 (indicates middle of a prosodic word)
                        prosodyString = "0 "                
                    } 
                } 
                # Case for one-syllable word that is not a non-stressed pronoun/determiner
                else
                {
                    # Check if this is the beginning of a prosodic word
                    if (Phonetic_word == 1) 
                    { 
                        # Assign prosody code 3 (indicates stressed prosodic word)
                        prosodyString = "3 " 
                    } 
                    else 
                    {
                        # Assign prosody code 1 (stressed syllable in middle of prosodic word)
                        prosodyString = "1 " 
                        Phonetic_word = 1 # Update Phonetic_word to indicate this word was prosodic
                    }   
                }
            } 
            # Case for multi-syllable words
            else 
            {
                prosodyString = "" # Initialize/clear the prosody string
                # Loop through each syllable in the word
                for (i in 1: Num_Sil)
                {
                    # If current syllable is stressed and it's the first syllable in the prosodic word
                    if ((Stressed_syllable == i) && (Phonetic_word == 1) && i == 1)
                    {
                        # Assign prosody code 3 (stressed first syllable in prosodic word)
                        prosodyString = paste(prosodyString, "3 ", sep = "")
                        Phonetic_word = 0 # Reset Phonetic_word flag
                    }
                    # If current syllable is stressed but not the first syllable in the prosodic word
                    else if ((Stressed_syllable == i) && (Phonetic_word == 0))
                    {
                        # Assign prosody code 1 (stressed syllable in middle of prosodic word)
                        prosodyString = paste(prosodyString, "1 ", sep = "")
                    }
                    # If current syllable is not stressed and it's the first syllable in the prosodic word
                    else if ((Stressed_syllable != i) && (Phonetic_word == 1))
                    {
                        # Assign prosody code 2 (unstressed first syllable in prosodic word)
                        prosodyString = paste(prosodyString, "2 ", sep = "")
                        Phonetic_word = 0 # Reset Phonetic_word flag
                    }
                    # If current syllable is not stressed and not the first syllable in the prosodic word
                    else if ((Stressed_syllable != i) && (Phonetic_word == 0))
                    {
                        # Assign prosody code 0 (unstressed syllable in middle of prosodic word)
                        prosodyString = paste(prosodyString, "0 ", sep = "")
                    }
                }
                Phonetic_word = 1 # Reset Phonetic_word flag for the next word
            }   
        }
        
        # Write the prosody string to the output file and append it
        cat(prosodyString, file = Output_file, append = TRUE) 
        
        # Return the updated Phonetic_word flag and the prosody string
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

    text = "Police are hunting for answers on what caused an explosion that injured twenty nine people in New York's Chelsea neighborhood shortly before a second suspicious device was found nearby."
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

    ## [1] "2 1 3 3 0 2 1 0 3 3 3 2 0 1 0 3 3 0 3 0 3 3 0 3 3 3 3 0 3 0 0 3 0 2 1 2 1 0 2 1 0 2 1 3 3 2 1 4 "

The `sample_text_log.txt` file will contain the log of the number of
syllables and stress of each word, and the corresponding phonological
encoding.

    output_log_file_text = readLines("sample_text_log.txt")

    output_log_file_text

    ##   [1] "[1] \"Word: police |  Entry: police\""                                                                                          
    ##   [2] "[1] \"numSyl: 2 |  Stress: 2\""                                                                                                 
    ##   [3] "[1] \"Prosodic encoding for 'police' is: 2 1 .  Phonetic Word= 1\""                                                             
    ##   [4] "[1] \"--------------------------------------------------------------------\""                                                   
    ##   [5] "[1] \"Word: are |  Entry: are\""                                                                                                
    ##   [6] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##   [7] "[1] \"Prosodic encoding for 'are' is: 3 .  Phonetic Word= 1\""                                                                  
    ##   [8] "[1] \"--------------------------------------------------------------------\""                                                   
    ##   [9] "[1] \"Word: hunting |  Entry: hunting\""                                                                                        
    ##  [10] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [11] "[1] \"Prosodic encoding for 'hunting' is: 3 0 .  Phonetic Word= 1\""                                                            
    ##  [12] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [13] "[1] \"Word: for |  Entry: for\""                                                                                                
    ##  [14] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [15] "[1] \"Prosodic encoding for 'for' is: 2 .  Phonetic Word= 0\""                                                                  
    ##  [16] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [17] "[1] \"ENTRY AND WORD DO NOT MATCH!!!!! \""                                                                                      
    ##  [18] "[1] \"Word: answers |  Entry: answer\""                                                                                         
    ##  [19] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [20] "[1] \"Prosodic encoding for 'answers' is: 1 0 .  Phonetic Word= 1\""                                                            
    ##  [21] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [22] "[1] \"Word: on |  Entry: on\""                                                                                                  
    ##  [23] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [24] "[1] \"Prosodic encoding for 'on' is: 3 .  Phonetic Word= 1\""                                                                   
    ##  [25] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [26] "[1] \"Word: what |  Entry: what\""                                                                                              
    ##  [27] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [28] "[1] \"Prosodic encoding for 'what' is: 3 .  Phonetic Word= 1\""                                                                 
    ##  [29] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [30] "[1] \"Word: caused |  Entry: caused\""                                                                                          
    ##  [31] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [32] "[1] \"Prosodic encoding for 'caused' is: 3 .  Phonetic Word= 1\""                                                               
    ##  [33] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [34] "[1] \"Word: an |  Entry: an\""                                                                                                  
    ##  [35] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [36] "[1] \"Prosodic encoding for 'an' is: 2 .  Phonetic Word= 0\""                                                                   
    ##  [37] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [38] "[1] \"Word: explosion |  Entry: explosion\""                                                                                    
    ##  [39] "[1] \"numSyl: 3 |  Stress: 2\""                                                                                                 
    ##  [40] "[1] \"Prosodic encoding for 'explosion' is: 0 1 0 .  Phonetic Word= 1\""                                                        
    ##  [41] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [42] "[1] \"Word: that |  Entry: that\""                                                                                              
    ##  [43] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [44] "[1] \"Prosodic encoding for 'that' is: 3 .  Phonetic Word= 1\""                                                                 
    ##  [45] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [46] "[1] \"Word: injured |  Entry: injured\""                                                                                        
    ##  [47] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [48] "[1] \"Prosodic encoding for 'injured' is: 3 0 .  Phonetic Word= 1\""                                                            
    ##  [49] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [50] "[1] \"Word: twenty |  Entry: twenty\""                                                                                          
    ##  [51] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [52] "[1] \"Prosodic encoding for 'twenty' is: 3 0 .  Phonetic Word= 1\""                                                             
    ##  [53] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [54] "[1] \"Word: nine |  Entry: nine\""                                                                                              
    ##  [55] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [56] "[1] \"Prosodic encoding for 'nine' is: 3 .  Phonetic Word= 1\""                                                                 
    ##  [57] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [58] "[1] \"Word: people |  Entry: people\""                                                                                          
    ##  [59] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [60] "[1] \"Prosodic encoding for 'people' is: 3 0 .  Phonetic Word= 1\""                                                             
    ##  [61] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [62] "[1] \"Word: in |  Entry: in\""                                                                                                  
    ##  [63] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [64] "[1] \"Prosodic encoding for 'in' is: 3 .  Phonetic Word= 1\""                                                                   
    ##  [65] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [66] "[1] \"Word: new |  Entry: new\""                                                                                                
    ##  [67] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [68] "[1] \"Prosodic encoding for 'new' is: 3 .  Phonetic Word= 1\""                                                                  
    ##  [69] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [70] "[1] \"Word: york |  Entry: York\""                                                                                              
    ##  [71] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [72] "[1] \"Prosodic encoding for 'york' is: 3 .  Phonetic Word= 1\""                                                                 
    ##  [73] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [74] "[1] \"Calling Merriam Webster dictionary because phonetic spelling was not given or  dictionary entry may have been incorrect\""
    ##  [75] "[1] \"Entry not found in backup Webster dictionary\""                                                                           
    ##  [76] "[1] \"Word was not found in dictionary.  Revert to original function\""                                                         
    ##  [77] "[1] \"Word: chelsea\""                                                                                                          
    ##  [78] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [79] "[1] \"Prosodic encoding for 'chelsea' is: 3 0 .  Phonetic Word= 1\""                                                            
    ##  [80] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [81] "[1] \"Word: neighborhood |  Entry: neighborhood\""                                                                              
    ##  [82] "[1] \"numSyl: 3 |  Stress: 1\""                                                                                                 
    ##  [83] "[1] \"Prosodic encoding for 'neighborhood' is: 3 0 0 .  Phonetic Word= 1\""                                                     
    ##  [84] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [85] "[1] \"Word: shortly |  Entry: shortly\""                                                                                        
    ##  [86] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [87] "[1] \"Prosodic encoding for 'shortly' is: 3 0 .  Phonetic Word= 1\""                                                            
    ##  [88] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [89] "[1] \"Word: before |  Entry: before\""                                                                                          
    ##  [90] "[1] \"numSyl: 2 |  Stress: 2\""                                                                                                 
    ##  [91] "[1] \"Prosodic encoding for 'before' is: 2 1 .  Phonetic Word= 1\""                                                             
    ##  [92] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [93] "[1] \"Word: a |  Entry: a\""                                                                                                    
    ##  [94] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ##  [95] "[1] \"Prosodic encoding for 'a' is: 2 .  Phonetic Word= 0\""                                                                    
    ##  [96] "[1] \"--------------------------------------------------------------------\""                                                   
    ##  [97] "[1] \"Word: second |  Entry: second\""                                                                                          
    ##  [98] "[1] \"numSyl: 2 |  Stress: 1\""                                                                                                 
    ##  [99] "[1] \"Prosodic encoding for 'second' is: 1 0 .  Phonetic Word= 1\""                                                             
    ## [100] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [101] "[1] \"Word: suspicious |  Entry: suspicious\""                                                                                  
    ## [102] "[1] \"numSyl: 3 |  Stress: 2\""                                                                                                 
    ## [103] "[1] \"Prosodic encoding for 'suspicious' is: 2 1 0 .  Phonetic Word= 1\""                                                       
    ## [104] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [105] "[1] \"Word: device |  Entry: device\""                                                                                          
    ## [106] "[1] \"numSyl: 2 |  Stress: 2\""                                                                                                 
    ## [107] "[1] \"Prosodic encoding for 'device' is: 2 1 .  Phonetic Word= 1\""                                                             
    ## [108] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [109] "[1] \"Word: was |  Entry: was\""                                                                                                
    ## [110] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ## [111] "[1] \"Prosodic encoding for 'was' is: 3 .  Phonetic Word= 1\""                                                                  
    ## [112] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [113] "[1] \"Word: found |  Entry: found\""                                                                                            
    ## [114] "[1] \"numSyl: 1 |  Stress: 1\""                                                                                                 
    ## [115] "[1] \"Prosodic encoding for 'found' is: 3 .  Phonetic Word= 1\""                                                                
    ## [116] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [117] "[1] \"Word: nearby |  Entry: nearby\""                                                                                          
    ## [118] "[1] \"numSyl: 2 |  Stress: 2\""                                                                                                 
    ## [119] "[1] \"Prosodic encoding for 'nearby' is: 2 1 .  Phonetic Word= 1\""                                                             
    ## [120] "[1] \"--------------------------------------------------------------------\""                                                   
    ## [121] "[1] \"Prosodic encoding for '.' is: 4.  Phonetic Word= 1\""                                                                     
    ## [122] "[1] \"--------------------------------------------------------------------\""

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
