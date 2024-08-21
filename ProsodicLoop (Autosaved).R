
ProsodicLoop=function(n){

Parse_Token <- function(Word, Num_Sil,Stressed_syllable,Output_file,Phonetic_word)
{       
	prosodyString= ""
     if (Num_Sil > 0) # fix to >
     {
         	if (Num_Sil == 1)
     	{	#### Check if HER and IT, YOU are object or subject pronoums : also look for This That (No one)
     		if (is.element(toupper(Word), c("ME", "US", "THEM", "HIM", "HER", "YOU", "IT", "A", "AN","OF", "THE"))) #, "FOR", "NOR", "BUT", "OR", "YET", "SO"))) #### One syllable and these are non-stressed pronouns that join main prosodic words
     		#ARE ALL ONE SYLLABLE CONJUNCTIONS PART OF THIS??????
     		{
                if (Phonetic_word == 1) ## This means it is the beginning of a prosodic word, so it is a 2
                { 
                    prosodyString="2 "
                    Phonetic_word <- 0
			    } else ## middle of prosodic word, so a 0
                    prosodyString="0 "     			
     		} 
     		else
     		{
                if (Phonetic_word == 1) ## This means it is the beginning of a prosodic word, so it is a 2
                { 
                   prosodyString= "3 " 
                    #Phonetic_word <- 0. not correct, because we know this word was a prosodic word itself, so move to next
			    } else ## middle of prosodic word, so a 0
                    prosodyString="1 " # for me 
                    Phonetic_word=1
     			
     		}
     		
     	
     	}else 
     	{
     		prosodyString= ""
     		for(i in 1: Num_Sil)
     		{
     			if ((Stressed_syllable == i) && (Phonetic_word == 1) && i==1)
                     {
                     	prosodyString= paste(prosodyString, "3 ", sep="")
                     	Phonetic_word = 0
                     	}
     			else if ((Stressed_syllable == i) && (Phonetic_word == 0))
     				prosodyString= paste(prosodyString, "1 ", sep="")
     			else if ((Stressed_syllable != i) && (Phonetic_word == 1))
     			{
     				prosodyString= paste(prosodyString, "2 ", sep="")
     				Phonetic_word = 0
     			}
				else if ((Stressed_syllable != i) && (Phonetic_word == 0))
     				prosodyString= paste(prosodyString, "0 ", sep="")
     		}
     		Phonetic_word = 1
     	}	
     	
    }
    
    cat(prosodyString ,file=Output_file, append=TRUE) 
    return(c(Phonetic_word, prosodyString))
 }








##################################################################################################
##################################################################################################
######## START HERE
##################################################################################################
##################################################################################################
   source('DictionaryExtractionV2.R')
   source('findNumSyl.R')
   library(spacyr)
   library(tokenizers)
   library(english)
   #file <-  paste("Data/US_Decoding/US_", n, ".txt", sep="")
   #file <-  paste("Data/UK_Decoding/UK_", n, ".txt", sep="")
   file <-  paste("Data/AUS_Decoding/AUS_", n, ".txt", sep="")
   #file <-  paste("Data/Poetry_Decoding/Poetry_", n, ".txt", sep="")

   Output_file <- paste(substr(file,1,nchar(file)-4),"_output.txt",sep="")
   Output_log <- paste(substr(file,1,nchar(file)-4),"_log.txt",sep="")

   close(file(Output_file, open="w")) #clears output file 
   text <- readChar(file, file.info(file)$size)
      
   sentences= tokenize_sentences(text) 
   sink(Output_log)
   #turn numbers to text #########################################################
   for( i in 1: length(sentences[[1]]))
   {
   	sentence= unlist(tokenize_words(sentences[[1]][i]))
   	sentence=gsub(",", "", sentence )
   	numbers= sapply(lapply(sentence, as.numeric), is.na)
   	indecies= which(!numbers)
    if(length(indecies)>0)
    	for (j in indecies)
    	{
    		sentence[j]= toString(as.english(as.numeric(sentence[j])))
    	}
    		
    sentences[[1]][i]= toString(paste(paste(sentence, collapse=' '), ".", sep=""))
   }
   ################################################################################   
      
   for( i in 1:length(sentences[[1]]))  # runs through all sentances of text 
   {
 
   	Phonetic_word=1
    tokenize= spacy_parse(sentences[[1]][i], pos=TRUE)
    for(j in 1: dim(tokenize)[1]) #runs through words in each sentance 
    {
    		word= tokenize$token[j]
    		pos=  tokenize$pos[j]
            
    		if(pos!="PUNCT" && pos!="SYM" && !is.element(tolower(word), c("'s","’s", "n't",  "n’t", "’t", "'t", "'ll", "’ll", "’d", "'d", "’n", "'n", "’ve", "'ve", "’re", "'re" )))
    		{
   				sink(Output_log, type="output",split=TRUE, append=TRUE)
    			sAs= sAsENG(word, pos) #syllables AND stress
    			#sAs= sAsBRIT(word, pos) #uncomment this for BRIT
    		    sink()
    			if(is.na(sAs[1])) #if word not found in dictionary 
    			{
    				result=WordAppendeges(word) #returns vector (Old word, new word, Syllables in prefix and suffix, syllables prefix, syllables suffix)
  					Num_Sil <- as.numeric(result[3])
 					Word<-result[2]
 					fullWord=result[1]
  
  					Word=separateSilentE(Word, eDictionary)
  					
  							for(k in 1: length(Word))
  								Num_Sil= Num_Sil + findNumSil(Word[k])
  								
  						if(Num_Sil<1)
  						Num_Sil=1

  						sAs[1]= Num_Sil
  						sAs[2]=1
  					sink(Output_log,append=TRUE, type="output",split=TRUE) 
  			print("Word was not found in dictionary.  Revert to original function")
   						print(paste("Word:", result[1]))
						print(paste("numSyl:", Num_Sil,  "|  Stress:", 1))

    					 sink()
  				}

    			}
    			if(!is.na(sAs[1]) && sAs[1]=="abb")#if abbreviation, each letter is its own syllable and stress
    			{ 
    				sink(Output_log,append=TRUE, type="output",split=TRUE)
    				print(paste("Word:", word, "entered as abbreviation"))	          
					
					 sink()
				

    		for(k in 1: nchar(word))
    		{
    			cat("3 ",file=Output_file, append=TRUE)
    		}
    		Phonetic_word=1
    		sink(Output_log, type="output",split=TRUE, append=TRUE)
    					print(paste("Prosodic encoding for '", word, "' is: ",gsub(",", "",toString(rep("3 ",nchar(word)))), ".  Phonetic Word= ",Phonetic_word,  sep=""))
    					print("--------------------------------------------------------------------")
    					 sink()
    			}else if(pos!="PUNCT" && pos!="SYM" && !is.element(tolower(word), c("'s","’s", "n't",  "n’t", "’t", "'t", "'ll", "’s" )))
    			{
    				prosody=Parse_Token(word, sAs[1],sAs[2], Output_file, Phonetic_word )
    				Phonetic_word= prosody[1]
    				encoding= prosody[2]
    				sink(Output_log, type="output",split=TRUE, append=TRUE)
     				print(paste("Prosodic encoding for '", word, "' is: ",encoding,".  Phonetic Word= ",Phonetic_word, sep=""))
    				print("--------------------------------------------------------------------")
    				sink()

    		 }else if(is.element(word, c(".", "?", "!")))
    		{
    			 cat("4 ",file=Output_file, append=TRUE)
    			 Phonetic_word=1 
    			 sink(Output_log, type="output",split=TRUE, append=TRUE)
    			 print(paste("Prosodic encoding for '", word, "' is: 4.  Phonetic Word= ",Phonetic_word, sep=""))
    			 print("--------------------------------------------------------------------")
    			 sink()
    				
    		}
    		
    }
  } 	
  }
    
     #make sure correct files and correct dictionaries!!!!
   

   for(i in 9:9){
   	ProsodicLoop(i)
   }
  
    
   






