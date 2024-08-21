
#http://teacher.scholastic.com/reading/bestpractices/vocabulary/pdf/prefixes_suffixes.pdf
#https://dictionary.cambridge.org/us/grammar/british-grammar/word-formation/prefixes
#https://www.thoughtco.com/common-suffixes-in-english-1692725


setClass("prefix",slots=list(prefix="character", numSil="numeric" ))
p1<- new("prefix",prefix="ANTI", numSil=2)
#p2<- new("prefix",prefix="DE", numSil=1) #Dean doesnt work 
p3<- new("prefix",prefix="DIS", numSil=1)
p4<- new("prefix",prefix="UN", numSil=1)
p5<- new("prefix",prefix="EN", numSil=1)
p6<- new("prefix",prefix="EM", numSil=1)
p7<- new("prefix",prefix="FORE", numSil=1)
p8<- new("prefix",prefix="IN", numSil=1)
p9<- new("prefix",prefix="IM", numSil=1)
p10<- new("prefix",prefix="IL", numSil=1)
p11<- new("prefix",prefix="IR", numSil=1)
p12<- new("prefix",prefix="INTER", numSil=2)
p13<- new("prefix",prefix="MID", numSil=1)
p14<- new("prefix",prefix="MIS", numSil=1)
p15<- new("prefix",prefix="NON", numSil=1)
p16<- new("prefix",prefix="OVER", numSil=2)
p17<- new("prefix",prefix="PRE", numSil=1)
#p18<- new("prefix",prefix="RE", numSil=1)
p19<- new("prefix",prefix="SEMI", numSil=2)
p20<- new("prefix",prefix="SUB", numSil=1)
p21<- new("prefix",prefix="SUPER", numSil=2)
p22<- new("prefix",prefix="TRANS", numSil=1)
p23<- new("prefix",prefix="UN", numSil=1)
p24<- new("prefix",prefix="UNDER", numSil=2)
p25<- new("prefix",prefix="AUTO", numSil=2)
p26<- new("prefix",prefix="DOWN", numSil=1)
p27<- new("prefix",prefix="EXTRA", numSil=2)
p28<- new("prefix",prefix="MEGA", numSil=2)
p29<- new("prefix",prefix="OUT", numSil=1)
p30<- new("prefix",prefix="POST", numSil=1)
p31<- new("prefix",prefix="PRE", numSil=1)
p32<- new("prefix",prefix="PRO", numSil=1)
p33<- new("prefix",prefix="TELE", numSil=2)
p34<- new("prefix",prefix="ULTRA", numSil=2)
p35<- new("prefix",prefix="UP", numSil=1)

prefixArray<-c(p1,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31, p32,p33,p34,p35)


####################################################################
setClass("suffix",slots=list(suffix="character", numSil="numeric" ))
s1<- new("suffix",suffix="ACY", numSil=2)
s2<- new("suffix",suffix="AL", numSil=1)
s3<- new("suffix",suffix="ANCE", numSil=1)
s4<- new("suffix",suffix="IED", numSil=1)
s5<- new("suffix",suffix="DOM", numSil=1)
#s6<- new("suffix",suffix="ED", numSil=0)
s7<- new("suffix",suffix="FORE", numSil=1)
s8<- new("suffix",suffix="OR", numSil=1)
s9<- new("suffix",suffix="ISM", numSil=1)
s10<- new("suffix",suffix="IST", numSil=1)
s11<- new("suffix",suffix="ITY", numSil=2)
s12<- new("suffix",suffix="TY", numSil=1)
s13<- new("suffix",suffix="MENT", numSil=1)
s14<- new("suffix",suffix="NESS", numSil=1)
s15<- new("suffix",suffix="SHIP", numSil=1)
s16<- new("suffix",suffix="SION", numSil=1)
s17<- new("suffix",suffix="TION", numSil=1)
s18<- new("suffix",suffix="ATE", numSil=1)
s19<- new("suffix",suffix="EN", numSil=1)
s20<- new("suffix",suffix="IFY", numSil=2)
s21<- new("suffix",suffix="FY", numSil=1)
s22<- new("suffix",suffix="IZE", numSil=1)
s23<- new("suffix",suffix="ISE", numSil=1)
s24<- new("suffix",suffix="ABLE", numSil=2)
s25<- new("suffix",suffix="IBLE", numSil=2)
s26<- new("suffix",suffix="AL", numSil=1)
s27<- new("suffix",suffix="ESQUE", numSil=2)
s28<- new("suffix",suffix="FUL", numSil=1)
s29<- new("suffix",suffix="IC", numSil=1)
s30<- new("suffix",suffix="ICAL", numSil=2)
s31<- new("suffix",suffix="OUS", numSil=1) #Problem with something like ious, which would be another suffix where sometimes the I is part of 2 syllables like nutritious vs studious.  Note, may have this problem with some other prefix/suffix //RESOLVE: CAN ONLY HAPPEN WITH DOUBLE VOWEL!!!!!
s32<- new("suffix",suffix="ISH", numSil=1)
s33<- new("suffix",suffix="IVE", numSil=1)
s34<- new("suffix",suffix="LESS", numSil=1)
s35<- new("suffix",suffix="EY", numSil=1)
s36<- new("suffix",suffix="EN", numSil=1)
s37<- new("suffix",suffix="ER", numSil=1)
s38<- new("suffix",suffix="EST", numSil=1)
s39<- new("suffix",suffix="IC", numSil=1)
s40<- new("suffix",suffix="LY", numSil=1)
#s41<- new("suffix",suffix="S", numSil=0) ## dont know for sure 
#s42<- new("suffix",suffix="ES", numSil=1) ## dont know for sure
#s43<- new("suffix",suffix="EST", numSil=1)
#s44<- new("suffix",suffix="IC", numSil=1)
s45<- new("suffix",suffix="ING", numSil=1)

suffixArray<-c(s1,s2,s3,s4,s5,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30,s31, s32,s33,s34,s35, s36,s37,s38,s39,s40, s45)



#how we will reference: prefixArray[[1]]@prefix to get the prefix 

WordAppendeges <- function(Word)
{
# we will go through word, find prefixes and suffixes that appear in the word, remove the prefixes and suffixes from the word, return the old word, the new word without the prefixes and suffixes, and the number of syllables that pertain to removes prefixes/suffixes.  We will return a vector with ("word", "oldWord", numSil)

oldWord<- Word 		#save word string 
newWord<-Word 		#initalize for later
prefix<-""  	        #initalize for later
suffix<-""             #initialize for later
prefixIndex<-0
suffixIndex<-0
numSil<-0

# find longest prefix

 for (ind in 1 : length(prefixArray))
   {
   	prefixLength= nchar(prefixArray[[ind]]@prefix) #length of prefix
   	currentPrefix<-prefixArray[[ind]]@prefix #save prefix string
   	if(toupper(substr(Word,1,prefixLength))==currentPrefix) #if prefix is in word
   	  {
   	  	if(nchar(currentPrefix)>nchar(prefix)) #if prefix is longer than last found prefix ensures we find longest prefix, "un" 	     										#vs "under"
   	  	{
   	  	prefix<-currentPrefix
   	  	prefixIndex<-ind
   	  	}
   	  }
   }
   
   # find longest suffix

 for (ind in 1 : length(suffixArray))
   {
   	suffixLength= nchar(suffixArray[[ind]]@suffix) #length of suffix
   	currentSuffix<-suffixArray[[ind]]@suffix #save suffix string
   	if(toupper(substr(Word,nchar(Word)-suffixLength +1,nchar(Word)))==currentSuffix) #if suffix is in word
   	  {
   	  	if(nchar(currentSuffix)>nchar(suffix)) #if suffix is longer than last found prefix ensures we find longest suffix 	     										      
   	  	{
   	  	suffix<-currentSuffix
   	  	suffixIndex<-ind
   	  	}
   	  }
   }
   
newWord<- substr(Word,nchar(prefix)+1,nchar(Word)-nchar(suffix))
prefixSyl=0
suffixSyl=0
if(prefixIndex!=0)
{
	prefixSyl=prefixArray[[prefixIndex]]@numSil
	numSil<-numSil+prefixArray[[prefixIndex]]@numSil
}

if(suffixIndex!=0)
{
	suffixSyl=suffixArray[[suffixIndex]]@numSil
	numSil<-numSil+suffixArray[[suffixIndex]]@numSil
}
prefixSylSuffixSyl=c(prefixSyl, suffixSyl)

results=c(oldWord,newWord,numSil,prefixSylSuffixSyl)
return(results)

}

## ------------ Program Silaba --------------- ##
##
## Last Update: 
## 
##
## 
## 

### CODE:
## 2 = non-streesed syllable - beggining of prosodic word
## 0 = non-stressed syllable - not beginning of prosodic word
## 1 = stressed syllable - not beginning  of prosodic word
## 3 = stressed syllable - beginning of prosodic word
## 4 = end of sentence = period








	# 1.	Separate prefixes and suffixes from root words.
# examples:  pre-view, work-ing, re-do, end-less, & out-ing

# 2. Are two (or more) consonants next to each other?
# Divide between the 1st and 2nd consonants.
# examples:  buf-fet, des-sert, ob-ject, ber-ry, & pil-grim

# Never split 2 consonants that make only 1 sound when pronounced together and aren't the same letter (i.e., 'ff').
# examples:  th, sh, ph, th, ch, & wh

# 3. Is the consonant surrounded by vowels?
# Does the vowel have a long sound?  (Like the 'i' in line)
# Divide before the consonant.
# examples:  ba-by, re-sult, i-vy, fro-zen, & Cu-pid
# Does the vowel have a short sound?  (Like the 'i' in mill)
# Divide after the consonant.
# examples:  met-al, riv-er, mod-el, val-ue, & rav-age

# 4. Does the word end with 'ckle'?
# Divide right before the 'le.'
# examples:  tack-le, freck-le, tick-le, & buck-le

# 5. Does the word end with 'le' (not 'ckle')?
# Is the letter before the 'le' a consonant?
# Divide 1 letter before the 'le.'
# examples:  ap-ple, rum-ble, fa-ble, & ta-ble
# Is the letter before the 'le' a vowel?
# Do nothing.
# examples:  ale, scale, sale, file, & tile



#https://en.wiktionary.org/wiki/Category:English_words_by_number_of_syllables use this to test later

#http://www.businessenglishresources.com/learn-english-for-business/teachers-section/mini-lessons/pronunciation-lessons-pronunciation-changes-in-words-that-are-both-nouns-and-verbs/.  syllable stress on a word that can act as a verb or noun is different, consider "the SUSpect" and "I susPECT you" could be the presence of "the" before to indicate noun.
#Point out the pattern — the stress goes on the first syllable if it is a noun and the second syllable if it is a verb.
#need to figure out how we will account for words ending in a silent e that are a non-ending subword, likeliness, suspenseful, etc.  could hardcode a few 
#http://www.syllablecount.com/syllable/rules/       syllable counting algorithm 
 
 
#Idea: check every sequence of three consecutive vowels, to see if two of the three vowels are emphasized (will account normally for only 1 vowel sound)
## ------------ Parse_Token --------------- ##
##
## Description:
## 
## 
##
##Begin copy paste at eDictionary

findNumSil<- function(Word)
{
	 Vowels <- c("a","e","i","o","u","y","A","E","I","O","U","Y")

	result=WordAppendeges(Word) #returns vector (Old word, new word, Syllables in prefix/suffix)
  Num_Sil <- as.numeric(result[3])
  Word<-result[2]
  fullWord=result[1]
  
  

   ## Count the number of Vowelss
     for (Ind in 1 : nchar(Word))
        if (is.element(substr(Word,Ind,Ind),Vowels)) 
            Num_Sil <- Num_Sil + 1
       
   ## subtract doulble vowels and Y+vowel, ...
     for (Ind in 1:nchar(Word)-1)       ###### ditongo. 
        if (is.element(toupper(substr(Word,Ind,Ind+1)), c("OA","OI","OO", "OY", "OU", "AU", "UA", "UE", "UO", "UI", "EE","EI","EU", "EY", "YA","YE","YO","YI","YU","AI","AY")))  #INCLUDED AY FOR SATURDAY, ESSAYWRITER, PLAY, ETC.
              Num_Sil <- Num_Sil - 1 
              
     # add ui back for counterintuitively, but not fruit
      for (Ind in 1:nchar(Word)) 
      {
               if (is.element(toupper(substr(Word,Ind,Ind+1)), c("UI")) && is.element(toupper(substr(Word,Ind-1,Ind-1)), c("T")))
               	 Num_Sil <- Num_Sil + 1  	
               if (is.element(toupper(substr(Word,Ind,Ind+2)), c("UIN")) && !is.element(toupper(substr(Word,Ind-1,Ind-1)), c("Q", "G"))) ## look for more guin or quin of UIN as one sillable
               	 Num_Sil <- Num_Sil + 1  	
      
               	 
       if (is.element(toupper(substr(Word,Ind,Ind+2)), c("UIC")) && is.element(toupper(substr(Word,Ind-1,Ind-1)), c("S"))) ## Suicide, may have other words, double vowel UI
               	 Num_Sil <- Num_Sil + 1  	
      
      }
              
              
              
      ##SINCE I ADDED "OA" ABOVE, NEED TO FIGURE OUT WHY OASIS SOUNDS THE WAY IT DOES, OASIS IS HARDCODED  
 for (Ind in 1 : nchar(Word)-4)
        if (is.element(toupper(substr(Word,Ind,Ind+4)),c("OASIS","OASES")))
            Num_Sil <- Num_Sil + 1             
              
    #add back one syllable for "EYO", "AYE" SPECIAL CASE WITH DIPTHONG Beyond, subtracted because double subtraction on three consecutive vowels 
      for (Ind in 1:nchar(Word)) 
               if (is.element(toupper(substr(Word,Ind,Ind+2)), c("EYO","AYE")))
               	 Num_Sil <- Num_Sil + 1  
               	 
  for (Ind in 1:nchar(Word))  #Subtract a syllable in dipthong "IA" in certain cases especially vs phobia
                if (is.element(toupper(substr(Word,Ind,Ind+3)), c("CIAL","TIAL","GIAT","TIAN","RIAG")))
               	 Num_Sil <- Num_Sil -1 
	
               
             
   ##Subtract one syllable if word ends with "ay"
      #if (is.element(toupper(substr(Word,nchar(Word)-1,nchar(Word))), c("AY")))
          #  Num_Sil <- Num_Sil - 1 
          #PUT AY IN DIPTHONG SECTION
        
     ### "EA" sounds         
     ##Check IDEA or exceptions  
     ## if it ends in (d or r or rh, k, v, n, c) + ea - we should separate the syllable = +0
     ## if it ends in (any other letter)+ea = -1
     #nausea vs sea?
     for (Ind in 1:(nchar(Word)))   
      {
        if (is.element(toupper(substr(Word,Ind,Ind+1)), c("EA")))
         {
           if (!is.element(toupper(substr(Word,Ind-1,Ind-1)), c("D", "R", "K", "V", "N", "C"))){
           			if (!is.element(toupper(substr(Word,Ind-2,Ind-1)), c("RH")))
           			{
              				Num_Sil <- Num_Sil - 1
              				}
              			}
 	
           #SEAN, DEAN CASE WHERE WORD ENDS WITH N
           if((Ind+2<=nchar(Word)) && toupper(substr(Word,Ind+2,Ind+2))=="N" && is.element(toupper(substr(Word,Ind-1,Ind-1)), c("D", "R", "K", "V", "N", "C")))
              	Num_Sil <- Num_Sil - 1
         }
       }
       
              
      for (Ind in 1:(nchar(Word)-3))  #great, can be fixed in upper for loop, possibly in revised code version
     {   
      if(is.element(toupper(substr(Word,Ind,Ind+3)), c("REAT")))
       	{
       		Num_Sil <- Num_Sil - 1
       		
          }
     }
     
      for (Ind in 1:nchar(Word))
      {   #IF READY IS IN THE WORD, Num_Sil <- Num_Sil - 1     ALREADY, READYMADE, BREADY	
        if (is.element(toupper(substr(Word,Ind,Ind+3)), c("READ"))) 
             Num_Sil <- Num_Sil - 1 
         if (is.element(toupper(substr(Word,Ind,Ind+2)), c("AYA")))  #Layaway, note: words like betrayal, portrayal, get filtered through "al" suffix,
             Num_Sil <- Num_Sil + 1 
      }



	 for(Ind in 1:nchar(Word)-3)
	 if (is.element(toupper(substr(Word,Ind,Ind+3)), c("LIKE")) && Ind+3!=nchar(Word)) #if the word contains "like" and the like does not occur at the end of the word, where the silent "e" condition fails to subtract a syllable 
              Num_Sil <- Num_Sil - 1
              
              
      #this is commented because of new "IE" section with science, etc below         
     #for (Ind in 1:nchar(Word))       
       # if (is.element(toupper(substr(Word,Ind,Ind+2)), c("IEU"))) # Lieutenant for example
        #     Num_Sil <- Num_Sil - 1 ## -1 because "EU" was already subtracted 
        
   ## subtract "tion", "sion", ...
      for (Ind in 1:nchar(Word))       
        if (is.element(toupper(substr(Word,Ind,Ind+3)), c("TION", "SION")))
            Num_Sil <- Num_Sil - 1
      

### DEAL with this in the silent e case in the middle of words      
#      for (Ind in 1:nchar(Word))  #facebook, silent e in the middle of words not accounted for 
#      {      
#        if (is.element(toupper(substr(Word,Ind,Ind+2)), c("ACE"))&& ((Ind+2)!=nchar(Word)))
#            Num_Sil <- Num_Sil - 1
#      }


    ## subtract ? Need to check all these words: Patience vs Science, conscience, hierarchical, trierarchies, furrie, carrier, piercing, tier, inconscient, susie, frier
	if (!toupper(substr(Word,1,5)) == "SCIEN") ## If it starts with "SCIEN" then we wont subtract "IE"
	{
      for (Ind in 1:(nchar(Word)-2)) ## the loop only goes to nchar(Word)-2 because ending in ie will go to slient e case  
        if (is.element(toupper(substr(Word,Ind,Ind+1)), c("IE")))
        {     
            Num_Sil <- Num_Sil - 1   
        
          if (is.element(toupper(substr(Word,Ind-1,Ind+2)), c("RIER", "RIER")) && !is.element(toupper(substr(Word,Ind-2,Ind+2)), c("TRIER", "BRIER")))   ### add back +1
              Num_Sil <- Num_Sil + 1   

          if (is.element(toupper(substr(Word,Ind-2,Ind+2)), c("SCIEN")) && !is.element(toupper(substr(Word,Ind-3,Ind+2)), c("NSCIEN")))   ### add back + 1
              Num_Sil <- Num_Sil + 1   
        }
    }
          
   ## subtract silent "e"
    if (toupper(substr(Word,nchar(Word),nchar(Word))) == "E" && toupper(substr(Word,nchar(Word)-1,nchar(Word)-1)) != "E")   #DONT SUBTRACT FOR DOUBLE E COMITTEE,PEE, GUARENTEE check for ending dipthongs
             Num_Sil <- Num_Sil - 1
              
     if ((toupper(Word) == "SHE") || (toupper(Word) == "HE") || (toupper(Word) == "BE")||(toupper(Word) == "WE"))   ## these are syllables themselves
              Num_Sil <- Num_Sil + 1
     if (toupper(substr(Word,nchar(Word)-1,nchar(Word))) == "LE")   ## able, uncle, ... not silent e
              Num_Sil <- Num_Sil + 1
     if (toupper(substr(Word,nchar(Word)-2,nchar(Word))) == "ILE")   ## automobile, ... is silent e
              Num_Sil <- Num_Sil - 1
     if (is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))), c("THE", "PHE", "SUE","LUE", "DUE", "RUE", "NUE")))   ## the, catastrophe, epitome, anemone ... not silent e
              Num_Sil <- Num_Sil + 1
      if (is.element(toupper(substr(Word,nchar(Word)-3,nchar(Word))), c("TOME", "MONE")))   ## epitome, anemone ... not silent e
              Num_Sil <- Num_Sil + 1

   
   ## subtract silent "ed" (THERE ARE EXCEPTIONS! NEED TO CHECK!)
     if (toupper(substr(Word,nchar(Word)-1,nchar(Word))) == "ED")   
              Num_Sil <- Num_Sil - 1
     if (toupper(Word) == "LED")    ## this is syllable itself
              Num_Sil <- Num_Sil + 1
     if (toupper(substr(Word,nchar(Word)-2,nchar(Word))) == "TED")   # interested
       Num_Sil <- Num_Sil + 1
     if (toupper(substr(Word,nchar(Word)-2,nchar(Word))) == "DED")   # decided
       Num_Sil <- Num_Sil + 1
     
  ## subtract "mes" (THERE ARE EXCEPTIONS! NEED TO CHECK!) (ie plumes, times, rhymes, etc)
   ##  if (toupper(substr(Word,nchar(Word)-2,nchar(Word))) == "MES")   
     ##  Num_Sil <- Num_Sil - 1
    ## if (toupper(substr(Word,nchar(Word)-2,nchar(Word))) == "GRES")   
    ##   Num_Sil <- Num_Sil - 1
    
    ##accounts for words that had silent "e" at the end, but were made plural. Need to include attendees, boxes, imagines, performances etc
    if (toupper(substr(Word,nchar(Word)-1,nchar(Word))) == "ES" && !is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word)-2)),c("C", "X", "E"))) #NOT CES    
        Num_Sil <- Num_Sil - 1
     
  ## subtract "mes" (THERE ARE EXCEPTIONS! NEED TO CHECK!) (ie plumes, times, rhymes, etc)
     if (toupper(substr(Word,nchar(Word)-3,nchar(Word))) == "MENT")   
       Num_Sil <- Num_Sil - 1   

   for (Ind in 4:nchar(Word)) # Queue, Dequeue words
        if (is.element(toupper(substr(Word,Ind-3,Ind)), c("UEUE"))) 
             Num_Sil <- Num_Sil + 1 ## +1 because silent e was subtratcted
   
   for (Ind in 6:nchar(Word)) #words like Worcestershire where the division of cester is worce-ster-shire
        if (is.element(toupper(substr(Word,Ind-5,Ind)), c("CESTER"))) 
             Num_Sil <- Num_Sil - 1 
             
#   if (toupper(substr(fullWord,nchar(fullWord)-2,nchar(fullWord))) == "TED")  #NEED TO FIND ALL CONSONANTS WITH ED THAT CREATED E VOWEL SOUND, uses full word becaue "ed" suffix was chopped off in suffix removal process 
#       Num_Sil <- Num_Sil + 1   


 #HARD CODED WORDS
 if (toupper(Word)=="COLONEL") 
             Num_Sil <- 2
 if (toupper(Word)=="WEDNESDAY") 
             Num_Sil <- 2
 if (toupper(Word)=="TRUANCY") 
             Num_Sil <- 3
 if (toupper(Word)=="NUANCE")  ## NUAN, quantify - if there is an N after UA, then add back, except for QUA???
             Num_Sil <- 2
 if (toupper(Word)=="CIAO") 
             Num_Sil <- 1            
             
return(Num_Sil)

}

eDictionary=read.csv2("silentEWords.txt", header=FALSE)
separateSilentE<- function(Word, eDictionary)
{
	eWord=""
	wordVector=c(Word)
	containsE=FALSE
	for(i in 1: length(eDictionary$V1))
	{
		eWord=as.character(eDictionary$V1[i])
		if(toupper(substr(Word,1,nchar(eWord)))==toupper(eWord) && nchar(eWord)!=nchar(Word))
		{
			containsE=TRUE
			break
		}
	}
	if(containsE)
	{
		wordVector=c(substr(Word,1,nchar(eWord)), substr(Word,nchar(eWord)+1, nchar(Word) ))
	}
	return(wordVector)
}