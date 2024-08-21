#libraries
library(rvest)
library(magrittr)
library(stringr)
library(tokenizers)


initialPat=function(pos)
 {
 	initialPattern=""
 	if(pos=="DET")
 	{
 		initialPattern="definite article"
 	}
 	if(pos=="ADJ")
 	{
 		initialPattern="adjective"
 	}
 	
 	if(pos=="ADV")
 	{
 		initialPattern="adverb"
 	}
 	
 	if(pos=="VERB")
 	{
 		initialPattern="verb"
 	}
 	
 	if(pos=="NOUN")
 	{
 		initialPattern="noun"
 	}
 	
 	if(pos=="PUNCT")
 	{
 		return() #no point of punctuation
 	}
 	
 	if(pos=="PROPN")
 	{
 		initialPattern=c("noun", "geographical name (1)", "biographical name (1)")
 	}
 	
 	if(pos=="AUX")
 	{
 		initialPattern="verb"
 	}
 	
 	if(pos=="SCONJ")
 	{
 		initialPattern="conjunction"
 	}
 	
 	if(pos=="CCONJ")
 	{
 		initialPattern="conjunction"
 	}
 	
 	if(pos=="PRON")
 	{
 		initialPattern="pronoun"
 	}
 	
 	if(pos=="ADP")
 	{
 		initialPattern="preposition"
 	}
 	
 	if(pos=="PART")
 	{
 		initialPattern=c("preposition", "adverb", "adjective")
 	}
 	if(pos=="X")
 	{
 		initialPattern="noun"
 	}

 	
 	return(initialPattern)

 }



#function to clean cleanPOS string 
cleanPOS= function(Array, Word){
	array=Array
	Vowels=c("a", "e", "i", "o", "u", "y")
	for(i in 1:length(array)){
		entry=strsplit(array[i], split="")
		if(!grepl(paste(Vowels, collapse="|"), array[i]) || grepl(tolower(Word), array[i]))
		{
			array[i]=""
		}
		comma=unlist(gregexpr(",", array[i]))[1]
		if(comma!=-1)
		{
			array[i]= substr(array[i], 1,comma-1)
		}
	}
	array=trimws(array[which(array!="")])
	return(array)
}

cleanPRON=function(string)
{
	doubleSlashInd= unlist(gregexpr("\\", string, fixed=TRUE))
	n=length(doubleSlashInd)/2
	Pron= vector("character", n)
	counter=1
	for(i in 1: n)
	{
		substr=substr(string,doubleSlashInd[counter], doubleSlashInd[counter+1])
		
		#if multiple pronounciations for given POS , take first
		if(grepl(",", substr)){
			index=unlist(gregexpr(",",substr))
			substr=substr(substr,1, index-1)
		}
		Pron[i]= substr
		counter= counter+2
	}
	Pron= trimws(gsub("\\", "", Pron, fixed=TRUE))
	
	return(Pron)
}






#sAs =syllables And stress

readUrl <- function(url) {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully

            #message("This is the 'try' part")

            read_html(url)
            # The return value of `read_html(url)` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            #message(paste("URL does not seem to exist:", url))
            #message("Here's the original error message:")
            #message(cond)
            # Choose a return value in case of error
            print("Dictionary could not locate webpage (error)")
            return(1)
        },
        warning=function(cond) {
            message(paste("URL caused a warning:", url))
            message("Here's the original warning message:")
            message(cond)
            # Choose a return value in case of warning
            return(NULL)
        },
        finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
            #message(paste("Processed URL:", url))
            #message("Some other message at the end")
        }
    )    
    return(out)
}



sAsENG= function(Word, type)
{
	url=paste("https://www.merriam-webster.com/dictionary/",tolower(Word), sep="")
	html= readUrl(url)  #try catch 
	if(typeof(html)=="double") #check if word exists
	{
		return(NA)
	}
	entry= html %>% html_nodes(".hword ") %>% html_text() %>% strsplit(split="\n")%>% unlist()
	entry=entry[1]
	
	if(length(unlist(tokenize_words(entry)))>1) # angeles returning East Los Angeles
	{
		return(NA)
	}
	    	
	#find html portion which contains all word part of speeches for the dictionary
	cleanPos= html %>% html_nodes(".col-lg-12  ") %>% html_text() %>% strsplit(split="\n")%>% 			unlist()
	
	#clean the string
	POS=0
	if(!is.null(cleanPos))
	{
		POS= cleanPOS(cleanPos, entry) 
	}

	#find html portion with all associated pronounciations 
	cleanPron= html %>% html_nodes(".col ") %>% html_text() %>% strsplit(split="\n")%>% unlist()%>% 		paste(collapse="")
	
	PRON= cleanPRON(cleanPron)
	if(is.na(PRON[1]) &&sum(POS=="abbreviation")==0){ #no pronounciation given and NOT abbreviation   #will get first entry from cambridge dictionary, planB
		print("Calling Cambridge dictionary because phonetic spelling was not given")
		return(sAsENGBackUp(Word, type))
		}
	
	loc=1
	wordType= initialPat(type)
	if(!POS==0)
	{
	found=FALSE
	i=1
	if(length(wordType)>1)
	{
		while(!found && i<=length(wordType))
		{
			loc= which(POS==wordType[i])
			if(length(loc)>0)
				found=TRUE
			else
				i=i+1	
		}
	}
	if(length(wordType)>1 && found==FALSE)
	{
		return(NA)
	}
	
	loc= which(POS==wordType[i])
	}

	
	
	if(length(loc)>1)
		loc=loc[1]
	if(length(loc)==0 && sum(POS=="abbreviation")>=1) #put here to check for abbreviation in case abbreviation can be both a word and an abbreviation
		return("abb")
	if(length(loc)==0 && length(PRON)!=0)
		loc=1
	if(loc>length(PRON)) #if there isnt a pronounciation for the POS
		loc=1
	#proper phonetic spelling to find info from
	decodeStr= PRON[loc]
	also=unlist(gregexpr("also",decodeStr)) #account for 2 pronounciations in phonetic spelling
	if(also!=-1)
	{
		decodeStr= substr(decodeStr,1, also-1)
	}
	
	## Address "couple" and "family". 
	
	colonLocs = regexpr(";", decodeStr) 
	if(colonLocs[1]!=-1)
	{
		decodeStr = substr(decodeStr, 1, colonLocs[1]-1)
		decodeStr = trimws(decodeStr)
	}
	
	##
	
	
	numSyl= str_count(decodeStr, "-")+1
	p1=unlist(gregexpr("-​)", decodeStr))[1]
	p2=unlist(gregexpr("(-", decodeStr, fixed=TRUE))[1]
	if(p1!=-1 ||p2!=-1 ) # case with (e-) in phonetic spelling
	{
		if(p1==-1)
		p1=NULL
		if(p2==-1)
		p2=NULL
		
		numSyl=numSyl-length(p1)-length(p2)
	}
	stressPos= unlist(gregexpr("ˈ",decodeStr))
 	substrB4Stress=substr(decodeStr, 1,stressPos)
 	numDashB4=unlist(gregexpr("-", substrB4Stress))
 	
 	if(numDashB4[1]==-1){
 		numDashB4=0
 		
 	}else{
 		numDashB4=length(numDashB4)
 	}
 		
 	Stress= numDashB4+1

#####################################
#check if word and entry are the same
#####################################
	
	if(Word!=entry)
	{
     if (is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("TED","DED", "FUL", "ING"))) #escalated, suspected, decided, suspenseful
      numSyl=numSyl+1
      
        if (is.element(toupper(substr(Word,nchar(Word)-3,nchar(Word))),c("CHES"))) #stretches, matches
      numSyl=numSyl+1 
      
      if(is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("IER")) && is.element(toupper(substr(entry, nchar(entry), nchar(entry))), c("Y")))
      numSyl=numSyl+1 
	}
	
	if(tolower(entry)!=tolower(Word))
	print("ENTRY AND WORD DO NOT MATCH!!!!!")
	
    print(paste("Word:", Word, "|  Entry:", entry))
    print(paste("numSyl:", numSyl,  "|  Stress:", Stress))
    
    closeAllConnections()
	return(c(numSyl, Stress))
	
}

cleanBrit= function(Array, Word){
	outputArray= vector("character", 4) #will have entry, POS, uk/us, pronunciation
	array=Array
	Vowels=c("a", "e", "i", "o", "u", "y")
	for(i in 1:length(array)){
		
		if(grepl(tolower(Word), tolower(array[i])) && !length(unlist(gregexpr("/", array[i])))>1){ #try to grab POS and uk/us 
			
			outputArray[1]= Word #because detected in if condition
			array[i]= gsub(Word, "", array[i]) #word already matched in if stmt
			strt= unlist(regexpr("[", array[i], fixed=TRUE))[1]-1
			if(!strt<0)
			{
				end= unlist(regexpr("]", array[i], fixed=TRUE))[1]+1
				array[i]= paste(substr(array[i], 1, strt), substr(array[i], end, 						nchar(array[i])))
			}
			if(grepl(" us ", array[i])[1])
			{
				outputArray[3]= "us"
				strtRM= unlist(regexpr(" us ", array[i]))[1]
				array[i]= paste(substr(array[i],1, strtRM), substr(array[i],      						strtRM+4,nchar(array[i])))
				outputArray[2]= trimws(array[i])
			}
			if(grepl(" uk ", array[i]))
			{
				outputArray[3]= "uk"
				strtRM= unlist(regexpr(" uk ", array[i]))[1]
				array[i]= paste(substr(array[i],1, strtRM), substr(array[i],      						strtRM+4,nchar(array[i])))
				outputArray[2]= trimws(array[i])
			}
			
			
		}else if(grepl("/", array[i])){ # grab pronunciation
			strt=regexpr("/", array[i], fixed=TRUE)
			end=regexpr(", ", array[i], fixed=TRUE)-1
			if(end[1]<=0){
				end=unlist(gregexpr("/", array[i], fixed=TRUE))[2]
			}
			
			outputArray[4]= substr(array[i], strt[1], end)
			break
		}
	}
	return(outputArray)
}


sAsBRIT= function(Word, type)
{
	url=paste("https://dictionary.cambridge.org/us/dictionary/english/",tolower(Word), sep="")
	html= readUrl(url) #try catch 
	
	if(typeof(html)=="double")
	{
		return(NA)
	}

	#cleanScript has everything we need, POS and pronunciation
	cleanScript= html %>% html_nodes(".dpos-h") %>% html_text() %>% strsplit(split="\n") 
	
	entries= html %>% html_nodes(".dhw") %>% html_text()
	entry=entries[1]
	#find proper entry 
	if(!is.na(entry)&& Word!=entry)
	{
		if(!is.na(which(entries==Word)[1])){
			entry=entries[which(entries==Word)[1]]
		}
	}
	
	if(length(cleanScript)==0 ||length(unlist(tokenize_words(entries[1])))>1)
	{
		print("Calling Merriam Webster dictionary because phonetic spelling was not given or  dictionary entry may have been incorrect")
		return(sAsBRITBackUp(Word, type))
	}
	
	found=FALSE
	
	for(i in 1: length(cleanScript))
	{
		cleanScript[[i]]=cleanBrit(cleanScript[[i]], entry)
		if(cleanScript[[i]][1]!="" && cleanScript[[i]][4]!="" ){
			found=TRUE
		}
	}
	k=2	
	while(!found && k<=length(entries))
	{
		cleanScript= html %>% html_nodes(".dpos-h") %>% html_text() %>% 						strsplit(split="\n") 
		entry=entries[k]
		#find proper entry 
	
		for(i in 1: length(cleanScript))
		{
			cleanScript[[i]]=cleanBrit(cleanScript[[i]], entry)
			if(cleanScript[[i]][1]!="" && cleanScript[[i]][4]!="")
			{
				found=TRUE
			}
		}

	 k=k+1	
	}	
	pos=type
	type= initialPat(type)
	ukTruth= vector("numeric", length(cleanScript))
	posTruth= vector("numeric", length(cleanScript))	
	pronTruth= vector("numeric", length(cleanScript))	
	for(i in 1: length(cleanScript))
	{
		if(!is.na(cleanScript[[i]][2]))
		{
			if(is.element(cleanScript[[i]][2], type))
			{
				posTruth[i]=1
			}
			if(cleanScript[[i]][3]=="uk")
			{
				ukTruth[i]= 1
			}
			if(!is.na(cleanScript[[i]][4]) && cleanScript[[i]][4]!="")
			{
				pronTruth[i]= 1
			}
		}
	}
	totTruth=ukTruth + posTruth+pronTruth
	
	if(sum(totTruth)==0)
	{
		print("Calling Merriam Webster dictionary because phonetic spelling was not given or  dictionary entry may have been incorrect")
		return(sAsBRITBackUp(Word, pos))
	}
	#look for entry with proper part of speech and uk 
	index= which(totTruth==2)[1]
	#if cant find entry with both uk
	if(is.na(index[1]))
	{
		#use entry with correct POS if exists
		if(sum(posTruth)!=0){
		index= which(posTruth==1)[1]
		} else{
		#if POS isnt there, use first entry of dictionary
		index=1
		}
	}
	#now proper entry is stored in index
	decodeStr=""
	if(!is.na(cleanScript[[index]][4]))
		index= which(pronTruth==1)[1]
	{
		decodeStr= gsub("/", "", cleanScript[[index]][4], fixed=TRUE)
	}
	if(length(decodeStr)==0)
	{
		return(NA) #will revert to old for number of syllables
	}

	sylIndicators=c(".", "·", "ˌ", "ˈ")
	stress=1
 	numSyls=1
 	for(i in 1:nchar(decodeStr))
 	{
 		if(is.element(substr(decodeStr, i, i), sylIndicators))
 		{
 			numSyls=numSyls+1
 			if(substr(decodeStr, i, i)=="ˈ"){
 				stress=numSyls
 			}
 					 					
 		}
 	}
	if(substr(decodeStr, 1, 1)=="ˈ")
 			{
 				numSyls=numSyls-1
 				stress=stress-1
 			}
 	if(substr(decodeStr, 1, 1)=="ˌ") #phonetic spelling for installation
 			{
 				numSyls=numSyls-1
 				stress=stress-1
 			}		
 			
 	if(Word!=entry)
	{
     if (is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("TED","DED", "FUL", "ING"))) #escalated, suspected, decided, suspenseful
      numSyls=numSyls+1
      
        if (is.element(toupper(substr(Word,nchar(Word)-3,nchar(Word))),c("CHES"))) #stretches, matches
      numSyls=numSyls+1 
      
      if(is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("IER")) && is.element(toupper(substr(entry, nchar(entry), nchar(entry))), c("Y")))
      numSyls=numSyls+1 
	}
	
	if(tolower(entry)!=tolower(Word))
	print("ENTRY AND WORD DO NOT MATCH!!!!! ")
	
    print(paste("Word:", Word, "|  Entry:", entry))
    print(paste("numSyl:", numSyls,  "|  Stress:", stress))				
	closeAllConnections()
	return(c(numSyls, stress))
}

sAsUser= function(Word){
	print(paste("The word:'", Word, "' could not be found in the dictionary.  User input is necessary.", sep=""))
	numSyls= readline(prompt=paste("Enter the number of syllables for: '", Word, "'\n", sep=""))
	stress=  readline(prompt=paste("Enter the location of the word stress (which syllable it occurs on) for: '", Word, "'\n", sep=""))
	return(as.integer(c(numSyls, stress)))	
} 


sAsENGBackUp= function(Word, type)
{
url=paste("https://dictionary.cambridge.org/us/dictionary/english/",tolower(Word), sep="")
	html= readUrl(url) #try catch 
	
	if(typeof(html)=="double")
	{
		return(NA)
	}

	#cleanScript has everything we need, POS and pronunciation
	cleanScript= html %>% html_nodes(".dpos-h") %>% html_text() %>% strsplit(split="\n") 
	
	entries= html %>% html_nodes(".dhw") %>% html_text()
	entry=entries[1]
	#find proper entry 
	if(!is.na(entry)&& Word!=entry)
	{
		if(!is.na(which(entries==Word)[1])){
			entry=entries[which(entries==Word)[1]]
		}
	}
	
	if(length(cleanScript)==0 ||length(unlist(tokenize_words(entries[1])))>1)
	{
		return(NA)
	}
	
	found=FALSE
	
	for(i in 1: length(cleanScript))
	{
		cleanScript[[i]]=cleanBrit(cleanScript[[i]], entry)
		if(cleanScript[[i]][1]!="" && cleanScript[[i]][4]!="" ){
			found=TRUE
		}
	}
	k=2	
	while(!found && k<=length(entries))
	{
		cleanScript= html %>% html_nodes(".dpos-h") %>% html_text() %>% 						strsplit(split="\n") 
		entry=entries[k]
		#find proper entry 
	
		for(i in 1: length(cleanScript))
		{
			cleanScript[[i]]=cleanBrit(cleanScript[[i]], entry)
			if(cleanScript[[i]][1]!="" && cleanScript[[i]][4]!="")
			{
				found=TRUE
			}
		}

	 k=k+1	
	}	
	pos=type
	type= initialPat(type)
	posTruth= vector("numeric", length(cleanScript))	
	pronTruth= vector("numeric", length(cleanScript))	
	for(i in 1: length(cleanScript))
	{
		if(!is.na(cleanScript[[i]][2]))
		{
			if(is.element(cleanScript[[i]][2], type))
			{
				posTruth[i]=1
			}
			if(!is.na(cleanScript[[i]][4]) && cleanScript[[i]][4]!="")
			{
				pronTruth[i]= 1
			}
		}
	}
	totTruth=posTruth+pronTruth
	
	if(sum(totTruth)==0)
	{
		
		return(NA)
	}

	index= which(totTruth==2)[1]
	
	if(is.na(index[1]))
	{
		#use entry with correct POS if exists
		if(sum(posTruth)!=0){
		index= which(posTruth==1)[1]
		} else{
		#if POS isnt there, use first entry of dictionary
		index=1
		}
	}
	#now proper entry is stored in index
	decodeStr=""
	if(!is.na(cleanScript[[index]][4]))
		index= which(pronTruth==1)[1]
	{
		decodeStr= gsub("/", "", cleanScript[[index]][4], fixed=TRUE)
	}

	sylIndicators=c(".", "·", "ˌ", "ˈ")
	stress=1
 	numSyls=1
 	for(i in 1:nchar(decodeStr))
 	{
 		if(is.element(substr(decodeStr, i, i), sylIndicators))
 		{
 			numSyls=numSyls+1
 			if(substr(decodeStr, i, i)=="ˈ"){
 				stress=numSyls
 			}
 					 					
 		}
 	}
	if(substr(decodeStr, 1, 1)=="ˈ")
 			{
 				numSyls=numSyls-1
 				stress=stress-1
 			}
 	if(substr(decodeStr, 1, 1)=="ˌ") #phonetic spelling for installation
 			{
 				numSyls=numSyls-1
 				stress=stress-1
 			}		
 			
 	if(Word!=entry)
	{
     if (is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("TED","DED", "FUL", "ING"))) #escalated, suspected, decided, suspenseful
      numSyls=numSyls+1
      
        if (is.element(toupper(substr(Word,nchar(Word)-3,nchar(Word))),c("CHES"))) #stretches, matches
      numSyls=numSyls+1 
      
      if(is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("IER")) && is.element(toupper(substr(entry, nchar(entry), nchar(entry))), c("Y")))
      numSyls=numSyls+1 
	}
	
	if(tolower(entry)!=tolower(Word))
	print("ENTRY AND WORD DO NOT MATCH!!!!! ")
	
    print(paste("Word:", Word, "|  Entry:", entry))
    print(paste("numSyl:", numSyls,  "|  Stress:", stress))
				
	closeAllConnections()
	return(c(numSyls, stress))

}

sAsBRITBackUp= function(Word, type)
{
	url=paste("https://www.merriam-webster.com/dictionary/",tolower(Word), sep="")
	html= readUrl(url)  #try catch 
	if(typeof(html)=="double") #check if word exists
	{
		return(NA)
	}
	entry= html %>% html_nodes(".hword ") %>% html_text() %>% strsplit(split="\n")%>% unlist()
	entry=entry[1]
	
	if(length(unlist(tokenize_words(entry)))>1) # angeles returning East Los Angeles
	{
		return(NA)
	}
	    	
	#find html portion which contains all word part of speeches for the dictionary
	cleanPos= html %>% html_nodes(".col-lg-12  ") %>% html_text() %>% strsplit(split="\n")%>% 			unlist()
	
	#clean the string
	POS=0
	if(!is.null(cleanPos))
	{
		POS= cleanPOS(cleanPos, entry) 
	}
	#find html portion with all associated pronounciations 
	cleanPron= html %>% html_nodes(".col ") %>% html_text() %>% strsplit(split="\n")%>% unlist()%>% 		paste(collapse="")
	
	PRON= cleanPRON(cleanPron)
	if(is.na(PRON[1]) &&sum(POS=="abbreviation")==0){ #no pronounciation given and NOT abbreviation   #will get first entry from cambridge dictionary, planB
		print("Entry not found in backup Webster dictionary")
		return(NA)
		}

	
	loc=1
	wordType= initialPat(type)
	if(!POS==0)
	{
	found=FALSE
	i=1
	if(length(wordType)>1)
	{
		while(!found && i<=length(wordType))
		{
			loc= which(POS==wordType[i])
			if(length(loc)>0)
				found=TRUE
			else
				i=i+1	
		}
	}
	if(length(wordType)>1 && found==FALSE)
	{
		return(NA)
	}
	
	loc= which(POS==wordType[i])
	}
	
	
	if(length(loc)>1)
		loc=loc[1]
	if(length(loc)==0 && sum(POS=="abbreviation")>=1) #put here to check for abbreviation in case abbreviation can be both a word and an abbreviation
		return("abb")
	if(length(loc)==0 && length(PRON)!=0)
		loc=1
	if(loc>length(PRON)) #if there isnt a pronounciation for the POS
		loc=1
	#proper phonetic spelling to find info from
	decodeStr= PRON[loc]
	also=unlist(gregexpr("also",decodeStr)) #account for 2 pronounciations in phonetic spelling
	if(also!=-1)
	{
		decodeStr= substr(decodeStr,1, also-1)
	}
	colonLocs = regexpr(";", decodeStr) 
	if(colonLocs[1]!=-1)
	{
		decodeStr = substr(decodeStr, 1, colonLocs[1]-1)
		decodeStr = trimws(decodeStr)
	}
	
	
	numSyl= str_count(decodeStr, "-")+1
	p1=unlist(gregexpr("-​)", decodeStr))[1]
	p2=unlist(gregexpr("(-", decodeStr, fixed=TRUE))[1]
	if(p1!=-1 ||p2!=-1 ) # case with (e-) in phonetic spelling
	{
		if(p1==-1)
		p1=NULL
		if(p2==-1)
		p2=NULL
		
		numSyl=numSyl-length(p1)-length(p2)
	}
	stressPos= unlist(gregexpr("ˈ",decodeStr))
 	substrB4Stress=substr(decodeStr, 1,stressPos)
 	numDashB4=unlist(gregexpr("-", substrB4Stress))
 	
 	if(numDashB4[1]==-1){
 		numDashB4=0
 		
 	}else{
 		numDashB4=length(numDashB4)
 	}
 		
 	Stress= numDashB4+1

#####################################
#check if word and entry are the same
#####################################
	
		if(Word!=entry)
	{
     if (is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("TED","DED", "FUL", "ING"))) #escalated, suspected, decided, suspenseful
      numSyl=numSyl+1
      
        if (is.element(toupper(substr(Word,nchar(Word)-3,nchar(Word))),c("CHES"))) #stretches, matches
      numSyl=numSyl+1 
      
      if(is.element(toupper(substr(Word,nchar(Word)-2,nchar(Word))),c("IER")) && is.element(toupper(substr(entry, nchar(entry), nchar(entry))), c("Y")))
      numSyl=numSyl+1 
	}
	
	if(tolower(entry)!=tolower(Word))
	print("ENTRY AND WORD DO NOT MATCH!!!!!")
	
	print("From Merriam Webster Dictionary...")
    print(paste("Word:", Word, "|  Entry:", entry))
	print(paste("numSyl:", numSyl,  "|  Stress:", Stress))
	

	closeAllConnections()
	return(c(numSyl, Stress))
	
}



