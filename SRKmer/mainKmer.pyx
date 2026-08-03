cpdef myFunction(object numIn, yqList, centList, yqSubList, centSubList):
	
	
	cdef int total_reads = 0, x, k, kmer=24, readLength
	cdef str readName, read
	cdef list baseCall

	#For each read in the cram file
	for read1 in numIn.fetch(until_eof=True):

		#Add 1 to the total reads
		total_reads+=1
		#Grab the read name
		#readName = str(read1).split('\t')[0]

		#Skip records with no stored sequence or no base qualities.
		#Line 22 below parses the quality array out of str(read1); when a
		#record has QUAL '*' that field is the literal text 'None' and the
		#parse raises IndexError. Secondary and supplementary alignments
		#are the usual source of such records.
		if read1.query_sequence is None or read1.query_qualities is None:
			continue

		#Grab the read sequence
		read = str(read1).split('\t')[9]
		readLength = len(read)

		for x in range(0, readLength - kmer + 1, 1):
			if read[x:x+kmer] in yqList:
				baseCall = [int(qscore.strip()) for qscore in str(read1).split("\t")[10].split('[')[1].split(']')[0].split(",")]
				if len([k for k in baseCall[x:x+kmer] if k >=20]) == 24:
					yqList['lengths']+= int(readLength)
					yqList[read[x:x+kmer]].append(str(str(read1).split('\t')[0])+"_"+str(readLength))
					yqList['total_reads']+=1
					if read[x:x+kmer] in yqSubList:
						yqSubList[read[x:x+kmer]].append(str(str(read1).split('\t')[0])+"_"+str(readLength))
						yqSubList['lengths']+= int(readLength)
						yqSubList['total_reads']+=1
					else:
						pass
					break
				else:
					continue
			
			elif read[x:x+kmer] in centList:
				baseCall = [int(qscore.strip()) for qscore in str(read1).split("\t")[10].split('[')[1].split(']')[0].split(",")]
				if len([k for k in baseCall[x:x+kmer] if k >=20]) == 24:
					centList['total_reads']+=1
					centList[read[x:x+kmer]].append(str(str(read1).split('\t')[0])+"_"+str(readLength))
					centList['lengths']+= int(readLength)
					if read[x:x+kmer] in centSubList:
						centSubList[read[x:x+kmer]].append(str(str(read1).split('\t')[0])+"_"+str(readLength))
						centSubList['lengths']+= int(readLength)
						centSubList['total_reads']+=1
					else:
						pass
					break
				else:
					continue
			else:
				continue

	return (total_reads, yqList ,centList, yqSubList, centSubList)
