# Funzioni relative solo al progetto A

def inizializzazione(fileInput):
	"""	Inizializza le strutture dati per il progetto A.

	Tali struttura sono necessarie all'elaborazione degli introni, al loro
	allineamento e alla produzione del file di output.

	Argomenti:
		nomeFile		:	nome del file di annotazione in input

	Ritorna:
		dictTranscript	:	dizionario informazioni sui Trascritti
		dictGenes		:	dizionario informazioni sui Geni
		dictEsoni		:	dizionario Esoni per TrascrittoID
		dictIntroni		:	dizionario Introni per TrascrittoID     
		dictGeneChr		:	dizionario coordinate Geni per Cromosomi
	"""
	
	dictTranscript 	= {}
	dictGenes 		= {}
	dictEsoni 		= {}
	dictIntroni 	= {}
	dictGeneChr 	= {}		

	# - Filtraggio file di annotazione in input
	# - Calcolo delle coordinate dei geni nei cromosomi
	lines, dictGeneChr = filtraFileDiAnn(fileInput)
	
	# Indici all'interno del dizionario degli esoni
	idx_starts 	= 0
	idx_ends	= 1
	idx_strand 	= 2
	
	# Indici all'interno del dizionario dei Geni
	idx_transcripts = 2


	# Creazione dei dizionari utili alla risoluzione del problema A
	for riga in lines:
		cromosoma		= riga[0]
		start_esone		= riga[3]
		end_esone 		= riga[4]
		strand 			= riga[6]
		geneName 		= riga[11]
		transcriptName 	= riga[12]
		
		TranscriptID 	= riga[9]
		GeneID 			= riga[8]
	
		# Creazione del dizionario dei transcritti
		dictTranscript[TranscriptID] = [transcriptName, GeneID]
		
		# Creazione del dizionario dei geni
		if not dictGenes.has_key(GeneID):								# Se il GeneID non e' presente..
			dictGenes[GeneID] = [geneName, cromosoma, [TranscriptID]]	# ..nel dizionario (come key)
		
		elif TranscriptID not in dictGenes[GeneID][idx_transcripts]:	# Se il GeneID e' presente ma non lo e'..
			dictGenes[GeneID][idx_transcripts].append(TranscriptID)     # ..il TranscriptID si aggiunge questo alla lista
		
		# Creazione del dizionario degli esoni
		if not dictEsoni.has_key(TranscriptID):						    	# Se il TranscriptID non e' presente.. 
			dictEsoni[TranscriptID] = [[start_esone],[end_esone],strand]    # ..nel dizionario (come key)
		
		else:
			dictEsoni[TranscriptID][idx_starts].append(start_esone)		# Il TranscriptID e' gia' presente quindi..
			dictEsoni[TranscriptID][idx_ends].append(end_esone)			# ..si aggiunge l'esone alla lista degli esoni
			

	# Creazione del dizionario degli introni
	for TranscriptID in dictEsoni:
		esoniPerTranscript = len(dictEsoni[TranscriptID][idx_starts])	# Si valuta il nr di esoni per TranscriptID corrente
		
		if int(esoniPerTranscript) > 1:
			start_introni 	= []										# Si preparano le variabili necessarie
			end_introni 	= []
			start_esoni 	= []													
			end_esoni 		= []
			
			# Si considera lo strand relativo al TranscriptID
			if dictEsoni[TranscriptID][idx_strand] == '+':				# Strand positivo -> esoni scritti in ordine crescente
				strand = True
				start_esoni = dictEsoni[TranscriptID][idx_starts]
				end_esoni 	= dictEsoni[TranscriptID][idx_ends]
				
			else:														# Strand negativo -> esoni scritti in ordine inverso..
				strand = False											# ..e per comodita' sono invertiti in ordine crescente
				start_esoni = dictEsoni[TranscriptID][idx_starts][::-1] 	  
				end_esoni 	= dictEsoni[TranscriptID][idx_ends][::-1]

			# Calcolo delle regioni introniche
			i = 0
			while i < int(esoniPerTranscript) - 1:						# Per ogni coppia di esoni
				if (int(start_esoni[i+1]) - int(end_esoni[i])) > 2:		# Se la regione tra due esoni consecutivi e' > 2..
					start_introni.append(int(end_esoni[i]) + 1)			# ..(considerando che gli estremi dell'introne sono..
					end_introni.append(int(start_esoni[i+1]) - 1)		# ..interni a quelli dei due esoni consecutivi correnti)
				
				i += 1
			
			if not strand:												# Si mantiene traccia del fatto che derivano da un..
				start_introni.reverse()									# ..TranscriptID con strand negativo..
				end_introni.reverse()									# ..(si inverte l'ordine degli introni)
			
			dictIntroni[TranscriptID] = [start_introni, end_introni]


	return [dictTranscript, dictGenes, dictEsoni, dictIntroni, dictGeneChr]

# ---------------------------------------------------------------------

def filtraFileDiAnn(nomeFile):
	"""	Filtra il file di annotazione inserito dall'utente.

	Estrae solamente gli esoni e per ogni gene, trova i suoi estremi
	all'interno di ciascun cromosoma al quale appartiene.

	Argomenti:
		nomeFile		:	nome del file di annotazione in input

	Ritorna:
		lines			:	lista delle righe del file di annotazione
							relative solamente agli esoni
		dictGeneChr		:	dizionario con inizio e fine per ogni
							corrispondenza trovata nel file di
							annotazione
	"""

	#---------------------
	# Creazione di una lista dove ogni elemento e' una riga del file 
	# Ogni elem e' una lista di informazioni divise per colonne 
	#
	# formato di un elemento di lines:
	#
	#	POSIZIONE 			CONTENUTO
	#		0					cromosoma
	#		3					start
	#		4					end
	#		6					strand
	#		8					gene_id
	#		9					transcript_id
	#		10					exon_number
	#		11					gene_name
	#		12					transcript_name	
	#


	stringa 	= '\texon\t'
	lines 		= []
	dictGeneChr = {}
	
	# Indici per il file di annotazione
	idx_cromosoma = 0
	idx_geneName  = 11
	idx_start     = 3
	idx_end       = 4
	
	for x in open(nomeFile):
		riga = x.strip(';\n').replace('; ','\t').split('\t')
				
		# Creazione del dizionario dei gene_name per ogni cromosoma
		key_geneChr = riga[idx_geneName] + '\t' + riga[idx_cromosoma]
		if not dictGeneChr.has_key(key_geneChr):
			dictGeneChr[key_geneChr] = [riga[idx_start], riga[idx_end]]	# Calcolo dello start e end del gene
																		# (serve per output collassato su introni)
		else:
			# Si aggiona il valore dello start del gene se si trova un 
			# valore piu' piccolo
			if int(dictGeneChr[key_geneChr][0]) > int(riga[idx_start]):
				dictGeneChr[key_geneChr][0] = riga[idx_start]
				
			# Si aggiorna il valore dell'end del gene se si trova un
			# valore piu' grande
			if int(dictGeneChr[key_geneChr][1]) < int(riga[idx_end]):	
				dictGeneChr[key_geneChr][1] = riga[idx_end]
		
		# Si filtra il file considerando solamente le regioni 
		# di tipo "exon"
		if stringa in x:
			lines.append(riga)


	return [lines, dictGeneChr]

