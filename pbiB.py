from os import path
from os import system

cartella = '%s/'														# Cartella relativa ad ogni gene inserito dall'utente

# Header per il file in formato BedGraph
#
intestazione = """browser position %s
browser hide all
browser pack refGene encodeRegions
browser full altGraph
track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20
"""

# ---------------------------------------------------------------------

def inizializzazione(fileInput, geneNames):
	"""	Inizializza le strutture dati per il progetto B.

	Tali struttura sono necessarie all'elaborazione degli introni, 
	al loro allineamento e alla produzione del file di output.

	Argomenti:
		fileInput		:	nome del file di annotazione in input
		geneNames		:	dizionario dei geni inseriti dall'utente

	Ritorna:
		[
		dictTranscript,	:	dizionario informazioni sui Trascritti
		dictGenes,		:	dizionario informazioni sui Geni
		dictEsoni,		:	dizionario Esoni per TrascrittoID
		dictIntroni,    : 	dizionario Introni per TrascrittoID     
		dictGeneChr		:	dizionario coordinate Geni per Cromosomi
		]
	"""
	dictTranscript 	= {}
	dictGenes 		= {}
	dictEsoni 		= {}
	dictIntroni 	= {}
	dictGeneChr 	= {}

	# - Filtraggio file di annotazione in input per 'exon' e 
	#   per nome gene
	# - Calcolo delle coordinate dei geni nei cromosomi
	lines, dictGeneChr = filtraFileDiAnn(fileInput, geneNames)
	
	
	# Indici all'interno del dizionario degli esoni
	idx_starts 	= 0
	idx_ends 	= 1
	idx_strand 	= 2
	
	# Indici all'interno del dizionario dei Geni
	idx_transcripts = 2


	# Creazione dei dizionari utili alla risoluzione del problema B
	for riga in lines:
		cromosoma 		= riga[0]
		start_esone 	= riga[3]
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
			dictGenes[GeneID][idx_transcripts].append(TranscriptID)		# ..il TranscriptID questo si aggiunge alla lista
		
		# Creazione del dizionario degli esoni
		if not dictEsoni.has_key(TranscriptID):						     # Se il TranscriptID non e' presente.. 
			dictEsoni[TranscriptID] = [[start_esone],[end_esone],strand] # ..nel dizionario (come key)
		
		else:
			dictEsoni[TranscriptID][idx_starts].append(start_esone)			 # Il TranscriptID e' gia' presente quindi..
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
					end_introni.append(int(start_esoni[i+1]) - 1)		#..interni a quelli dei due esoni consecutivi correnti)
				i += 1
			
			if not strand:												# Si mantiene traccia del fatto che derivano da un..
				start_introni.reverse()									# ..TranscriptID con strand negativo..
				end_introni.reverse()									# ..(si inverte l'ordine degli introni)
		
			dictIntroni[TranscriptID] = [start_introni, end_introni]


	# Si eliminano i geni che non presentano regioni introniche:
	# 	- dalla lista di tutti i geni si rimuovono quelli che 
	#     hanno introni;
	#	- dal dizionario si rimuovono quelli rimasti nella lista.
	tuttiIGeni = geneNames.keys()
	for TranscriptID in dictIntroni:
		geneID   = dictTranscript[TranscriptID][1]
		nomeGene = dictGenes[geneID][0]
		
		if nomeGene in tuttiIGeni:
			tuttiIGeni.remove(nomeGene)


	for nomeGene in tuttiIGeni:
		del geneNames[nomeGene]
		print 'Il gene %s non presenta regioni introniche.' % nomeGene


	return [dictTranscript, dictGenes, dictEsoni, dictIntroni, dictGeneChr]

#  --------------------------------------------------------------------

def filtraFileDiAnn(fileInput, geneNames):
	"""	Filtra il file di annotazione inserito dall'utente.

	Estrae solamente gli esoni e per ogni gene, trova i suoi estremi
	all'interno di ciascun cromosoma al quale appartiene.

	Argomenti:
		nomeFile		: 	nome del file di annotazione in input

	Ritorna:
		[
		lines,			: 	lista delle righe del file di annotazione
							relative solamente agli esoni
		dictGeneChr		: 	dizionario con inizio e fine  del gene
							per ogni corrispondenza trovata nel file
							di annotazione
		]
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
	
	for x in open(fileInput):
		riga = x.strip(';\n').replace('; ','\t').split('\t')

		if not geneNames.has_key(riga[idx_geneName]):
			continue
				
		# Creazione del dizionario dei gene_name per ogni cromosoma
		key_geneChr = riga[idx_geneName] + '\t' + riga[idx_cromosoma]
		if not dictGeneChr.has_key(key_geneChr):
			dictGeneChr[key_geneChr] = [riga[idx_start], riga[idx_end]]
		
		else:	
			# Si aggiona il valore dello start del gene se si trova un 
			# valore piu' piccolo
			if int(dictGeneChr[key_geneChr][0]) > int(riga[idx_start]):
				dictGeneChr[key_geneChr][0] = riga[idx_start]
				
			# Si aggiorna il valore dell'end del gene se si trova un
			# valore piu' grande
			if int(dictGeneChr[key_geneChr][1]) < int(riga[idx_end]):	
				dictGeneChr[key_geneChr][1] = riga[idx_end]
		
		
		# Si filtra il file considerando solamente le regioni di 
		# tipo "exon"
		if stringa in x:
			lines.append(riga)


	return [lines, dictGeneChr]

# ---------------------------------------------------------------------

def unisciFile(fileInput, geneNames, fileOut):
	"""
	Unisce i file chiamati 'fileInput' dalle cartelle dei geni.

	La funzione unisce file con nomi uguali in cartelle di geni diversi.
	
	Argomenti:
		fileInput		:	nome di ogni file da unire
		geneNames		:	dizionario con i geni scelti dall'utente
		fileOut			:	file di output
	"""

	nuovoFile = open(fileOut, 'w')

	for gene in geneNames.values():
		f = open(cartella % gene + fileInput).read()					# Apertura del file in write nella cartella del gene..
		nuovoFile.write(f)												# ..corrente


	nuovoFile.close()

#  --------------------------------------------------------------------

def sortRegioni(tupla):
	""" Restituisce la chiave di ordinamento per unisciRegioni.
	
	Viene utilizzata per ordinare 3 liste rispetto alla prima
	creando una tupla (l'ordinamente avviene su base numerica).
	"""
	return int(tupla[0])

#  --------------------------------------------------------------------

def unisciRegioni(dictEsoni, dictIntroni):
	"""	Unisce le strutture dati di regioni esoniche e introniche.

	Argomenti:
		dictEsoni		:	dizionario Esoni per TrascrittoID
		dictIntroni   	:	dizionario Introni per TrascrittoID

	Ritorna:
		dictEsIn		:	dizionario Esoni+Introni per TrascrittoID
	"""
	idx_starts	=	0
	idx_ends	=	1
	idx_tipo	=	2
	
	dictEsIn = {}

	# Inizializzazione del tipo di regioni considerate
	# Servono per stampare il file in formato '.gtf'
	# (relativo ad introni ed esoni)
	esone = True
	introne = False

	for transcriptID in dictEsoni:
		starts_esoni 	= dictEsoni[transcriptID][idx_starts]			# Si salvano le liste relative a..
		ends_esoni		= dictEsoni[transcriptID][idx_ends]				# ..start e end degli esoni

		nrEsoni = len(starts_esoni)
		tipo_esoni		= [esone] * nrEsoni								# Inizializzazione lista di tipo 'esone' 

		if transcriptID in dictIntroni:									# Se il transcript_id ha regioni introniche..
			starts_introni	= dictIntroni[transcriptID][idx_starts]		# ..si salvano le liste relative a start e end..
			ends_introni	= dictIntroni[transcriptID][idx_ends]		# ..degli introni
			
			nrIntroni = len(starts_introni)
			tipo_introni	= [introne] * nrIntroni						# Inizializzazione lista di tipo 'introne'
			
		else:															# Se il transcript non ha regioni introniche..
			starts_introni 	= []										# ..si cancellano le liste per non avere residui..												
			ends_introni	= []										# ..dai transcript precedenti
			tipo_introni	= []

		starts_esoni.extend(starts_introni)								# Concatenazione liste di start di esoni e introni
		ends_esoni.extend(ends_introni)									# Concatenazione liste di end di esoni e introni
		tipo_esoni.extend(tipo_introni)									# Concatenazione liste di tipo di esoni e introni


		# Ordinamento su base numerica delle tre liste rispetto 
		# alla prima (start) ordinata
		regioni = sorted(zip(starts_esoni, ends_esoni, tipo_esoni), key=sortRegioni)

		# Conversione da tupla a tre liste come argomenti nel
		# dizionario dictEsIn
		dictEsIn[transcriptID] = [[], [], []]
		dictEsIn[transcriptID][idx_starts] 	= [t[idx_starts] for t in regioni]
		dictEsIn[transcriptID][idx_ends] 	= [t[idx_ends] for t in regioni]
		dictEsIn[transcriptID][idx_tipo] 	= [t[idx_tipo] for t in regioni]


	return dictEsIn

#  --------------------------------------------------------------------

def caricaReadsEsIn(fileInput):
	"""	Carica in memoria le reads di Esoni e Introni.

	Argomenti:
		fileInput		:	nome del file '.gtf' totale

	Ritorna:
		dictReadsEsIn	:	dizionario con gene e per ogni gene 
							un dizionario di cromosomi con gli 
							intervalli mappati
	"""

	idx_gene 	= 4 
	idx_chrom 	= 0
	idx_start	= 1
	idx_end		= 2
	idx_reads	= 6

	dictReadsEsIn = {}

	lines = [x.strip('\n').split('\t') for x in open(fileInput)]
	
	for riga in lines:
		geneName 	= riga[idx_gene]
		chrom		= riga[idx_chrom]
		start		= riga[idx_start]
		end			= riga[idx_end]
		reads		= riga[idx_reads]

		if not geneName in dictReadsEsIn:
			dictReadsEsIn[geneName] = {}
			dictReadsEsIn[geneName][chrom] = [False, [start], [end], [reads]]	# Il primo campo indica se il cromosoma ha almeno..
		 																		# ..una regione con reads
		elif chrom not in dictReadsEsIn[geneName]:
			dictReadsEsIn[geneName][chrom] = [False, [start], [end], [reads]]
		
		else:
			dictReadsEsIn[geneName][chrom][idx_start].append(start)
			dictReadsEsIn[geneName][chrom][idx_end].append(end)
			dictReadsEsIn[geneName][chrom][3].append(reads)

		i = len(dictReadsEsIn[geneName][chrom][3])
		if int(dictReadsEsIn[geneName][chrom][3][i-1]) != 0:
			dictReadsEsIn[geneName][chrom][0] = True					# Indica se c'e' almeno una regione esonica/intronica..
																		# ..che mappa delle reads

	# Si eliminano i cromosomi che non hanno mappato reads 
	# ne' su introni ne' su esoni (primo value del dizionario = FALSE)
	geneKeys = dictReadsEsIn.keys()
	for geneName in geneKeys:
		chromKeys = dictReadsEsIn[geneName].keys()
		for chrom in chromKeys:
			if not dictReadsEsIn[geneName][chrom][0]:
				del dictReadsEsIn[geneName][chrom]
				# Si eliminano i geni che non hanno piu' cromosomi
				#
				if not dictReadsEsIn[geneName]:
					del dictReadsEsIn[geneName]
					print 'Il gene %s non presenta cromosomi' % geneName,
					print 'con reads mappanti.' 


	return dictReadsEsIn

#  --------------------------------------------------------------------

def stampaBedGraph(dictReadsEsIn, dictGeneChr, fileOut, geneNames):
	"""	Stampa i file '.bedGraph' per ogni cromosoma per ogni gene. 
	
	I file vengono creati uno per ogni cromosoma che presenta reads 
	mappanti, nella cartella del gene che appartiene a tale cromosoma.

	Argomenti:
		dictReadsEsIn	:	dizionario con gene e per ogni gene 
							un dizionario di cromosomi con gli 
							intervalli mappati
		dictGeneChr		:	dizionario con inizio e fine per ogni
							corrispondenza trovata nel file di
							annotazione
	"""

	keyF 			= '%s\t%s'											# Formato key del dizionario dictGeneChr
	coordinateF 	= '%s:%s-%s'										# Formato coordinate del gene all'interno del cromosoma
	rigaF			= '%s\t%s\t%s\t%s\n'								# Formato della riga del file BedGraph
	bedGraphF		= '%s%s.bedGraph'									# Formato nome del file BedGraph

	idx_start	= 1
	idx_end		= 2
	idx_reads	= 3

	for geneName in dictReadsEsIn:
		if geneName not in geneNames:									# Se il gene non presenta regioni introniche..
			continue													# ..non se ne stampa il bedGraph

		cod = geneNames[geneName]
		for chrom in dictReadsEsIn[geneName]:
			startG, endG = dictGeneChr[keyF % (geneName, chrom)]
			pos_browser = coordinateF %(chrom, startG, endG)

			nuovoFile = open(cartella % cod + bedGraphF % (chrom, fileOut), 'w') # Apertura del file nella cartella relativa al gene
			nuovoFile.write(intestazione % pos_browser)							 # Stampa della header all'inizio del file

			for i in range(0, len(dictReadsEsIn[geneName][chrom][idx_start])):
				nuovoFile.write(rigaF % (chrom, \
										 dictReadsEsIn[geneName][chrom][idx_start][i], \
										 dictReadsEsIn[geneName][chrom][idx_end][i], \
										 dictReadsEsIn[geneName][chrom][idx_reads][i]))

			nuovoFile.close()
			

def stampaGTFEsIn(dictTranscript, dictGenes, dictInput, fileOut, geneNames):
	"""	Stampa file in formato '.gtf'.
	
	Il contenuto dei file varia in base al terzo parametro
	(in base al fatto che dictInput contenga introni filtrati o non)

	Argomenti:
		dictTranscript,	:	dizionario informazioni sui Trascritti
		dictGenes,		:	dizionario informazioni sui Geni
		dictInput,		:	dizionario relativo alle regioni da
							includere nel file '.gtf'
		fileOut,		:	base per il nome dei file di output
		geneNames		:	dizionario con i geni scelti dall'utente

	Ritorna:
		File nel formato '.gtf':
		chrx | start | end | intron/exon_number | gene_name | trascript_name
	"""

	stringaGTF	= 	'%s\t%s\t%s\t%s\t%s\t%s\n'							# Formato della riga da stampare nel file
	exonF		= 	'exon_number "%d"'									# Formato della stringa di tipo exon (True)
	intronF		=	'intron_number "%d"'								# Formato della stringa di tipo intron (False)
	
	# Indici all'interno del dizionario dei transcript
	idx_transcriptName = 0
	idx_geneID         = 1
	
	# Indici all'interno del dizionari dei geni
	idx_geneName  = 0
	idx_cromosoma = 1

	# Indici all'interno del dizionario degli introni e degli esoni
	idx_start = 0
	idx_end   = 1
	idx_tipo  = 2	

	# Tipo di regioni
	esone   = True
	introne = False


	# Apertura e preparazione dei file da scrivere (un file gtf con
	# esoni/introni per ogni gene e uno totale con tutte le regioni per tutti
	# i geni passati dall'utente
	files = {}															  	
	for gene in geneNames:												  
		cod = geneNames[gene]
		
		# Avendo tanti geni, ad ogni nome di gene si associa la relativa
		# cartella del gene corrente tra quelli passati dall'utente
		if not path.exists(cartella % cod):
			system('mkdir ' + cartella % cod)
		files[gene] = open(str(cartella % cod + fileOut), 'w')
		
		
	# File contenente le regioni esoniche/introniche di tutti i geni
	# passati dall'utente (serve per allineamento)
	fileGtf = open(str(fileOut), 'w')							  

	for transcriptID in dictInput:
		geneID 			= dictTranscript[transcriptID][idx_geneID]
		cromosoma		= dictGenes[geneID][idx_cromosoma]
		geneName		= dictGenes[geneID][idx_geneName]
		transcriptName 	= dictTranscript[transcriptID][idx_transcriptName]
		
		# Inizializzazione del numero di esone/introne da stampare nel file
		nrEs 			= 1
		nrIn 			= 1
		
		for i in range(0, len(dictInput[transcriptID][idx_start])):
			start		= dictInput[transcriptID][idx_start][i]
			end			= dictInput[transcriptID][idx_end][i]
			tipo		= dictInput[transcriptID][idx_tipo][i]

			if tipo == esone:
				regione = exonF % (nrEs)								# Stampa della stringa in formato exon
				nrEs   += 1
				
			else:
				regione = intronF % (nrIn)								# Stampa della stringa in formato intron
				nrIn   += 1
				
			strGtf = stringaGTF % (cromosoma, str(start), str(end), 
								   regione, geneName, transcriptName)	# Creazione della riga del file
			
			if geneName in geneNames:									# Se il gene presenta regioni introniche..
				files[geneName].write(strGtf)							# ..si stampa il file gtf relativo alle proprie..
																		# ..regioni introniche nella propria cartella
			fileGtf.write(strGtf)
				
				
	if geneNames:
		for gene in files:
			files[gene].close()

	fileGtf.close()
