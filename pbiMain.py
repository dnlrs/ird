import sys
from os import path
from os import system
from string import * 
import subprocess
# Importare funzioni relative al progetto A e al progetto B
# dai file pbiA.py e pbiB.py ripettivamente
import pbiA
import pbiB

"""
PROGETTO 1 e 2 	-------------------------------------------------------
							INTRON RETENTION
 --- Strutture usate --------------------------------------------------

 - DIZIONARIO DEI TRANCRIPT {dictTranscript}

 lo stesso tr_id non appartiene a transcript_name diversi
 lo stesso tr_id non appartiene a geni_id diversi
 lo stesso tr_id non appartiene a geni_name diversi
 lo stesso tr_id non appartiene a cromosomi diversi

	key 	= transcript_id
	value	= [transcript_name, gene_id]



 - DIZIONARIO DEI GENI {dictGenes}

 lo stesso gene_id non appartiene a due gene_name diversi
 lo stesso gene_id non appartiene a due cromosomi diversi

	key 	= 	gene_id
	value	= 	[gene_name, chromosome, [transcript_id]]



 - DIZIONARIO DEGLI ESONI {dictEsoni}

 lo stesso transcript_id ha piu' esoni al suo interno
 nello stesso transcript_id lo strand non cambia
	
 	key 	=	transcript_id
	value	=	[[start_esoni], [end_esoni], strand_transcript_id]



 - DIZIONARIO DEGLI INTRONI {dictIntroni}

 lo stesso transcript_id ha piu' introni al suo interno

	key		= 	transcript_id
	value	=	[[start_introni], [end_introni]]



 - DIZIONARIO DEGLI INTRONI FILTRATI {dictIntroniFiltrati}

 per ogni gene, per ogni transcript_id del gene, si confronta ogni 
 introne del transcript_id con tutti gli esoni del gene e si cercano
 gli introni effettivi per il gene (si eliminano le sovrapposizioni
 degli introni con gli esoni)

	key 	= 	transcript_id
	value 	=	[[start_introni_filtrati], [end_introni_filtrati]]



 - DIZIONARIO DEGLI ESTREMI DEL GENE PER OGNI CROMOSOMA {dictGeneChr}

 lo stesso gene_name puo' far parte di piu' cromosomi
 lo stesso gene_name all'interno di un cromosoma ha piu' gene_id
 
	key		= 	(gene_name + '\t' + chromosome)
	value	=	[start_gene, end_gene]



 - DIZIONARIO DELLE READS {dictReads}

 	key		=	(start_introne + '\t' + end_introne)
	value	=	numero di reads per l'intervallo intronico



 - DIZIONARIO PER LA STRUTTURA DEL FILE DI OUTPUT FINALE {dictOutput}

	(ha un dizionario interno)

	key 	= 	(gene_name + '\t' + chromosome)
	value	=	{	key		=	(start_introne + '\t' + end_introne)
					value	=	[lista di transcript_name per
								 il gene e il cromosoma della chiave
								 esterna]
				}


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 CALCOLO DELLE COORDINATE INTRONICHE 


 Calcolo dell'intervallo intronico se lo strand e' positivo:

   		   start  end
 esone1		2      5
 esone2		9     11
 esone3		15    18

 introni da [fine_esone1+1, inizio_esone2-1]
 regioni introniche:
 introne1   6     8
 introne2   12   14


 Calcolo dell'intervallo intronico se lo strand e' negativo:

		 	   start  end
 esone1		15    18
 esone2		9     11
 esone3		2      5

 In questo caso, si inverte l'ordine degli esoni e si ricade nel
 caso precedente.

 Condizione per avere un introne: tra fine di un esone e inizio 
 del successivo ci sono almeno due basi.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Inizializzazioni dei nomi usati per nominare i file -----------------

pA					=	'pA_'											# Indica file del progetto A
pB					=	'pB_'											# Indica file del progetto B
cartella 			= 	'%s/'											# Il nome corrisponde al nome del gene (solo progetto B)
																		# (aggiunto dopo)
																				
fileIntr			= 	'introni.gtf'									# Nome file gtf con tutti gli introni
fileIntrFilt		=	'introniFiltrati.gtf'							# Nome file gtf con gli introni filtrati
fileEsIn			=	'introniEsoni.gtf'								# Nome file gtf di introni ed esoni(serve per BedGraph)
fileEsInFilt		=	'introniEsoniFiltrati.gtf'						# Nome file gtf di introni filtrati ed esoni
																		# (serve per BedGraph)

reads				=	'reads_mapped.txt'								# Nome file dopo mappatura reads su introni
readsFilt			=	'reads_mapped_filtrati.txt'						# Nome file dopo mappatura reads su introni filtrati
readsEsIn			=	'reads_mapped_introniEsoni.txt'					# Nome file dopo mappatura reads su introni ed esoni
readsEsInFilt		=	'reads_mapped_introniEsoniFiltrati.txt'			# Nome file dopo mappatura reads su introni filtrati
																		# ed esoni

output				=	'output.txt'									# Nome file di output collassato su introni
outputFilt			=	'outputFiltrato.txt'							# Nome file di output collassato su introni filtrati

fileBedGraph		=	''												# il nome corrisponde al cromosoma (aggiunto dopo)
fileBedGraphFilt	=	'_filtrati' 									# il nome corrisponde al cromosoma (aggiunto dopo) con
																		# indicazione di considerare introni filtrati


stringaGTF 			= 	'%s\t%s\t%s\tintron_number "%d"\t%s\t%s\n'		# Formato della stringa da inserire nel file '.gtf'
coverageBed			=	'bedtools coverage -abam %s -b %s > %s'			# Comando di subprocess da richiamare per mappare le
																		# reads sul file fornito come secondo parametro (.gtf)

#----------------------------------------------------------------------

# NEL FILE pbiMain.py sono presenti solo le funzioni comuni ad
# entrambi i progetti A e B. 

def stampaGTF(dictTranscript, dictGenes, dictInput, fileOut, geneNames={}):
	"""	Stampa file in formato '.gtf'.
	
	Il contenuto dei file varia in base al terzo parametro
	(introni filtrati o non)

	Argomenti:
		dictTranscript	:	dizionario informazioni sui Trascritti
		dictGenes		:	dizionario informazioni sui Geni
		dictInput		:	dizionario relativo alle regioni da
							includere nel file '.gtf'
		fileOut			:	base per il nome dei file di output
		geneNames		:	dizionario con i geni scelti dall'utente

	Ritorna:
		File nel formato '.gtf':
		chrx | start | end | intron_number* | gene_name | trascript_name
	"""
	# Indici all'interno del dizionario dei transcript
	idx_transcriptName = 0
	idx_geneID         = 1
	
	# Indici all'interno del dizionari dei geni
	idx_geneName  = 0
	idx_cromosoma = 1

	# Indici all'interno del dizionario degli introni e degli esoni
	idx_start = 0
	idx_end   = 1


	# Apertura e preparazione dei file da scrivere
	if geneNames:																# Se geneNames non e' vuoto allora stiamo eseguendo..
		files = {}															  	# ..il progetto B
		for gene in geneNames:												  
			cod = geneNames[gene]												# Nome del gene inserito dall'utente
			# Avendo tanti geni, ad ogni nome di gene si associa la 
			# relativa cartella contenente i file di output  
			#
			if not path.exists(cartella % cod):									# Se la cartella per il gene corrente non esiste..
				system('mkdir ' + cartella % cod)								# ..si crea..
			files[gene] = open(str(cartella % cod + fileOut), 'w')				# ..altrimenti si apre nella cartella il file in write
		
	fileGtf = open(str(fileOut), 'w')											# Si crea comunque il file '.gtf' totale per tutti i
																				# geni (da usare per mappare le reads)

	for transcriptID in dictInput:
		geneID 			= dictTranscript[transcriptID][idx_geneID]
		cromosoma		= dictGenes[geneID][idx_cromosoma]
		geneName		= dictGenes[geneID][idx_geneName]
		transcriptName 	= dictTranscript[transcriptID][idx_transcriptName]
		
		for i in range(0, len(dictInput[transcriptID][idx_start])):				# Per tutte le regioni introniche trovate
			start		= dictInput[transcriptID][idx_start][i]
			end			= dictInput[transcriptID][idx_end][i]
			
			strGtf = stringaGTF % (cromosoma, str(start), str(end), int(i+1),	# Creazione della stringa in formato '.gtf'
								   geneName, transcriptName)
			
			if geneNames:														# Se si sta eseguendo il progetto B..
				files[geneName].write(strGtf)									# si stampano file separati per ogni gene

			fileGtf.write(strGtf)												# Per progetto A si stampa il file con tutti i geni
																				# Per progetto B si stampa il file totale dei geni..
																				# ..inseriti dall'utente
	if geneNames:
		for gene in files:
			files[gene].close()

	fileGtf.close()

#  --------------------------------------------------------------------

def filtraIntroni(dictTranscript, dictGenes, dictEsoni, dictIntroni):
	""" Filtra le false regioni introniche nei diversi Trascritti.

	Il filtraggio viene fatto per gene,	eliminando i tratti in cui
	gli introni si sovrappongono ad esoni.

	Argomenti:
		dictTranscript	:	dizionario informazioni sui Trascritti
		dictGenes		:	dizionario informazioni sui Geni
		dictEsoni		:	dizionario Esoni per TrascrittoID
		dictIntroni    	: 	dizionario Introni per TrascrittoID

	Ritorna:
		dictIntroniFiltrati	: dizionario Introni filtrati per TrID
	"""
	# Indici all'interno del dizionario degli introni e degli esoni
	idx_starts 	= 0
	idx_ends 	= 1
	
	# Solo per il dizionario degli esoni abbiamo lo strand
	idx_strand 	= 2
	
	# Indici all'interno del dizionario dei Geni
	idx_transcripts = 2
	
	dictIntroniFiltrati	= {}
	for geneID in dictGenes:
		for transcriptID in dictGenes[geneID][idx_transcripts]:
			
			if dictIntroni.has_key(transcriptID):
				numeroIntroni = len(dictIntroni[transcriptID][idx_starts])
				
				# si considera lo strand delle regioni introniche
				if dictEsoni[transcriptID][idx_strand] == '+':
					strandIntroni = True  
				else:
					strandIntroni = False

					
				starts_introni_filtrati = []							      	# Liste introni filtrati
				ends_introni_filtrati = []
				starts_introni = dictIntroni[transcriptID][idx_starts]	      	# Liste introni non filtrati per transcriptID corrente
				ends_introni = dictIntroni[transcriptID][idx_ends]

				
				i = 0
				flagSplitted = False									      	# Flag che indica se l'introne e' splittato da un esone 
				while i < numeroIntroni:								      	# Per tutti gli introni relativi ad un transcript_id 
					flagIntrOK = True									      	# Flag che indica se inserire l'introne nel file finale 
	
					start_introne = int(starts_introni[i])				      	# Inizializzazione delle coordinate..
					end_introne   = int(ends_introni[i])				      	# ..dell'introne corrente
					
					starts_esoni = []									      	# Inizializzazione delle liste contenenti gli esoni..					 
					ends_esoni   = []									      	# ..con i quali l'introne deve essere confrontato
					for innerTranscriptID in dictGenes[geneID][idx_transcripts]:
						
						# si considera lo strand delle regioni esoniche
						if dictEsoni[innerTranscriptID][idx_strand] == '+':
							strandEsoni = True
							starts_esoni = dictEsoni[innerTranscriptID][idx_starts] # Inizializzazione delle coordinate..
							ends_esoni 	 = dictEsoni[innerTranscriptID][idx_ends]	# ..dell'esone corrente
							
						else:														# Se strand e' negativo si inverte l'ordine degli esoni
							strandEsoni = False
							starts_esoni = dictEsoni[innerTranscriptID][idx_starts][::-1]		
							ends_esoni 	 = dictEsoni[innerTranscriptID][idx_ends][::-1]						      
							
						
						start_introne_precedente = 0					      	# Varibili che servono per quando un esone..
						end_introne_precedente   = 0						  	# ..e' completamente compreso nell'introne

						#	Casi di interazione esone-introne analizzati:
						# 
						#	Caso 1: esone				    |---------------|
						# 		 	introne |--------------|
						#
						#
						#	Caso 2: esone			|---------------|
						#			introne		|------------------------|
						#
						#
						# 	Caso 3: esone			|---------------|
						#			introne				|----------------|
						#
						#
						# 	Caso 4: esone			|---------------|
						#			introne		|-----------|
						#
						#
						# 	Caso 5: esone 			|---------------------|
						#			introne					|-------|

						# Legenda parametri commenti:
						#
						# 	fi = fine introne 
						# 	fe = fine esone
						# 	ii = inizio introne 
						# 	ie = inizio esone
						#

						for j in range(0, len(starts_esoni)):
							
							if int(start_introne) < int(ends_esoni[j]):			# Se ii > fe: l'introne non interseca l'attuale esone
																				# Se il contrario si continuano controlli

																				# Caso 1
								if int(end_introne) < int(starts_esoni[j]):	    # Essendo ii < fi, se ii < fe e se anche fi < ie:
									break										# l'introne si trova tra l'esone precedente
																				# e l'esone corrente e non interseca nulla

																					   # Caso 2
								else: 												   # A questo punto: ii < fe e fi > ie				
									if int(end_introne) > int(ends_esoni[j]):		   # Se fi > fe:
										if int(start_introne) < int(starts_esoni[j]):  # Se anche ii < ie: esone contenuto in introne
																					   # L'introne viene splittato in due parti:
																					   # -una a sinistra dell'esone(introne_precedente)
																					   # -una alla destra dell'esone
											start_introne_precedente = int(start_introne)	
											end_introne_precedente   = int(starts_esoni[j]) - 1
											start_introne            = int(ends_esoni[j]) + 1	# parte alla destra dell'esone
											# end_introne										# end_introne non viene aggiornato

																								# Caso 3
										else:													# Se fi > fe e anche ii > ie:
																								# l'introne inizia dopo l'esone
																								# (ma finisce dopo)
											start_introne            = int(ends_esoni[j]) + 1	# Rimane la parte a destra dell'esone

																								# Caso 4
									else:														# Se fi < fe:
										if int(start_introne) < int(starts_esoni[j]):			# Se ii < i:
											end_introne = int(starts_esoni[j]) - 1				# Rimane la parte a sinista dell'esone	

																								# Caso 5
										else:													# Invece se ii > ie (e fi < fe):
											flagIntrOK = False									# introne non viene considerato		
							
							if int(start_introne_precedente) != 0:				# Se l'introne e' stato splittato dall'esone
								flagSplitted = True
								numeroIntroni += 1
								starts_introni.append(start_introne_precedente)	# Si appende alla lista degli introni da analizzare..
								ends_introni.append(end_introne_precedente)		# ..la parte precedente dell'introne (quella successiva..
								start_introne_precedente = 0					# ..continua ad essere elaborata 
								end_introne_precedente   = 0

	
					if flagIntrOK:												
						starts_introni_filtrati.append(start_introne)			# Si aggiunge l'introne alla lista di introni..
						ends_introni_filtrati.append(end_introne)				# ..da inserire nel file degli introni filtrati
					
					i += 1
	
	
				numeroIntroniFiltrati = len(starts_introni_filtrati)
				if numeroIntroniFiltrati > 0:
					if flagSplitted:										  	# Si puo' ordinare perche'..
						starts_introni_filtrati.sort()						  	# ..gli introni non si sovrappongono per lo stesso..
						ends_introni_filtrati.sort()						  	# ..trascritto: rimane la corrispondenza start-end..
																				# ..ordinando separatamente le liste
					if not strandIntroni:
						starts_introni_filtrati.reverse()					  	# Si mantiene l'informazione sullo strand in base..
						ends_introni_filtrati.reverse()						  	# ..allo strand del trascritto di origine dell'introne
					
					dictIntroniFiltrati[transcriptID] = [starts_introni_filtrati, \
														 ends_introni_filtrati]
		
	return dictIntroniFiltrati

#  --------------------------------------------------------------------

def estraiReads(fileInput, fileBam, fileOut, geneNames={}):
	""" Esegue la chiamata del tool di mappatura delle reads su file gtt:
		tool Coverage di bamtools
		
	A seconda dei file che vengono passati come parametri, la funzione
	esegue la mappatura su regioni introniche filtrate, non filtrate e
	sull'intera sequenza del gene.
	
	Argomenti:
		fileInput		:		nome del file gtf in in input
		fileBam			:		file di allineamento
		fileOut			:		file di output con le reads mappate
	"""

	subprocess.call(coverageBed	% (fileBam, fileInput, fileOut), shell=True)

#  --------------------------------------------------------------------

def caricaReads(fileInput):
	"""	Crea il dizionario delle reads dal file di input.

	Argomenti:
		fileInput		:		file di mappatura delle reads

	Ritorna:
		dictReads		:		dizionario delle reads
	"""
	# Nel file di mappatura reads di tipo '.gtf' le reads sono in
	# posizione idx_reads
	idx_start	= 1
	idx_end		= 2
	idx_reads 	= 6

	dictReads  	= {}
	keyF = '%s\t%s'														# Formato delle key del dizionario dictReads
	
	for x in open(fileInput):
		riga = x.strip('\n').split('\t')
		keyEstremi = keyF % (riga[idx_start], riga[idx_end])
	
		if int(riga[idx_reads]) > 0:									# Si aggiunge l'introne solo se il numero di..
			dictReads[keyEstremi] = riga[idx_reads]						# ..reads mappate e' maggiore di 0


	return dictReads

#  --------------------------------------------------------------------

def collapsed(dictTranscript, dictGenes, dictIntroni, dictGeneChr, fileInput, fileOut):
	"""	Esegue le elaborazioni necessarie per il file di output.

	Argomenti:
		dictTranscript,	:	dizionario informazioni sui Trascritti
		dictGenes,		:	dizionario informazioni sui Geni
		dictIntroni,   	:	dizionario Introni per TrascrittoID
		dictGeneChr		:	dizionario coordinate Geni per Cromosomi

	Ritorna:
		Valore booleano che indica se sono state trovate regioni con
		intron retention nei geni analizzati.
	"""
	dictOutput = {}
	keyF = '%s\t%s'														# Formato delle key per il dizionario dictOutput

	dictReads  = caricaReads(fileInput)

	if not dictReads :
		print 'Non sono state trovate reads mappanti per i geni inseriti.'
		return False

	#   Creazione della struttura dati di output {dictOutput}
	#	(con un dizionario interno)
	#
	#	key 	= 	(gene_name+ '\t' + chromosome)
	#	value	=	{	
	#					key		=	(start_introne+ '\t' + end_introne)
	#					value	=	[lista di transcript_name per
	#								 il gene e il cromosoma della
	#								 chiave esterna]
	#				}
	#
	# Accumulo per ogni introne dei trascritti dove esso e' presente.
	#
	for transcriptID in dictIntroni:									# Per ogni transcript_id del file di annotazione..
																		# ..se ha introni..
		transcriptName = dictTranscript[transcriptID][0]				# ..si risale al suo gene_id (gene_name) e.. 
		key_geneID     = dictTranscript[transcriptID][1]				# ..al cromosoma univoci (dai dizionari creati..
																		# ..durante l'inizializzazione)
		key_geneChr    = keyF % (dictGenes[key_geneID][0], \
								 dictGenes[key_geneID][1])
		
		for i in range(0, len(dictIntroni[transcriptID][0])):			  	# Per ogni introne del transcript_id corrente
			key_estremi = keyF % (str(dictIntroni[transcriptID][0][i]), \
								  str(dictIntroni[transcriptID][1][i]))	  	# Si crea la chiave sull'intervallo 
			
			if dictReads.has_key(key_estremi):	
						
				# Se nell'intervallo corrente sono state mappate 	
				# reads si aggiunge l'intervallo nel dizionario 
				# dictOutput nella posizione corrispondente al 
				# gene_name e cromosoma (relativi al transcript_id 
				# al quale l'intervallo stesso appartiene)
				#
				if not dictOutput.has_key(key_geneChr):			
					dictOutput[key_geneChr] = {}					
					dictOutput[key_geneChr][key_estremi] = [transcriptName]
					
				else:
					if not dictOutput[key_geneChr].has_key(key_estremi):
						dictOutput[key_geneChr][key_estremi] = [transcriptName]
						
					else:
						if transcriptName not in dictOutput[key_geneChr][key_estremi]:
							dictOutput[key_geneChr][key_estremi].append(transcriptName)


	# Richiamo della funzione per la stampa del file di output
	stampaOutput(dictGeneChr, dictReads, dictOutput, fileOut)
	return True

# ---------------------------------------------------------------------

def stampaOutput(dictGeneChr, dictReads, dictOutput, fileOut):
	"""	Stampa il file di output collassato sugli introni con reads.

	Formato file di output:
	gene_name(chr'':start_gene-end_gene)	[intron#(start_intron-end_intron):[transcripts]:#reads]
	
	Argomenti:
		dictGeneChr, 		: 		dizionario coordinate Geni per Cromosomi
		dictReads,			:		dizionario delle reads
		dictOutput			:		dizionario con struttura dati finale 
		fileOut				:		file di output finale
	"""
	GeneChrF = '%s(%s:%s)'												# Formato prima parte della riga di output(gene e chr)
	IntroneF = '\tintron%d(%s-%s):%s:%s'								# Formato seconda parte della riga di output..
																		# ..(introni, trascritti e reads)
	output = open(fileOut, 'w')
	
	for geneChr in dictOutput:											# Per ogni combinazione "gene_name cromosoma"
		geneName, chromosome = geneChr.split('\t')	
		strGene = GeneChrF % (geneName, chromosome, \
							  '-'.join(dictGeneChr[geneChr]))
		
		output.write(strGene)											# Si scrive la prima parte della riga relativa..
																		# ..al gene e cromosoma correnti
																				
		intronNumber = 1												# Contatore per gli introni appartenenti al gene..
																		# ..e al cromosoma correnti
		for intron in sorted(dictOutput[geneChr].keys()):
			strTranscripts = ', '.join(dictOutput[geneChr][intron])		# Si crea la stringa dei trascritti_name per..
			startIntrone, endIntrone = intron.split('\t')				# ..l'intervallo intronico corrente
			nReads = dictReads[intron]
			
			strIntrone = IntroneF %(int(intronNumber), startIntrone, \
									endIntrone, strTranscripts, nReads)
			
			output.write(strIntrone)
			intronNumber += 1

		
		output.write('\n')

	
	output.close()

# -------------------------------------------------------------------------

def progettoA(fileAnn, fileBam, scelta):
	"""	Esegue il progetto A.

	Argomenti:
		fileAnn		:	file di annotazione
		fileBam		:	file di alineamento
		scelta		:	scelta effettuata dall'utente
	"""
	# Strutture dati
	#-----------------------
	dictTranscript 		= {}
	dictGenes 			= {}
	dictEsoni 			= {}
	dictIntroni 		= {}
	dictIntroniFiltrati = {}
	dictGeneChr 		= {}
	#-----------------------
	
	system('clear')
	""" Inizializzazione """
	print '\nInizializzazione delle strutture dati',
	print 'per il progetto A in corso...'
	dictTranscript, \
	dictGenes, \
	dictEsoni, \
	dictIntroni, \
	dictGeneChr = pbiA.inizializzazione(fileAnn)								

	if int(scelta) == 1:
		dictI 	= dictIntroni
		fileI 	= fileIntr
		readsI 	= reads
		outputI	= output

	elif int(scelta) == 2:
		# Per ogni funzione il dizionario di input degli introni
		# e' quello filtrato ottenuto con la funzione seguente:
		""" Filtraggio degli introni """
		print 'Filtraggio degli introni in corso...'
		dictIntroniFiltrati = filtraIntroni(dictTranscript, dictGenes, \
											dictEsoni, dictIntroni)
		print 'Filtraggio completato.\n'

		dictI 	= dictIntroniFiltrati
		fileI 	= fileIntrFilt
		readsI 	= readsFilt
		outputI	= outputFilt

	# Operazioni comuni alla scelta 1) e 2)
	#
	""" Stampa del file '.gtf' """
	print 'Stampa del file \'.gtf\' in corso...' 
	stampaGTF(dictTranscript, dictGenes, dictI, pA + fileI)
	print 'E\' stato creato il file %s.\n' % (pA + fileI)

	""" Esecuzione della mappatura reads"""
	print 'Mappatura delle reads in corso...'
	estraiReads(pA + fileI, fileBam, pA + readsI)
	print 'E\' stato creato il file %s.\n' % (pA + readsI)

	""" Creazione del file di output """
	print 'Creazione del file di output in corso...'
	collapsed(dictTranscript, dictGenes, dictI, dictGeneChr, \
			  pA + readsI, pA + outputI)
	print 'E\' stato creato il file %s.\n' % (pA + outputI)

	# Pulizia memoria
	#-----------------------
	del dictTranscript		
	del dictGenes
	del dictEsoni
	del dictIntroni
	del dictIntroniFiltrati
	del dictGeneChr
	#-----------------------

# ---------------------------------------------------------------------

def progettoB(fileAnn, fileBam, scelta):
	"""
	Esegue il progetto B.

	Argomenti:
		fileAnn		:	file di annotazione
		fileBam		:	file di alineamento
		scelta		:	scelta effettuata dall'utente
	"""
	# Strutture dati
	#-----------------------
	dictTranscript 		= {}
	dictGenes 			= {}
	dictEsoni 			= {}
	dictIntroni 		= {}
	dictIntroniFiltrati = {}
	dictEsIn			= {}
	dictGeneChr 		= {}
	geneNames 			= {}
	dictReadsEsIn		= {}
	#-----------------------
	
	system('clear')
	
	""" Acquisizione dei geni da elaborare """
	print 'Inserisci la lista dei nomi dei geni',
	print 'separati da spazio:\n',
	geniInseriti = raw_input('>> ').strip('\n').split(' ')

	nomeGeneF = 'gene_name "%s"'										# Formato del nome del gene inserito dall'utente

	# Questo controllo evita interpretazioni errate dell'input
	# nel caso di uno spazio alla fine della riga o di un doppio spazio
	for gene in geniInseriti:
		if gene:
			geneNames[nomeGeneF % gene] = gene

	""" Inizializzazione """
	print 'Inizializzazione delle strutture dati',
	print 'per il progetto B in corso...\n'
	dictTranscript, \
	dictGenes, \
	dictEsoni, \
	dictIntroni, \
	dictGeneChr = pbiB.inizializzazione(fileAnn, geneNames)

	# Se dopo l'inizializzazione con l'eliminazione dei geni che
	# non presentano regioni introniche non rimangono geni da
	# analizzare, si esce dal progetto B.
	#
	if not geneNames:
		print 'Attenzione: nessuno dei geni inseriti presenta',
		print 'regioni introniche.'
		return None


	if int(scelta) == 1 or int(scelta) == 3:
		dictI 	= dictIntroni
		fileI 	= fileIntr
		readsI 	= reads
		outputI	= output


	elif int(scelta) == 2 or int(scelta) == 4:
		# Per ogni funzione il dizionario di input degli introni
		# e' quello filtrato ottenuto con la funzione seguente:
		""" Filtraggio degli introni """
		print 'Filtraggio degli introni in corso...'
		dictIntroniFiltrati = filtraIntroni(dictTranscript, dictGenes, \
											dictEsoni, dictIntroni)
		print 'Filtraggio completato.\n'
		
		dictI 	= dictIntroniFiltrati
		fileI 	= fileIntrFilt
		readsI 	= readsFilt
		outputI	= outputFilt 


	""" Procedimento:
			- si stampano N file GTF, uno per ogni gene, e
			  contemporaneamente se se ne stampa uno completo con tutti
			  i geni
			- si richiama bedtools per quest'unico file
			- si stampa il file di output 
	"""
	""" Stampa del file '.gtf' """
	print 'Stampa dei file \'.gtf\' per ogni gene in corso...'
	stampaGTF(dictTranscript, dictGenes, dictI, pB + fileI, geneNames)
	print 'Sono stati creati i file %s' % (pB + fileI),
	print 'in ciascuna cartella dei geni.\n'


	""" Esecuzione della mappatura reads """
	print 'Mappatura delle reads in corso...'
	estraiReads(pB + fileI, fileBam, pB + readsI)
	print 'E\' stato creato il file %s.\n' % (pB + readsI)


	""" Creazione del file di output """
	print 'Creazione del file di output in corso...'
	collapsed(dictTranscript, dictGenes, dictI, dictGeneChr, \
			  pB + readsI, pB + outputI)
	print 'E\' stato creato il file %s.\n' % (pB + outputI)



	if int(scelta) == 3:
		fileEsInI = fileEsIn
		readsEsInI = readsEsIn
		fileBedGraphI = fileBedGraph

	elif int(scelta) == 4:
		fileEsInI = fileEsInFilt
		readsEsInI = readsEsInFilt
		fileBedGraphI = fileBedGraphFilt
		
	else:
		# Pulizia memoria
		#----------------------
		del dictTranscript
		del dictGenes
		del dictEsoni
		del dictIntroni
		del dictIntroniFiltrati
		del dictEsIn
		del dictGeneChr
		del geneNames 
		del dictReadsEsIn
		#----------------------
		return



	""" Ricombinazione delle regioni introniche ed esoniche """
	dictEsIn = pbiB.unisciRegioni(dictEsoni, dictI)


	""" Stampa del file '.gtf' con introni ed esoni e tutti i geni """
	print 'Stampa del file \'.gtf\' con introni',
	print 'ed esoni di tutti i geni in corso...' 
	pbiB.stampaGTFEsIn(dictTranscript, dictGenes, dictEsIn, \
					   pB + fileEsInI, geneNames)
	print 'E\' stato creato il file %s.\n' % (pB + fileEsInI)


	""" Esecuzione della mappatura reads esoni + introni """
	print 'Mappatura delle reads su regioni introniche ed esoniche in corso...'
	estraiReads(pB + fileEsInI, fileBam, pB + readsEsInI)
	print 'E\' stato creato il file %s.\n' % (pB + readsEsInI)


	""" Caricamento in memoria delle reads esoni + introni """
	dictReadsEsIn = pbiB.caricaReadsEsIn(pB + readsEsInI)


	""" Stampa dei file bedGraph """
	print 'Stampa dei file \'.bedGraph\' in corso...'
	pbiB.stampaBedGraph(dictReadsEsIn, dictGeneChr, \
						fileBedGraphI, geneNames)
	print 'Sono stati creati i file \'.bedGraph\' per ogni gene.\n' 


	# Pulizia meoria
	#-----------------------
	del dictTranscript
	del dictGenes
	del dictEsoni
	del dictIntroni
	del dictIntroniFiltrati
	del dictEsIn
	del dictGeneChr
	del geneNames 
	del dictReadsEsIn
	#------------------------


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Menu iniziale nel quale scegli il progetto
#
# L'utenete puo' scegliere progetto A per ottenere il file di output 
# per tutti gli introni derivanti dal file di annotazione
# oppure il progetto B se vuole analizzare solo determinati geni
# e/o ottenere i file in formato BedGraph.
#
def menuProgetto():
	while True:
		print """Scegli progetto:
		
	a) Progetto A
	b) Progetto B

	e) Esci
"""
		scelta = raw_input(">> ").lower()
		if scelta in ['a', 'b', 'e']:
			return scelta
		else:
			system('clear')
			print "Scelta progetto: Scelta non valida, riprovare.\n"

# ---------------------------------------------------------------------

def menuA():
	print """Azioni per il progetto A:

	1) Esegui il progetto A con tutti gli introni
	2) Esegui il progetto A con gli introni filtrati

		
	c) Cambia progetto
	e) Esci
	"""
	scelta = raw_input(">> ")
	
	return scelta

# ---------------------------------------------------------------------

def menuB():
	print """Azioni per il progetto B:
		
	1) Esegui il progetto B con tutti gli introni
	2) Esegui il progetto B con gli introni filtrati

	3) Punto 1 + bedGraph
	4) punto 2 + bedGraph

		 
	c) Cambia progetto
	e) Esci
	"""
	scelta = raw_input(">> ")
	
	return scelta

# ---------------------------------------------------------------------

def inserisciFileDiAnn():
	""" Acquisisce il nome del file di annotazione.

	Se il file non esiste lo richiede finche' non
	ottiene un file valido.
	"""
	system('clear')
	while True:
		fileAnn = raw_input("Inserisci il nome del file di annotazione: ")
		if not path.exists(fileAnn):
			system('clear')
			print 'File non trovato!'
		else:
			return fileAnn

# ---------------------------------------------------------------------

def inserisciFileDiAll():
	""" Acquisisce il nome del file di allineamento (.bam).

	Se il file non esiste lo richiede finche' non
	ottiene un file valido.
	"""
	system('clear')
	while True:
		fileAll = raw_input("Inserisci il nome del file di allineamento (.bam): ")
		if not path.exists(fileAll):
			system('clear')
			print 'File non trovato!'
		else:
			return fileAll

# -------------------------------------------------------------------------

def main():

	system('clear')
	# Il file di annotazione viene richiesto all'utente							
	fileAnn 	= inserisciFileDiAnn()
	
	# Il file di allineamento viene richiesto all'utente
	fileBam 	= inserisciFileDiAll()

	system('clear')
	finito = False
	while not finito:
		progetto = menuProgetto()
		system('clear')

		if progetto == 'a':
			scelta = menuA()
			if scelta == '1' or scelta == '2':
				progettoA(fileAnn, fileBam, scelta)

			elif scelta.lower() == 'c':
				continue
				
			elif scelta.lower() == 'e':
				finito = True
				
			else:
				system('clear')
				print 'Scelta \'%s\' non valida per il progetto A!' % scelta


				
		elif progetto == 'b':
			scelta = menuB()
			if scelta in ['1', '2', '3', '4']:
				progettoB(fileAnn, fileBam, scelta)
				
			elif str(scelta).lower() == 'c':
				continue
				
			elif str(scelta).lower() == 'e':
				finito = True

			else:
				system('clear')
				print 'Scelta \'%s\' non valida per il progetto B!' % scelta



		elif progetto == 'e':
			finito = True
			
			
# Identificazione della funzione che costituisce il main
if __name__ == '__main__': main()
