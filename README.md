# Specifiche 

## Funzionalità

Il progetto A deve produrre:
- il file gtf con tutte le regioni introniche
- il file di output collassato sulle regioni introniche 

Il progetto B deve produrre:
- N file gtf degli N geni passati dall'utente che presentano regioni
  introniche (se il gene non presenta regioni introniche non viene 
  stampato il file gtf)
- il file di output (unico) collassato sulle regioni introniche dei geni 
  che sono stati passati dall'utente
  
Funzionalità aggiuntive:
- produzione degli output considerando gli introni filtrati, ossia 
  eliminando dalle regioni introniche quelle parti che in altri 
  trascritti costituiscono regioni esoniche (come falsi positivi)
- solo per progetto B: possibilità di stampare il file BedGraph 
  relativo al gene inserito dall'utente (qualora questo abbia regioni
  introniche perchè il problema è intron retention)
  
## Verifiche del corretto funzionamento del programma sviluppato

- Verifica di aver ottenuto il formato dei file in modo corretto
- Verifica di aver calcolato correttamente le regioni introniche
  (metodo grafico mappatura introni esoni per lo stesso trascritto)
- Verifica di aver effettuato in maniera corretta il filtraggio degli 
  introni (metodo grafico..)
- Verifica di aver ottenuto il file di output solo per introni che hanno
  un numero di reads diverso da 0 e che il file sia ottenuto in maniera
  coerente con i file precedenti
  

## Descrizione dei progetti

Il problema affrontato nei progetti è quello di andare a valutare la 
presenza dell'intron retention relativo ai geni contenuti nel file di 
annotazione fornito (Homo_sapiens.GRCh.37.60.chr.gtf).

Il progetto A si occupa di fornire il file di output che contiene le
informazioni relative al quantitativo di reads che vengono mappate
sui singoli introni nei vari geni per l'intero file di annotazione 
fornito. 

Il progetto B si occupa di fornire il file di output nello stesso 
formato del progetto A ma solo per i geni che l'utente vuole analizzare.
Inoltre, con il file in formato bedGraph per i vari geni, è possibile
visualizzare su genome Browser la mappatura delle reads.

Il problema che si può presentare è la presenza di falsi positivi: 
quando un introne in un altro trascritto è un esone (totalmente o in 
parte) un quantitativo alto di reads potrebbe non essere associato al
fenomeno di intron retention. Per questo motivo è utile produrre la 
stessa tipologia di file di output filtrando gli introni nei vari 
trascritti dello stesso gene. In questo modo le regioni introniche 
rimanenti sono tali per tutti i trascritti e, qualora si verificasse, 
sarebbe oppurtuno parlare di intron retention.



## Guida all'utente

Passi per un'esecuzione corretta del programma:
- Lancia il programmma nominato 'pbiMain.py' da terminale, comparirà una 
  schermata che richiede l'inserimento del nome del file di annotazione 
  e del file di allineamento da utilizzare.
- Inserisci il nome dei file richiesti, comparirà una schermata per 
  procedere con la scelta del progetto (A o B)o uscire dall'applicativo.
- Inserisci la lettera corrispondente all'opzione da eseguire.
- Se la scelta è di eseguire uno dei due progetti comparirà il relativo 
  menù di scelta per ottenere gli output desiderati. Se hai inserito una
  lettera non valida comparirà un opportuno messaggio di errore e la 
  precedente schermata per procedere.


Scelta di procedere con il progetto A:
- Scegli una delle opzioni riportate nel menù del progetto A per 
  proseguire, 'c' per cambiare progetto oppure 'e' per uscire. 
- Se la scelta non è valida comparirà un apposito messaggio di 
  errore e la schermata iniziale di 'scegli progetto' per procedere. 
  Se la scelta è valida comparirà una schermata che riporta l'esecuzione 
  corrente,   indicando i passi che si stanno svolgendo e i file di 
  output creati.
  In particolare, se la scelta è:
  - 1) Esegui il progetto A con tutti gli introni,  
    il programma crea un file 'pA_introni.gtf' per tutte le regioni 
    introniche ricavate dal file di annotazione fornito considerando   
    solo le righe dove compare 'exon' che verrà utilizzato per eseguire 
    l'allineamento delle reads. 
    Da quest'ultimo processo viene stampato un file 'pA_reads_mapped.txt' 
    che contiene le informazioni delle reads mappate sulle varie regioni 
    introniche.
    Il programma procede con la stampa del file 'pA_output.txt' contenete
    le informazioni della mappatura collassate per introni: ogni introne
    è associato ad un gene e ad un cromosoma e ad ogni introne possono
    essere associati più trascritti.
  - 2) Esegui il progetto A con gli introni filtrati,
    il programma crea un file 'pA_introniFiltrati.gtf' per le regioni
    introniche filtrate ricavate dal file di annotazione fornito 
    considerando solo le regioni dove compare 'exon' che verrà utilizzato
    per eseguire l'allineamento delle reads.
    Da quest'ultimo processo viene stampato il file 'pA_reads_mapped_filtrati.txt'
    che contiene le informazioni delle reads mappate sulle varie regioni 
    introniche filtrate.
    Il programma procede con la stampa del file 'pA_outputFiltrato.txt' 
    contenente le informazioni della mappatura collassate per introni
    filtrati: ogni introne è associato ad un gene e ad un cromosoma e ad
    ogni introne possono essere associati più trascritti.
    
    
  Alla fine dell'esecuzione del programma comparirà la schermata 
  iniziale della scelta del progetto.
 
Scelta di procedere con il progetto B:
- Scegli una delle opzioni riportate nel menù del progetto B per 
  proseguire, 'c' per cambiare progetto oppure 'e' per uscire.
- Se la scelta non è valida comparirà un apposito messaggio di 
  errore e la schermata iniziale di 'scegli progetto' per procedere. 
  Se la scelta è valida comparirà una schermata che richiede di inserire 
  il nome dei geni da analizzare.
- Inserisci il nome dei geni (differenziando tra lettere maiuscole e 
  minuscole) separati da uno spazio. Comparirà una schermata che riporta 
  l'esecuzione corrente, indicando i passi che si stanno svolgendo e i
  file di output creati. Se il gene inserito non contiene regioni 
  introniche verrà stampato un messaggio per indicare l'accaduto e
  quel gene viene escluso dalla lista dei geni da elaborare.
  In particolare, se la scelta è:
  - 1) Esegui il progetto B con tutti gli introni,
    dopo aver inserito il nome dei geni da analizzare verranno create 
    N cartelle, tante quanti sono i geni passati da linea di comando che 
    effettivamente presentano regioni introniche. Nelle varie cartelle, 
    per ogni gene da analizzare, verrà creato un file 'pB_introni.gtf' per 
    le regioni introniche ricavate dal file di annotazione fornito 
    considerando solo le regioni dove compare 'exon'.  Nella cartella
    principale verrà creato un unico file 'pB_introni.gtf' che contiene 
    le regioni introniche trovate per ciascun gene e che verrà utilizzato
    per eseguire l'allineamento sulle stesse.
    Nella cartella principale verrà stampato l'output ottenuto con 
    l'allineamento 'reads_mapped.txt'.
    Il programma procede con la stampa del file 'pB_output.txt' contenente 
    le informazioni della mappatura collassate per introni : 
    ogni introne è associato ad un gene e ad un cromosoma e ad
    ogni introne possono essere associati più trascritti.
  - 2) Esegui il progetto B con gli introni filtrati
	dopo aver inserito il nome dei geni da analizzare verranno create 
    N cartelle, tante quanti sono i geni passati da linea di comando che 
    effettivamente presentano regioni introniche. Nelle varie cartelle, 
    per ogni gene da analizzare, verrà creato un file 'pB_introniFiltrati.gtf' per 
    le regioni introniche filtrate ricavate dal file di annotazione fornito 
    considerando solo le regioni dove compare 'exon'.  Nella cartella
    principale verrà creato un unico file 'pB_introniFiltrati.gtf' che contiene 
    le regioni introniche trovate per ciascun gene e che verrà utilizzato
    per eseguire l'allineamento sulle stesse.
    Nella cartella principale verrà stampato l'output ottenuto con 
    l'allineamento 'reads_mapped_filtrati.txt'.
    Il programma procede con la stampa del file 'pB_outputFiltrati.txt' contenente 
    le informazioni della mappatura collassate per introni : 
    ogni introne è associato ad un gene e ad un cromosoma e ad
    ogni introne possono essere associati più trascritti.
  - 3) Scelta 1 + bedGraph,
    si esegue la scelta 1 e, in aggiunta, verranno stampati i file
    'pB_InEs.gtf' per le regioni introniche ed esoniche per ogni 
    gene nella cartella corrispondente al gene in esame e un unico file 
    'pB_introniFiltrati.gtf' nella cartella principale che verrà utilizzato 
    per eseguire l'allineamento.
    Verrà stampato l'output dell'allineamento 'reas_mappedInEs.txt' e i
    file in formato BedGraph per ogni gene nelle rispettive cartelle
    dei geni ('nome_gene.bedGraph').
  - 4) Scelta 2 + bedGraph,
    si esegue la scelta 2 e, in aggiunta, verranno stampati i file
    'pB_introniEsoni_filtrati.gtf' per le regioni introniche filtrate ed esoniche per ogni 
    gene nella cartella corrispondente al gene in esame e un unico file 
    'pB_InESFiltrati.gtf' nella cartella principale che verrà utilizzato 
    per eseguire l'allineamento.
    Verrà stampato l'output dell'allineamento 'reas_mappedInEs_filtrati.txt' e i
    file in formato BedGraph per ogni gene nelle rispettive cartelle
    dei geni ('nome_gene.bedGraph').


## Albero funzioni

```
main()
|
|----> inserisciFileDiAnn()
|----> inserisciFileDiAll()
|
|----> menuProgetto()
		|
		|----> menuA()
		|		|
		|		|----> progettoA()
		|				|
		|				|----> pbiA.inizializzazione()
		|				|		|
		|				|		|----> filtraFileDiAnn()
		|				|
		|				|----> *filtraIntroni()
		|				|----> stampaGTF()
		|				|----> estraiReads()
		|				|		|
		|				|		|----> callBedtools()
		|				|
		|				|----> collapsed()
		|						|
		|						|----> caricaReads()
		|						|----> stampaOutput()
		|
		|
		|----> menuB()
				|
				|---- progettoB()
						|
						|----> pbiB.inizializzazione()
						|		|
						|		|----> filtraFileDiAnn()
						|
						|----> *filtraIntroni()
						|----> stampaGTF()
						|----> estraiReads()
						|		|
						|		|----> callBedtools()
						|
						|----> collapsed()
						|		|
						|		|----> caricaReads()
						|		|----> stampaOutput()
						|
						|----> pbiB.unisciRegioni()
						|----> pbiB.stampaGTFEsIn()
						|----> estraiReads()
						|		|
						|		|----> callBedtools()
						|
						|----> pbiB.caricaReadsEsIn()
						|----> pbiB.stampaBedGraph()

* solamente se si lavora con gli introni filtrati.
```










