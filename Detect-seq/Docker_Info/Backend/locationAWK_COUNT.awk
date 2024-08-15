BEGIN {
    FS = " "
    start = ARGV[1]   # Primo argomento passato in input
    end = ARGV[2]     # Secondo argomento passato in input
    chromosome = ARGV[3]   # Terzo argomento passato in input
    name= ARGV[4]
    ARGV[1] = ARGV[2] = ARGV[3] = ARGV[4] = ""   # Rimuovi gli argomenti dalla lista degli argomenti
    sumCol = 0
    countNonZero= 0
    positions= ""
}

{
    # Estrae il numero del cromosoma dalla prima colonna nel formato "chr12_114438651"
    split($1, chrPos, "_")
    split(chrPos[1], chr, "chr")
    pos = substr(chrPos[2], 1, length(chrPos[2]))
    # Verifica se la posizione e il cromosoma rientrano nell'intervallo specificato
    if (chr[2] == chromosome && pos >= start && pos <= end) {
        positions = positions $1 "\n"   # Aggiungi la posizione all'array
        #print($1)
        #print($2)
        #print($3)
        #print("finedellaprimaparte")
        sumCol += $8
        if ($8 != 0) {
            countNonZero++
        }

    }

}
END {

    outputFileName = FILENAME
    sub(".pmat", "", outputFileName)
    outputFileName2 = outputFileName "_" name "_STAT.txt"
    outputFileName3 = outputFileName "_" name "_STAT_POS.txt"
    print "total number of mutation:", sumCol > outputFileName2
    print "# location with mutation :", countNonZero >> outputFileName2
    print positions > outputFileName3   # Stampa l'array nel file

}
