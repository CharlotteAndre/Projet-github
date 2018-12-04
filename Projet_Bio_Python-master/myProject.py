"""
22/11/2018
ANDRE Charlotte
CLAUDE Elsa
HUCTEAU Alexis
LAPORTE Antoine
BioPython
"""

#2.2.1 Statistiques basiques


def oneWord(seq,pos,wlen): #return a piece with a lenght = wlen of the string seq from the index start
    '''
    La fonction sert à copier une séquence dans une liste à partir d'une position start et sur une longueur donnée

        argument: wlen: longueur de la liste voulue
        pos: position voulue de départ
        seq:correspond a la séquence Fasta récuperer

        Return: la fonction return une liste appelée sequence contenant
        la séquence nucléotidique débutant a la position strat souhaité et d'une longueur wlen
    '''
    a=0
    sequence=''
    for i in range(len(seq)):
        if i>=pos and a<wlen :
            a=a+1
            sequence=sequence+str(seq[i])
    return sequence



def isCodonStart(seq,pos): #return True if the string seq has start codon
    '''
    La fonction sert à identifier le codon start de la séquence, en parcourant la séquence de 1 en 1 et stocke dans la variable
        codon une suite de 3 nucléotides.

        argument: seq:séquence Fasta
                pos: position de départ du codon de 3 nucléotides a tester

        return: La fonction retourne True si en parcourant la séquence,
        elle trouve un des codons start de la liste, sinon elle retourne False
    '''

    codon=oneWord(seq,pos,3)
    if codon in ('ATG','TTA','TTG','CTG','ATT','ATC','ATA','GTG'):
        return True
    else :
        return False

def isCodonStop(seq,pos): #return True if seq has a stop codon
    '''
    la fonction sert à identifier le codon stop de la séquence, en parcourant la séquence de 1 en 1 et stocke
            dans la variable codon une suite de 3 nucléotides

        argument: seq: Séquence fasta
                pos:position de départ du codon de 3 nucléotides à tester

        return: La fonction return True, si en parcourant la séquence
        elle trouve un des deux codon Stop, sinon elle retourne False
    '''

    codon=oneWord(seq,pos,3)
    if codon=='TAA' or codon=='TAG':
        return True
    else :
        return False


def isGene3(seq,version,threshold=90):
    '''
    La fonction parcours la séquence sur sa longueur -2. Elle cherche alors en appellant la fonction isCodonStart
        si il y a présence d'un codon sans cadre de lecture spécifique. Si elle identifie un codon alors
        elle parcours de 3 en 3 la séquence en commencant à la position i+threshold, qui correspond a la position
        du codon start plus un seuil minimal d'ORF,par défaut cette valeur de threshold est de 90. Une fois cette distance atteinte,
        la fonction va alors appelé isCodonStop afin d'identifier un codon stop dans l'ORF associé au codon start identifié, elle s'arrete au premier des qu'elle en trouve un.
        Une fois un ORF determine, la fonction va print la version de la séquence étudié (5' 3' ou complémentaire inversé),
        le cadre de lecture (0,1,2), la longueur j+2-i égale au début du codon start et à la fin du codon stop. De plus, elle print la suite nucléotidique du codon start et du codon stop
        et leur position, sinon elle retourne False

        arguement: seq: séquence Fasta étudiée
                    version: la version peut être 5'3' ou complementaire inversé
                    threshold: cela correspond à une longueur d'ORF minimal afin de limiter les faux positifs

        return: La fonction ne return rien, mais elle print des informations
    '''
    for i in range(len(seq)-2):
        if isCodonStart(seq,i) == True:
            for j in range(i+threshold,len(seq)-2,3):
                if isCodonStop(seq,j) == True:
                        print("=====================",version,": Frame :",i%3,"===================== ")
                        print("Length:",j+2-i,"pb")
                        print("Codon Start : "+oneWord(seq,i,3)+", Position : ",i,"; Codon Stop : "+oneWord(seq,j,3)+" position : ",j,"\n")
                        break

######### WARNING ###########
#getGeneticCode(NCBI_ID) pour les biais de codon entre les espèces
#findORF(seq, threshold,codeTable) = isGene3 mais prend en compte les biais de Codon


#Ce serait bien de save les séquences inv, comp et comp inv dans un fasta ? DONE
# -> Ajout des codons init' alternatifs
# -> Suppression d'un codon stop (non présent dans le biais de codon)
# -> Ajout Threshold (minimum length : 90pb) "various thresholds (No threshold, 90bp, 210bp, 300bp, 420bp, for example.)"
#############################

def openFasta(file):
    '''La fonction a pour objectif de récuperer le fichier Fasta présent dans le repoertoire
    courant, elle le converti en string, en enlevant les sauts de lignes et le stocke dans une variable

    argument: file: nom du fichier source

    return: La fonction return data, c'est la variable qui contient la séquence Fasta standardisé pour le programme.
    '''
    with open(file) as fasta:
        data = fasta.read()
        data = str.replace(data,"\n","")
    print (data)
    return data

def four_lectures(seq):
    a=''
    seq_inv=seq[::-1]
    seq_comp=''
    for i in seq:
        if i =='A':
            a='T'
        elif i=='T':
            a='A'
        elif i=='G':
            a='C'
        elif i=='C':
            a='G'
        seq_comp += a
    seq_comp_inv=seq_comp[::-1]
    return seq_inv,seq_comp,seq_comp_inv


def writeFasta(data,FICHIER):
    '''Cette fonction n'as pas d'utilité a proprement parlé dans l'éxécution du programme
    elle sert à ecrire dans un fichier txt, les séquences obtenu à la fin de l'éxécution de four-lectures'''
    fic=open(FICHIER, "w")
    for letter in data :
        fic.write(letter)
    fic.close()




if __name__=='__main__':
    print ("Quelle valeur de threshold voulez-vous ? Entrez un nombre et tappez sur ENTREE pour valider votre choix. La valeur par défaut est de 90pb.")
    threshold=int(input())
    data = openFasta("sequence.fasta")
    data_inv,data_comp,data_inv_comp=four_lectures(data)
    writeFasta(data_inv_comp,"fasta_inv_comp.txt")
    writeFasta(data_inv,"fasta_inv.txt")
    writeFasta(data_comp,"fasta_comp.txt")
    isGene3(data,"5'-3'",threshold)
    isGene3(data_inv_comp,"comp_3'-5'",threshold)
