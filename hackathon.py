filename="MN908947.txt"
f=open(filename,"r")
fp=open("DNA.txt","w")
ft=open("Patient.txt","w")
# Dynamic Programming implementation of LCS problem 
  
def lcs(X , Y): 
    # find the length of the strings 
    m = len(X) 
    n = len(Y) 
  
    # declaring the array for storing the dp values 
    L = [[None]*(n+1) for i in range(m+1)] 
  
   # """Following steps build L[m+1][n+1] in bottom up fashion 
   # Note: L[i][j] contains length of LCS of X[0..i-1] 
   # and Y[0..j-1]"""
    for i in range(m+1): 
        for j in range(n+1): 
            if i == 0 or j == 0 : 
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: 
                L[i][j] = L[i-1][j-1]+1
            else: 
                L[i][j] = max(L[i-1][j] , L[i][j-1]) 
  
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return L[m][n] 
#end of function lcs 
  
X = input("Enter the genetic Sequence of the patient : " )
a=[];b=[]
lines=f.readlines()
for Y in lines:
    k=lcs(X, Y)
    print ("Length of LCS is ",k,"Genetic Subsequence",Y)
    b.append(k)
    a.append(k)
b=sorted(b,reverse=True)

maximum=b[0]
for i in range(len(b)):
    if(maximum==a[i]):
       fp.write(lines[i])

ft.write(X)
fp.close()
f.close()
ft.close()










inputfile ="DNA.txt" 
fp = open(inputfile, "r") 
ft=open("Patient.txt","r")
lines= fp.readlines() 
pop=[];push=[]
  
  
def translate(seq): 
      
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    
    n=len(seq)-1
       
        

    
    if n%3 == 0: 
        for i in range(0, n, 3): 
            codon = seq[i:i + 3] 
            protein= protein+table[codon] 
    return protein 



for line in lines:
    k=translate(line)
    pop.append(k)
readl=ft.read()
t=translate(readl)

def lcs(X , Y): 
    # find the length of the strings 
    m = len(X) 
    n = len(Y) 
  
    # declaring the array for storing the dp values 
    L = [[None]*(n+1) for i in range(m+1)] 
  
   # """Following steps build L[m+1][n+1] in bottom up fashion 
   # Note: L[i][j] contains length of LCS of X[0..i-1] 
   # and Y[0..j-1]"""
    for i in range(m+1): 
        for j in range(n+1): 
            if i == 0 or j == 0 : 
                L[i][j] = 0
            elif X[i-1] == Y[j-1]: 
                L[i][j] = L[i-1][j-1]+1
            else: 
                L[i][j] = max(L[i-1][j] , L[i][j-1]) 
  
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1] 
    return L[m][n] 
#end of function lcs 
  
percentage=[]
for Y in pop:
    k=lcs(t, Y)
    percentage.append((k/len(t))*100)
for i in range(len(percentage)):
    print("Genetic Sequence",lines[i],"percentage",percentage[i])
ft.close()
fp.close()    
