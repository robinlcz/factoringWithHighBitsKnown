import time

def calckn(pourcentage, N):
    R.<X> = PolynomialRing(QQ)
    pr = pourcentage / 2
    polyk = (pr-X)**2 - 4*X*(X+1)*pr
    listroots = polyk.roots(ring=RR)
    (k,_) = listroots[len(listroots) - 1]
    k = ceil(k)
    a = true
    while(a):
        polyn = pr*X**2+(pr-k)*X+k*(k+1)
        listroots = polyn.roots(ring=RR)
        (n0,_) = listroots[0]
        (n1,_) = listroots[1]
        if ceil(n0) <= floor(n1):
            n = ceil(n0)
            while n<=floor(n1):
                if N**(pr*(n+1)/2+k*(k+1)/(2*n))*2**(n/4)*sqrt(n)<N**(k/2):
                    a = False
                    break
                n+=1
        if(a):
            k+=1
    return (k,n-k,n)

def solutions(pourcentage, N, pp, kk):
    Z.<X> = PolynomialRing(ZZ)
    P = X + pp
    L2 = P.coefficients()
    poly = 0
    (k,j,n) = calckn(pourcentage, N)          # On calcule k, j et n
    print("k = ", k, "j = ", j, "n = ", n)
    m = matrix(ZZ,n,n)
    PP=N^k+0*X^1                            # on pose le premier polynome X^k de la matrice m
    b = 10^kk


    for i in range(0,k):
        L=PP.coefficients()                 # on prend les coeff de PP
        #print("Coefficients de PP = ", L)
        for numliste in range(0, len(L)):
            m[i,n-1-numliste] = L[numliste] # on ajoute chaque coefficient a la matrice
        PP=(PP/N)*P                         # a chaque iteration, on divise par N et on multiplie par P (N^k-i P^i)


    L=PP.coefficients()
    for i in range(k, n):
        for numliste in range(0, len(L)):
            m[i, n-1-i+numliste] = L[len(L)-1-numliste] # on rajoute le polynôme X^j*P^k
    for i in range(0,n): 
        m[:,i] *=b^((n-1)-i)

    reduc = copy(m.LLL()[0])                # on execute LLL et on prend le vecteur avec la plus petite norme
    for i in range(0,n):
        reduc[i] /= b**(n-1-i)
    for i in range(0, n):
        poly+=(X^i)*reduc[n-1-i]
    return poly.roots()

def bruteforce(pp, prc, kk, brut):
    if prc<=brut:
        return solutions()
    i=0
    while(prc>brut):
        kk-=1
        prc=kk/nn*1.0
        i+=1
    print("nombre de chiffres à trouver par bruteforce: ", i)
    j=0
    while j<10^i:
        tmp=pp+j*10**(kk)
        sol=solutions(prc, N, tmp, kk)
        if len(sol)==1:
            (test,_)=sol[0]
            print(test)
            tmp+=test
            qq=N/tmp
            if tmp.is_prime() and N%tmp==0 and q.is_prime():
                print (tmp)
                return tmp
        j+=1


def complexityComparison(N, pourcentage,kk):
    j = kk
    i = 0
    l = 0
    pr = pourcentage
    while pr >= 0.5: 
        l += 1
        j -= 1
        pr = j/nn*1.0
    (_,_,n) = calckn(pr,N)
    tempComp = n**5*ln(n)**2;
    while j > 0:
        (_,_,n1) = calckn(j/nn*1.0,N)
        compWithBF = n1**5*ln(n1)**2*10**(kk-j)
        if(compWithBF < tempComp):
            tempComp = compWithBF
            i = j
        j -= 1
    return kk-(i+l)


nn= 100
kk= 49
print("Taille de p :", nn, ", Nombre de chiffres manquants : ", kk)
p=9430189963183715651008289116780446323417798987864230107952609005015226936179180274533260292767432569
q=6618416334439580571371987820009840341586493242076480010115674011556314969757091081187012602370641419
N=p*q
pp=(p//10**kk)*10**kk
print("N =", N) 
print("pp = ", pp)
pourcentage=kk/nn*1.0
print("Il manque donc ", pourcentage*100, "% de p")
tps1 = time.time()
#print("Racine entière de Q obtenue après notre algorithme :", solutions(pourcentage,N,pp,kk))
print("Valeur de p obtenue en utilisant la force brute :", bruteforce(pp, pourcentage, kk,0.47))
tps2 = time.time()
print("Nous avons donc p =", p)
print("Ce qui nous permet d'avoir q =", N/p)
print("Vérification de l'égalité N = pq : ", N == p*q)
print("Temps d'exécution en secondes :", tps2-tps1)
