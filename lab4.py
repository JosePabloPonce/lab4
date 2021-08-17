import math
import scipy.special as ss
import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import gf2matrix
import re
import random
from collections import Counter
import pandas

def test01(input, n):
        
    ones = input.count('1') #number of ones

    zeroes = input.count('0')    #number of zeros

    s = abs(ones - zeroes)  

    p = math.erfc(float(s)/(math.sqrt(float(n)) * math.sqrt(2.0))) #p-value

    success = ( p >= 0.01)  # success = true if p-value >= 0.01

    return [p, success]

def test02(input, n, M=32):
    # Compute number of blocks M = block size. N=num of blocks
    # N = floor(n/M)
    # miniumum block size 20 bits, most blocks 100
    # fieldnames = ['number','chisq','p-value', 'success']
    
    N = int(math.floor(n/M))
    
    if N > 99:
        N=99
        M = int(math.floor(n/N))

    if n < 100:
        # Too little data for test. Input of length at least 100 bits required
        return [0.0, 0.0, False]

    num_of_blocks = N

    block_size = M 

    proportions = list()

    for i in range(num_of_blocks):
        
        block = input[i*(block_size):((i+1)*(block_size))]
        
        ones = block.count('1')

        zeroes = block.count('0') 
        
        proportions.append(Fraction(ones,block_size))

    chisq = 0.0

    for prop in proportions:
        chisq += 4.0*block_size*((prop - Fraction(1,2))**2)
    
    p = ss.gammaincc((num_of_blocks/2.0),float(chisq)/2.0) # p-value
    
    success = (p>= 0.01)

    return [p, success]

def test03(input, n):
            
    ones = input.count('1') #number of ones

    zeroes = input.count('0')    #number of zeros

    prop = float(ones)/float(n)

    tau = 2.0/math.sqrt(n)

    vobs = 0.0

    if abs(prop-0.5) > tau:
        p = 0
    else:

        vobs = 1.0
        for i in range(n-1):
            if input[i] != input[i+1]:
                vobs += 1.0

        p = math.erfc(abs(vobs - (2.0*n*prop*(1.0-prop)))/(2.0*math.sqrt(2.0*n)*prop*(1-prop) ))
    
    success = (p>=0.01)


    return [p, success]


def test04(input, n):

    M8 = [0.2148, 0.3672, 0.2305, 0.1875]

    # Length of blocks
    M = 8 
                
    K = 3

    N = 16
            
    # Table of frequencies
    v = [0,0,0,0,0,0,0]

    for i in range(N): # over each block
        #find the longest run
        block = input[i*M:((i+1)*M)] # Block i
        
        run = 0
        longest = 0
        for j in range(M): # Count the bits.
            if block[j] == '1':
                run += 1
                if run > longest:
                    longest = run
            else:
                run = 0

        if longest <= 1:    v[0] += 1
        elif longest == 2:  v[1] += 1
        elif longest == 3:  v[2] += 1
        else:               v[3] += 1
    
    # Compute Chi-Sq
    chi_sq = 0.0
    for i in range(K+1):
        p_i = M8[i]
        upper = (v[i] - N*p_i)**2
        lower = N*p_i
        chi_sq += upper/lower
    # p-value
    p = ss.gammaincc(K/2.0, chi_sq/2.0)
    
    success = (p>=0.01)

    return [p, success]

def test05(input, n, M=32, Q=32):
    N = int(math.floor(n/(M*Q))) #Number of blocks
    
    if N < 38:
        print("  Number of blocks must be greater than 37")
        p = 0.0
        return [0]*9
        
    # Compute the reference probabilities for FM, FMM and remainder 
    r = M
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FR_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    r = M-1
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FRM1_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    LR_prob = 1.0 - (FR_prob + FRM1_prob)
    
    FM = 0      # Number of full rank matrices
    FMM = 0     # Number of rank -1 matrices
    remainder = 0
    for blknum in range(N):
        block = [None] * (M*Q)
        
        for i in range(M*Q):
            block[i] = int(input[blknum*M*Q + i],2)
            
        # Put in a matrix
        matrix = gf2matrix.matrix_from_bits(M,Q,block,blknum) 
        rank = gf2matrix.rank(M,Q,matrix,blknum)


        if rank == M: # count the result
            FM += 1
        elif rank == M-1:
            FMM += 1  
        else:
            remainder += 1

    chisq =  (((FM-(FR_prob*N))**2)/(FR_prob*N))
    chisq += (((FMM-(FRM1_prob*N))**2)/(FRM1_prob*N))
    chisq += (((remainder-(LR_prob*N))**2)/(LR_prob*N))
    
    p = math.e **(-chisq/2.0)

    success = (p >= 0.01)

    return [p, success]

def padding(input, n):
	while len(input) <n:
		input = '0' + input
	return input

def berlekamp_massey(input):
    n = len(input)

    b = '1' + '0'*(n-1)
    c = '1' + '0'*(n-1)

    L = 0
    m = -1
    N = 0
    while (N < n):
        #compute discrepancy
        d = int(input[N],2)
        if L>0:
            k = c[1:L+1]
            h = input[N-L:N][::-1]

            k = int(k,2)  #binary to integer

            h = int(h,2)    #binary to integer

            r = k&h #bitwise-and

            r = bin(r)[2:] 

            r = r.count('1')

            d = d ^ (r%2)

        if d != 0:  # If d is not zero, adjust poly
            t = c
            k = c[N-m:n]
            k = int(k, 2)
            h = b[0:n-N+m]
            h = int(h,2)
            k = k^h
            c = c[0:N-m] + padding(bin(k)[2:], n-N+m)
            # print(c)
            if (L <= (N/2)):
                L = N + 1 - L
                m = N
                b = t 
        N = N +1
    # Return length of generator and the polynomial
    return L , c[0:L]
    
def test10(input, n, patternlen=None):
    n = len(input)
    # Step 1. Choose the block size
    if patternlen != None:
        M = patternlen  
    else: 
        if n < 1000000:
            print("Error. Need at least 10^6 input")
            return [0]*9
        M = 512
    K = 6 
    N = int(math.floor(n/M))

    # Step 2 Compute the linear complexity of the blocks
    LC = list()
    for i in range(N):
        x = input[(i*M):((i+1)*M)]
        t = berlekamp_massey(x)[0]
        LC.append(t)

    # Step 3 Compute mean
    a = float(M)/2.0
    b = ((((-1)**(M+1))+9.0))/36.0
    c = ((M/3.0) + (2.0/9.0))/(2**M)
    mu =  a+b-c
    
    T = list()
    for i in range(N):
        x = ((-1.0)**M) * (LC[i] - mu) + (2.0/9.0)
        T.append(x)
        
    # Step 4 Count the distribution over Ticket
    v = [0,0,0,0,0,0,0]
    for t in T:
        if t <= -2.5:
            v[0] += 1
        elif t <= -1.5:
            v[1] += 1
        elif t <= -0.5:
            v[2] += 1
        elif t <= 0.5:
            v[3] += 1
        elif t <= 1.5:
            v[4] += 1
        elif t <= 2.5:
            v[5] += 1            
        else:
            v[6] += 1

    # Step 5 Compute Chi Square Statistic
    pi = [0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]
    chi_sq = 0.0
    for i in range(K+1):
        chi_sq += ((v[i] - (N*pi[i]))**2.0)/(N*pi[i])

    # Step 6 Compute P Value
    P = ss.gammaincc((K/2.0),(chi_sq/2.0))

    success = (P >= 0.01)
    
    return [P, success]


def test15(input, n):

    x = list()             # Convert to +1,-2
    for i in range(n):
        x.append(int(input[i],2)*2-1)

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)  
    # print(s)  
    sprime = [0]+s+[0] # Add 0 on each end

    # Count the number of cycles J
    J = 0
    for value in sprime[1:]:
        if value == 0:
            J += 1
            
    # Build the counts of offsets
    count = [0 for x in range(-9,10)]
    for value in sprime:
        if (abs(value) < 10):
            count[value] += 1

    # Compute P values
    success = True
    plist = list() # list of p-values for each state
    p_average = 0.0
    for x in range(-9,10): 
        if x != 0:
            top = abs(count[x]-J)
            bottom = math.sqrt(2.0 * J *((4.0*abs(x))-2.0))
            p = ss.erfc(top/bottom)

            # print("p[" + str(x) +"] = " + str(p))

            p_average +=p
            plist.append(p)
            if p < 0.01:
                success = False

    return [p, success]


def test09(input, n, patternlen=None, initblocks=None):

    # Step 1. Choose the block size
    if patternlen != None:
        L = patternlen  
    else: 
        ns = [904960,2068480,4654080,10342400,
              22753280,49643520,107560960,
              231669760,496435200,1059061760]
        L = 6
        if n < 387840:
            # Too little data. Inputs of length at least 387840 are recommended
            return [0] * 8
        for threshold in ns:
            if n >= threshold:
                L += 1 

    # Step 2 Split the data into Q and K blocks
    nblocks = int(math.floor(n/L))
    if initblocks != None:
        Q = initblocks
    else:
        Q = 10*(2**L)
    K = nblocks - Q
    
    # Step 3 Construct Table
    nsymbols = (2**L)
    T=[0 for x in range(nsymbols)] # zero out the table
    for i in range(Q):             # Mark final position of
        pattern = input[i*L:(i+1)*L] # each pattern
        idx = int(pattern, 2)
        T[idx]=i+1      # +1 to number indexes 1..(2**L)+1
                        # instead of 0..2**L
    # Step 4 Iterate
    sum = 0.0
    for i in range(Q,nblocks):
        pattern = input[i*L:(i+1)*L]
        j = int(pattern,2)
        dist = i+1-T[j]
        T[j] = i+1
        sum = sum + math.log(dist,2)
    
    # Step 5 Compute the test statistic
    fn = sum/K
       
    # Step 6 Compute the P Value
    # Tables from https://static.aminer.org/pdf/PDF/000/120/333/
    # a_universal_statistical_test_for_random_bit_generators.pdf
    ev_table =  [0,0.73264948,1.5374383,2.40160681,3.31122472,
                 4.25342659,5.2177052,6.1962507,7.1836656,
                 8.1764248,9.1723243,10.170032,11.168765,
                 12.168070,13.167693,14.167488,15.167379]
    var_table = [0,0.690,1.338,1.901,2.358,2.705,2.954,3.125,
                 3.238,3.311,3.356,3.384,3.401,3.410,3.416,
                 3.419,3.421]
                 
    # sigma = math.sqrt(var_table[L])
    sigma = abs((fn - ev_table[L])/((math.sqrt(var_table[L]))*math.sqrt(2)))
    P = math.erfc(sigma)

    success = (P >= 0.01)
    return [P, success]


def test12(input, n):
    
    m = int(math.floor(math.log(n,2)))-6
    if m < 2:
        m = 2
    if m >3 :
        m = 3
    
    Cmi = list()
    phi_m = list()
    for iterm in range(m,m+2):
        # Step 1 
        padded_input=input+input[0:iterm-1]
    
        # Step 2
        counts = list()
        for i in range(2**iterm):
            count = 0
            for j in range(n):
                if int(padded_input[j:j+iterm],2) == i:
                    count += 1
            counts.append(count)
    
        # step 3
        Ci = list()
        for i in range(2**iterm):
            Ci.append(float(counts[i])/float(n))
        
        Cmi.append(Ci)
    
        # Step 4
        sum = 0.0
        for i in range(2**iterm):
            if (Ci[i] > 0.0):
                sum += Ci[i]*math.log((Ci[i]/10.0))
        phi_m.append(sum)
        
    # Step 5 - let the loop steps 1-4 complete
    
    # Step 6
    appen_m = phi_m[0] - phi_m[1]

    chisq = 2*n*(math.log(2) - appen_m)

    # Step 7
    p = ss.gammaincc(2**(m-1),(chisq/2.0))
    
    success = (p >= 0.01)

    return [p, success]

# RANDOM EXCURSION TEST
def test14(input, n):

    # Convert to +1,-1
    x = list()
    for i in range(n):
        x.append(int(input[i],2)*2 -1 )

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)    
    sprime = [0]+s+[0] # Add 0 on each end
    

    # Build the list of cycles
    pos = 1
    cycles = list()
    while (pos < len(sprime)):
        cycle = list()
        cycle.append(0)
        while sprime[pos]!=0:
            cycle.append(sprime[pos])
            pos += 1
        cycle.append(0)
        cycles.append(cycle)
        pos = pos + 1
    
    J = len(cycles)  
    
    vxk = [['a','b','c','d','e','f'] for y in [-4,-3,-2,-1,1,2,3,4] ]

    # Count Occurances  
    for k in range(6):
        for index in range(8):
            mapping = [-4,-3,-2,-1,1,2,3,4]
            x = mapping[index]
            cyclecount = 0
            #count how many cycles in which x occurs k times
            for cycle in cycles:
                oc = 0
                #Count how many times x occurs in the current cycle
                for pos in cycle:
                    if (pos == x):
                        oc += 1
                # If x occurs k times, increment the cycle count
                if (k < 5):
                    if oc == k:
                        cyclecount += 1
                else:
                    if k == 5:
                        if oc >=5:
                            cyclecount += 1
            vxk[index][k] = cyclecount
    
    # Table for reference random probabilities 
    pikx=[[0.5     ,0.25   ,0.125  ,0.0625  ,0.0312 ,0.0312],
          [0.75    ,0.0625 ,0.0469 ,0.0352  ,0.0264 ,0.0791],
          [0.8333  ,0.0278 ,0.0231 ,0.0193  ,0.0161 ,0.0804],
          [0.875   ,0.0156 ,0.0137 ,0.012   ,0.0105 ,0.0733],
          [0.9     ,0.01   ,0.009  ,0.0081  ,0.0073 ,0.0656],
          [0.9167  ,0.0069 ,0.0064 ,0.0058  ,0.0053 ,0.0588],
          [0.9286  ,0.0051 ,0.0047 ,0.0044  ,0.0041 ,0.0531]]
    
    success = True
    plist = list()
    chi_sq = list()
    p_total = 0.0
    for index in range(8):
        #list of states
        mapping = [-4,-3,-2,-1,1,2,3,4] 
        x = mapping[index]
        chisq = 0.0
        for k in range(6):
            top = float(vxk[index][k]) - (float(J) * (pikx[abs(x)-1][k]))
            top = top*top
            bottom = J * pikx[abs(x)-1][k]
            chisq += top/bottom

        p = ss.gammaincc(5.0/2.0,chisq/2.0)
        p_total += p
        plist.append(p)
        chi_sq.append(chisq)
        if p < 0.01:
            success = False

    return [p, success]



def resumenestadistico(cadenadebits):
    print(test01(cadenadebits, len(cadenadebits)))
    print(test02(cadenadebits, len(cadenadebits)))
    print(test03(cadenadebits, len(cadenadebits)))
    print(test04(cadenadebits, len(cadenadebits)))
    print(test05(cadenadebits, len(cadenadebits)))
    print(test10(cadenadebits, len(cadenadebits)))
    print(test15(cadenadebits, len(cadenadebits)))
    print(test09(cadenadebits, len(cadenadebits)))
    print(test12(cadenadebits, len(cadenadebits)))
    print(test14(cadenadebits, len(cadenadebits)))

    

def resumenestadisticohistograma():
    lista =[]
    for i in range (0, 1000):
        cadenadebits = LCG()
        a= (test01(cadenadebits, len(cadenadebits)))
        b= (test02(cadenadebits, len(cadenadebits)))
        c= (test03(cadenadebits, len(cadenadebits)))
        d= (test04(cadenadebits, len(cadenadebits)))
        e= (test05(cadenadebits, len(cadenadebits)))
        f= (test10(cadenadebits, len(cadenadebits)))
        g= (test15(cadenadebits, len(cadenadebits)))
        h= (test09(cadenadebits, len(cadenadebits)))
        i= (test12(cadenadebits, len(cadenadebits)))
        j= (test14(cadenadebits, len(cadenadebits)))
        print(a)
        if(a[1] == False):
            lista.append('1')
            
        if(b[1] == False):
            lista.append('2')

        if(c[1] == False):
            lista.append('3')

        if(d[1] == False):
            lista.append('4')

        if(e[1] == False):
            lista.append('5')

        if(f[1] == False):
            lista.append('6')

        if(g[1] == False):
            lista.append('7')

        if(h[1] == False):
            lista.append('8')

        if(i[1] == False):
            lista.append('9')

        if(j[1] == False):
            lista.append('10')

    letter_counts = Counter(lista)
    df = pandas.DataFrame.from_dict(letter_counts, orient='index')
    df.plot(kind='bar')
    plt.show()


kbits = 1
def LCG(m=2**31-1, a=1103515245, b=12345):
    x = np.random.choice(m)
    bits=""
    for i in range(0, 10**6):
        x = (a*x + b) % m
        stb = ''.join(format(ord(i), '08b') for i in str(x))
        stb = stb[::-1]
        for i in range(0, kbits):
            bits+=stb[i]
    return bits

s1 = LCG()

resumenestadistico(s1)
resumenestadisticohistograma()