KARATSUBA_CUTOFF = 2**10  # note) this is an arbitrary value
from typing import List

def karatsuba_mul(a:int,b:int):
    # a = a_front*2**n_half + a_back
    # b = b_front*2**n_half + b_back
    # a*b = z1**2*n + z2*2**n_half + z3
    a_len = a.bit_length()
    b_len = b.bit_length()
    n = max(a_len,b_len)
    if n<=KARATSUBA_CUTOFF:
        return a*b
    n_half = n>>1
    a_front = a>>n_half
    b_front = b>>n_half
    a_back  = a&((1<<n_half)-1)
    b_back  = b&((1<<n_half)-1)
    z1 = karatsuba_mul(a_front,b_front)
    z3 = karatsuba_mul(a_back,b_back)
    z2 = karatsuba_mul(a_back+a_front,a_back+b_back) - z1 - z3
    return (z1<<n_half<<n_half) + (z2<<n_half) + z3

def get_size_power_of_2(a:int,b:int):
    # to compute FFT(cooley-tukey), size = 2**n
    size = 1
    a_len = a.bit_length()
    b_len = b.bit_length()
    thres = max(a_len, b_len)
    cnt=0
    while size<thres:
        size = size<<1
        cnt+=1
    return [size<<1,cnt+1]

def num_to_array(num:int,l:int,b:int):
    # slice num
    N=[]
    cut = (1<<l)-1
    for _ in range(b):
        N.append(num&cut)
        num = num >> l
    return N
    
def bit_reversal(x:List):
    rev = 0
    for i in range(len(x)):
        bit = i
        rev = 0
        for _ in range(len(x).bit_length()-1):
            rev = (rev<<1) + (i&1)
            i = i >> 1
        if bit<rev:
            tmp = x[bit]
            x[bit] = x[rev]
            x[rev] = tmp    
    return x

def mod_bit_op(a:int,l:int):
    # mod_bit_op(a,l) == a%(2^2l+1)
    # a1(2^(2l)+1)+a2 = a
    a1 = a>>(2*l)
    if a1 > 1<<(l<<1):
        a1 = mod_bit_op(a1,l)
    a_back = (a&((1<<(l<<1))-1))
    result = a_back - a1
    if result<0:
        return result+((1<<(l<<1))+1)
    return result

def fft_in_place(U:List,omegas:List,l:int,b:int):
    U = U[:]
    U = bit_reversal(U)
    iterations = len(U).bit_length()-1
    M = 2
    m = 1 # m : M = 2^m
    for _ in range(iterations):
        for i in range(0, len(U), M):
            g = 0
            for j in range(0, M >> 1):
                k = i + j + (M >> 1)
                even = U[i + j]
                odd = U[k]*omegas[g]
                # note) U[k]*omegas[g] should be calculated in SSA
                U[i + j] = mod_bit_op(even + odd,l)
                U[k] = mod_bit_op(even-odd,l)
                g = g + (len(U) >> m)
                g = g&(b-1)
        m+=1
        M <<= 1
    return U


def pointwise_multiplication(U:List,B:List,l:int):
    C = [0]*len(U)
    for i in range(len(C)):
        C[i] = mod_bit_op(U[i]*B[i],l)
        # note) U[i]*B[i] should be calculated in SSA
    return C

def array_to_num(Y:List, l:int):
    # carry out by l
    y = 0
    for u in reversed(Y):
        y = u + (y<<l)
    return y

def mod_by_b(U:List[int],W:List[int],b:int):
    logb = b.bit_length()-1
    U_with_pad = 0
    W_with_pad = 0
    for i in range(len(U)):
        U_with_pad = (U_with_pad<<(3*logb))|U[-1-i]&(b-1)
        W_with_pad = (W_with_pad<<(3*logb))|W[-1-i]&(b-1)
    Y_with_pad = karatsuba_mul(U_with_pad, W_with_pad)
    Y_mod_b = [0]*len(U)
    for i in range(len(Y_mod_b)):
        Y_mod_b[i] = Y_with_pad&(b-1)
        Y_with_pad = Y_with_pad>>(3*logb)
    return Y_mod_b

def chinese_remainder(Y_idft:List, Y_mod_b:List, b:int, l:int):
    Y = [0]*len(Y_idft)

    for i in range(len(Y)):        
        # Y[i] = (((Y_mod_b[i] - Y_idft[i])<<(l<<1)) + Y_mod_b[i])%(b*(2**(2*l)+1))
        delta = (Y_mod_b[i] - Y_idft[i])&(b-1)
        Y[i] = Y_idft[i]+delta+(delta<<(l<<1))
        # TODO
        # why it works? : https://math.stackexchange.com/questions/72771/sch%c3%b6nhage-strassen-multiplication
    return Y

def omega_to_list(diff_ceil, b, l):
    # omega = 1<<(4<<diff_ceil)
    omegas = [1]*b
    inverse_omegas = [1]*b
    log2_omega=0
    for i in range(1,b):
        omegas[i] = omegas[i-1]<<(4<<diff_ceil)
        log2_omega+=(4<<diff_ceil)
        # note) omegas[i]*inverse_omegas[i]=2**(4l)
        inverse_omegas[i] = mod_bit_op(1<<((l<<2)-log2_omega),l)
    return [omegas, inverse_omegas]

def multiplication(n1:int,n2:int):
    n,k = get_size_power_of_2(n1,n2) # n=2**k
    diff_ceil = k&1
    l = 1<<((k>>1)+diff_ceil) # l=2**k/2 ceiling
    b = 1<<(k>>1) # b=2**k/2 floor

    U = num_to_array(n1,l,b)
    W = num_to_array(n2,l,b)
    # print(list(map(lambda x:bin(x)[2:],reversed(U))))
    # print(list(map(lambda x:bin(x)[2:],reversed(W))))

    omegas, inverse_omegas = omega_to_list(diff_ceil,b,l)
    U_dft = fft_in_place(U,omegas,l,b)
    W_dft = fft_in_place(W,omegas,l,b)

    Y_dft = pointwise_multiplication(U_dft,W_dft,l)

    Y_idft = fft_in_place(Y_dft, inverse_omegas,l,b)
    Y_idft = list(map(lambda x:mod_bit_op(x*(1<<(l<<2)-(b.bit_length()-1)),l),Y_idft))
    Y_mod_b = mod_by_b(U,W,b)

    Y = chinese_remainder(Y_idft,Y_mod_b, b, l)
    y = array_to_num(Y, l)
    return y
11
n1 = 101010101010101152415241264561425162461241624162411111111111111
n2 = 10101152152523352523
# print(n1*n2)
# print(karatsuba_mul(n1,n2))
# n1 = 0b1
# n2 = 0b1
print(multiplication(n1,n2))
print(n1*n2)