max_power = 100

def mean(a:float, d:float, p:float):
    '''Computes the value of m_p(a,d)'''
    if p==float("inf"):
        return max(abs(a), abs(d))
    elif p==-float("inf"):
        return min(abs(a), abs(d))
    elif p==0:
        return a**(1/2) * d**(1/2)
    else:
        return (1/2*(a**p + d**p))**(1/p)
    
    
def additive_inverse(a:float, b:float, c:float, d:float):
    '''(a,b,c,d) -> (-a, -b, -c, -d)'''
    return -a, -b, -c, -d

def exchange_of_the_means(a:float, b:float, c:float, d:float):
    '''(a,b,c,d) -> (a,c,b,d)'''
    return a, c, b, d

def exchange_of_the_extremes(a:float, b:float, c:float, d:float):
    '''(a,b,c,d) -> (d,b,c,a)'''
    return d, b, c, a

def inversion_of_ratios(a:float, b:float, c:float, d:float):
    '''(a,b,c,d) -> (b,a,d,c)'''
    return b, a, d, c

def canonical_form(a:float, b:float, c:float, d:float):

    if sum(x<0 for x in [a, b, c, d])>=3:
        a, b, c, d = additive_inverse(a,b,c,d)
    if a > d:
        a, b, c, d = exchange_of_the_extremes(a,b,c,d)
    if b > c:
        a, b, c, d = exchange_of_the_means(a,b,c,d)
    if a > b:
        a, b, c, d = inversion_of_ratios(a, b, c, d)
    return a, b, c, d



def _DS_cicle(S:str, lower:float, upper:float, f):
    assert(f(lower)*f(upper)<=0)

    eps = 10**-6 if S == 'R' else 2
    
    while abs(upper-lower) > eps:
        middle = (lower+upper)/2

        if (S=='PO' or S=='O') and middle%2==0:
            middle += 1
        elif S == 'NO' and middle%2==0:
            middle -= 1
        elif S == 'E' and middle%2==1:
            middle = middle + 1 if middle!=-1 else middle-1
        elif S == 'E' and middle==0:
            middle = 2
            if lower==-2 and upper==2:
                eps = 4
        
        if f(middle)>0:
            upper = middle
        else:
            lower = middle

    return (lower, upper) 

def DS(a:float, b:float, c:float, d:float, S:str, lower:float, upper:float):
    f = lambda p : mean(a,d,p)-mean(b,c,p)
    lower, upper = _DS_cicle(S,lower,upper,f)
    if f(lower) == 0:
        return [lower]
    elif f(upper) == 0:
        return [upper]
    elif S == 'R':
        return [(lower+upper)/2]
    else:
        return []



def case_i_set_choice(a:float, b:float, c:float, d:float):
    '''Gives numbers l and u with m_l(a,d) < m_l(b,c) and m_u(a,d) > m_u(b,c)'''
    success_flag = True

    lower = -2
    while success_flag and mean(a,d,lower) >= mean(b,c,lower):
        lower *= 2
        if abs(lower)>=max_power:
            success_flag = False

    upper = 2
    while success_flag and mean(a,d,upper) <= mean(b,c,upper):
        upper *= 2
        if upper>=max_power:
            success_flag = False

    return (lower, upper, success_flag)

def case_i(a:float, b:float, c:float, d:float):
    '''0 < a < b <= c, d'''

    assert(0<a and 0<b and 0<c and 0<d)

    sol = []

    if a<b<=c<d and mean(a,d,0) != mean(b,c,0):
        #Search for the solution
        lower, upper, success = case_i_set_choice(a, b, c, d)
        if success:
            sol += DS(a,b,c,d,'R',lower,upper)
        
    return sol

def case_ii_b(a:float, b:float, c:float, d:float):
    '''a < 0 < b <= d <= c'''

    assert(a<0<b<=d<=c)

    sol = []

    if b<abs(a)<c:
        a_even, b_even, c_even, d_even = canonical_form(abs(a),b,c,d)
        lower_even, upper_even, success_even = case_i_set_choice(a_even,b_even,c_even,d_even)
        if success_even:
            sol += DS(a_even,b_even,c_even,d_even,'E',lower_even,upper_even)
    
    return sol

def case_ii_c_set_choice(a:float, b:float, c:float, d:float):
    '''Gives number u such that m_u(a,d)>m_u(b,c)'''
    success_flag = True

    lower = 1

    upper = 1
    while success_flag and mean(a,d,upper) <= mean(b,c,upper):
        upper *= 3
        if upper > max_power:
            success_flag = False
    return lower, upper, success_flag
    
def case_ii_c(a:float, b:float, c:float, d:float):
    '''a < 0 < b <= c <= d'''

    assert(a<0<b<=c<=d)

    sol = []

    if abs(a)<b:
        a_even, b_even, c_even, d_even = canonical_form(abs(a),b,c,d)
        lower_even, upper_even, success_even = case_i_set_choice(a_even,b_even,c_even,d_even)

        if success_even:
            sol += DS(a_even,b_even,c_even,d_even,'E',lower_even,upper_even)
    
    if abs(a)<d and mean(a,d,1) <= mean(b,c,1):
        lower_odd, upper_odd, success_odd = case_ii_c_set_choice(a,b,c,d)

        if success_odd: 
            sol += DS(a,b,c,d,'PO',lower_odd,upper_odd)
    
    return sol


def case_ii(a:float, b:float, c:float, d:float):
    '''a < 0 < b, c, d'''

    assert(a<0<b<=c and d>0)
    
    if d < b:
        return [-p for p in case_ii_c(1/a, 1/c, 1/b, 1/d)]
            
    elif b<=d<=c:
        return case_ii_b(a, b, c, d)
        
    else:
        return case_ii_c(a, b, c, d)
    

def case_iv_set_choice_odd(a:float, b:float, c:float, d:float):
    '''Gives numbers l and u such that m_l(a,d) < m_l(b,c) and m_u(a,d) > m_u(b,c)'''
    success_flag = True

    lower = -1
    while success_flag and mean(a,d,lower) >= mean(b,c,lower):
        lower *= 3
        if abs(lower) >= max_power:
            success_flag = False

    upper = 1
    while success_flag and mean(a,d,upper) <= mean(b,c,upper):
        upper *= 3
        if upper >= max_power:
            success_flag = False

    return lower, upper, success_flag

def case_iv(a:float, b:float, c:float, d:float):
    '''a <= b < 0 < c, d'''

    assert(a<=b<0 and c>0 and d>0)

    sol = []

    a_odd, b_odd, c_odd, d_odd = canonical_form(abs(a),abs(b),d,c)
    a_even, b_even, c_even, d_even = canonical_form(abs(a),abs(b),c,d)
    
    if a_even<b_even<=c_even<d_even:
        lower_even, upper_even, success_even = case_i_set_choice(a_even,b_even,c_even,d_even)

        if success_even:
            sol += DS(a_even,b_even,c_even,d_even,'E',lower_even,upper_even)
    
    if c_odd < d_odd:
        lower_odd, upper_odd, success_odd = case_iv_set_choice_odd(a_odd, b_odd, c_odd, d_odd)

        if success_odd:
            sol += DS(a_odd,b_odd,c_odd,d_odd,'O',lower_odd,upper_odd)
    
    return sol


def case_x(a:float, b:float, c:float, d:float):
    '''a <= d < 0 < b <= c'''

    assert(a<=d<0<b<=c)

    sol = []

    if b<abs(d)<=abs(a)<c or abs(d)<b<=c<abs(a):
        a, b, c, d = canonical_form(abs(a),b,c,abs(d))
        lower_even, upper_even, success_even = case_i_set_choice(a,b,c,d)

        if success_even:
            sol += DS(a,b,c,d,'E',lower_even,upper_even)
    
    return sol


def alg(a:float, b:float, c:float, d:float) -> list:
    assert(a!=0 and b!=0 and c!=0 and d!=0)

    a, b, c, d = canonical_form(a,b,c,d)
    m = max(abs(a), abs(b), abs(c), abs(d))

    a, b, c, d = a/m, b/m, c/m, d/m
    
    aux = []
    for p in [-float('inf'), 0, float('inf')]:
        if mean(a,d,p) == mean(b,c,p):
            aux += [p]
    
    if a==b and c==d:
        if a<0 and abs(a)==c:
            return aux + ['IR - (Z^*_- - 2Z)']
        return aux + ['IR']

    a_abs, b_abs, c_abs, d_abs = canonical_form(abs(a), abs(b), abs(c), abs(d))
    if a_abs==b_abs and c_abs==d_abs:
        return aux + ['2Z^*']

    if a<0 and b<0 and abs(a)==d and abs(b)==c:
        return aux + ['Z^*_+ - 2Z']

    elif a>0 and b>0 and c>0 and d>0:
        return aux + case_i(a, b, c, d)
    
    elif a<0 and b>0 and c>0 and d>0:
        return aux + case_ii(a, b, c, d)
    
    elif a<0 and b<0 and c>0 and d>0:
        return aux + case_iv(a, b, c, d)
    
    elif a<0 and b>0 and c>0 and d<0:
        return aux + case_x(a, b, c, d)
    
    else:
        return
    

if __name__ == "__main__":
    a = 1
    b = 2
    c = 3
    d = 4.1
    result = alg(a, b, c, d)
    print(f"{a}:{b}::^p{c}:{d} for p in {set(result)} ")