from random import uniform

def generate_quadruple(n:int):
    '''generates four non-null numbers between -10**n and 10**n'''
    return [uniform(-10**n, 0) if uniform(0,1)<1/2 else -uniform(-10**n, 0) for _ in range(4)]

def gen_dataset(n:int):
    f = open(f'tests/test-{n}.txt', "a")
    for _ in range(10**6):
        a, b, c, d = generate_quadruple(n)
        f.write(f'{a} {b} {c} {d}\n')

    f.close()


if __name__ == "__main__":
    n_max = 10
    for n in range(n_max):
        gen_dataset(n)