import random

target = 10

if __name__ == '__main__':
    ls = [(19.5 + random.random()) / 20 * target for n in range(3)]
    print('\t'.join(map(str, ls)))
