max_range_pow = 10


def generate_functions():
    filepath = 'functions.py'
    with open(filepath, 'w') as outfile:
        for i in range(max_range_pow):
            outfile.write('' +
                          f'def f_neg_pow{i}(dist): '
                          f'return (1.0 / dist) ** {i}\n'
                          )
        outfile.write('functions = ')
        outfile.write('[')
        for i in range(max_range_pow - 1):
            outfile.write(f'f_neg_pow{i}, ')
        outfile.write(f'f_neg_pow{max_range_pow - 1}')
        outfile.write(']')


if __name__ == '__main__':
    generate_functions()
