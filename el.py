els = {}
with open('elements', 'r') as f:
    for l in f.readlines():
        num, name = l.split()
        num = int(num)

        els[num] = name


def el_name(num_code):
    if num_code > 9999:
        print('Only works for isotopes with mass <= 99')
        return

    mass = num_code % 100
    Z = num_code // 100

    return f'{mass}{els[Z]}'


if __name__ == '__main__':
    from sys import argv

    if len(argv) < 2:
        print('Error: no input')
        exit(0)
    else:
        try:
            code = int(argv[1])
        except:
            print('Error: code needs to be a number')
        else:
            print(el_name(code))
