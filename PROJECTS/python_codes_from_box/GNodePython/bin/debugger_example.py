
def add(arg1, arg2):
    return arg1 + arg2

def sub(arg1, arg2):
    return arg1 - arg2

def div(arg1, arg2):
    return arg1 / arg2

def sum_over_difference(arg1, arg2):
    """Compute sum of arguments over difference of arguments."""
    arg_sum = add(arg1, arg2)
    arg_diff = sub(arg1, arg2)
    frac_sum_diff = div(arg_sum, arg_diff)
    return frac_sum_diff

if __name__ == '__main__':
    result = sum_over_difference(12., 3.)
    print result
