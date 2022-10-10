#cashier.py
prices = {'milk': 1.00, 'wine': 2.50, 'apples': 0.6}
 
def sum_bill(purchase):
    """Calculate the total amount to pay"""
    total = 0
    for item, quantity in purchase:
        total += prices[item]*quantity
    return total
 
#Testing the code
if __name__=='__main__':
    my_purchase = [('milk', 2), ('wine', 1),  
                   ('apples', 1.2)]
    bill = sum_bill(my_purchase)
 
    print 'I owe %.2f Euros' % bill
