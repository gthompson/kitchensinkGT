# count events of each mainclass / subclass
import pandas as pd
import os
SEISAN_DATA = os.path.join( os.getenv('HOME'),'DATA','MVO')
SEISAN_DB = 'MVOE_'

allcsv = os.path.join(SEISAN_DATA, 'reawav_%sall.csv' % SEISAN_DB)
dfall = pd.read_csv(allcsv)
print(len(dfall))

from collections import Counter
counter_of_classes = Counter(dfall['mainclass'])
print(counter_of_classes.most_common())

counter_of_subclasses = Counter(dfall['subclass'])
print(counter_of_subclasses.most_common())
