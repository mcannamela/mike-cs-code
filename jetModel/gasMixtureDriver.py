
#ensure the directory is on the path and import the module
import sys
modelDir = "C:\\Documents and Settings\\m.cannamela\\My Documents\\svnCheckout\\jetModel"
comparator = lambda s: s==modelDir
b = not any([comparator(s) for s in sys.path])
if b:
    sys.path+=[modelDir]


import gasMixture
reload(gasMixture)
a = gasMixture.gasMixture()
a.blend(['argon','helium'],[.68,.32])

