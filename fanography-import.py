import fanography
import pprint

fanos = fanography.application.fanos

hodge = dict()

for rho in fanos:
    for ID in fanos[rho]:
        X = fanos[rho][ID]
        hodge[(rho, ID)] = X.h12

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(hodge)
