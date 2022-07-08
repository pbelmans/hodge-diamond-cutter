import pprint
import fanography

fanos = fanography.application.fanos

hodge = {}

for rho in fanos:
    for ID in fanos[rho]:
        X = fanos[rho][ID]
        hodge[(rho, ID)] = X.h12

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(hodge)
