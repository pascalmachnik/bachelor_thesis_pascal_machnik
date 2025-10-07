# TCFit
## TrackCalib fitting library 

Fitting library used in TrackCalib2 package and based on iminuit, numpy, scipy, numba python libraries. 

### Basic usage

PDFs and their parameters can be defined as follows
```pycon
from parameter import Parameter
from basic_pdfs import Gaussian

mass = Parameter(name ="mass", limit = (2800, 3200))
mu = Parameter(name ="mu", value = 3000)
sigma = Parameter(name ="sigma", value = 10)

gaussian = Gaussian(name ="gaussian", observable = observable, mu = mu, sigma = sigma)
```

Unbinned and binned likelihood fit are implemented
```pycon
from likelihood import LikelihoodFit
from extensions import ExtendedPDF
from model import Model
from numpy.random import normal

data = normal(mu.value, sigma.value, 1000)

norm = Parameter(name="norm", value = 1000, limits = (900,1100))
ext_gauss = ExtendedPDF(name="extended_gaussian",pdf=gaussian,normalisation=norm)

model = Model(name="model",pdfs=[ext_gauss])
likelihood = LikelihoodFit(data = [data], model = model)
minuit     = likelihood.minimise(hesse = False, minos = True)
```

