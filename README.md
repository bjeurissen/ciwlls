# Introduction

This repository contains *constrained weighted linear least squares and nonlinear solvers* for problems that can be cast as y = exp(A*x).

This includes, but is not limited to, T2(\*) mapping, Diffusion Tensor Imaging (DTI), Diffusion Kurtosis Imaging (DKI), Q-space Trajectory Imaging (QTI), ...

# Usage

Modality-specific demonstrations can be found at: `demo_{t2,dti,dki,qti}.m`

The general work-flow is demonstrated below using QTI as an example:

1. Reading data (this example uses MRtrix .mif files):
```matlab
[y,vox,v2w,grad] = MRtrix.readToMatlab('dwi.mif'); % DWI data
mask = MRtrix.readToMatlab('mask.mif'); % processing mask
```

2. Applying a mask to avoid calculations in the background:
```matlab
y = Volumes.mask(y,logical(mask)); % sets all voxels outside the mask to NaN
```

3. Creating a model instance. Note that this can be any of T2, DTI, DKI and QTI:
```matlab
model = QTI(grad);
```
Note that the above constructor creates a model instance backed by a constrained iteratively reweighted linear least squares solver with default constraints, which is fine for most applications. The solver/constraints can be modified by supplying additional name/value pairs to the constructor (default in bold):

* to set the type of estimator to be used: <pre><code>'estimator', 'lls/<b>wlls</b>/nls'</code></pre>
* to set the type of initialization for WLLS weights (log of the raw data (data) or unweighted (ones)):  <pre><code>'init_weight', 'data/<b>ones</b>'</code></pre>
* to set the number of iterative reweightings for WLLS: <pre><code>'iter', <b>2</b></code></pre>
* to modify the constraints that are enforced (non-neg diffusivity, non-neg total kurtosis, non-neg isotropic kurtosis, non-neg anisotropic kurtosis, monotonic signal decay): <pre><code>'constraints', [<b>0 0 1 1 1</b>]</code></pre> Note that the vector of supported constraints will depend on the specific model, please consult ``<Modelname>.m`` to obtain a models-specific list.

5. Fitting your model:
```matlab
x = model.solve(y);
```

6. Obtaining scalar metrics:
```matlab
m = model.metrics(x);
```

# Reference

When using the estimators in this repository, please cite the associated paper:

J. Morez, F. Szczepankiewicz, A. J. den Dekker, F. Vanhevel , J. Sijbers, B. Jeurissen. Optimal experimental design and estimation for q-space trajectory imaging. *Hum. Brain Mapp.* doi:10.1002/hbm.26175.
