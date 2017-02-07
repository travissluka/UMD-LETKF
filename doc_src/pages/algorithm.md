title: LETKF Algorithm




The observation operator applied to each ensemble member. This is performed *outside* the LETKF library.
> \( \mathbf{y}^{b(i)} = H \mathbf{x}^{b(i)}\)

The core LETKF algorithm can conceptually be thought of accomplishing two things, first an analysis mean is computed

> \( \tilde{\mathbf{P}}^a = [(k-1)\mathbf{I}+(\mathbf{Y}^b)^T \mathbf{R}^{-1} \mathbf{Y}^b ]^{-1} \)

> \( \bar{\mathbf{w}}^a = \tilde{\mathbf{P}}^a (\mathbf{Y}^b)^T \mathbf{R}^{-1} (\mathbf{y}^o - \bar{\mathbf{y}}^b ) \)

> \( \bar{\mathbf{x}}^a = \bar{\mathbf{x}}^b + \mathbf{X}^b\bar{\mathbf{w}}^a\)

and then the analysis ensemble perturbations \( \mathbf{X}^a \)  are calculated

> \( \mathbf{W}^a = [ (k-1) \tilde{\mathbf{P}}^a ] ^{1/2} \)

> \( \mathbf{X}^a = \mathbf{X}^b\mathbf{W}^a\)