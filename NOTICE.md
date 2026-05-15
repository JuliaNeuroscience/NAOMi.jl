# NOTICE

NAOMi.jl is a Julia port of the **NAOMi-Sim** simulator originally written
in MATLAB by Alex Song and Adam Charles, distributed at
<https://bitbucket.org/adamshch/naomi_sim> under the MIT License
(Copyright 2021 Alex Song, Adam Charles).

The Julia code in this repository is a clean-room reimplementation written
from upstream documentation and algorithmic descriptions; no upstream
`.m` or `.cpp` source has been vendored into this repository. The
algorithms, parameter defaults, and overall pipeline structure are
nevertheless derived directly from the upstream work and authorship credit
for those design decisions belongs to Alex Song and Adam Charles.

When citing NAOMi.jl, please also cite the original publication:

> Song, A., Gauthier, J. L., Pillow, J. W., Tank, D. W., and Charles, A. S.
> (2021). Neural anatomy and optical microscopy (NAOMi) simulation for
> evaluating calcium imaging methods. *Journal of Neuroscience Methods*,
> 358, 109173.

## Upstream MIT License

```
Copyright 2021 Alex Song, Adam Charles

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

NAOMi.jl itself is distributed under the MIT License (see `LICENSE`),
Copyright 2026 Tim Holy and contributors.
