---
layout: default
---

# MFT (Magentic force theory)
Calculate spin exchange coupling parameters J from ([OpenMX](http://www.openmx-square.org/))DFT/Wannier hamiltonians via linear-response theory.
We support the following features.
- J in momentum space.
- Orbital resolved J
- Local axis redefinition for orbital resolved J
- Full [OpenMX](http://www.openmx-square.org/) DFT Hamiltonian and Wannier Hamiltonian ([OpenMX](http://www.openmx-square.org/)/Wannier90)

![Logo](Logo.svg)

# Quick-Start

### Install Julia ([https://julialang.org/](https://julialang.org/))

Julia was designed from the beginning for high performance. Julia programs compile to efficient native code for multiple platforms via LLVM.
[See Julia benchmark compared to C, Fortran, Pythons, Matlab & more...](https://julialang.org/benchmarks/).

 * For Linux system: 

[Using Julia auto-installer for Linux](https://github.com/abelsiqueira/jill)
 ```bash
 JULIA_INSTALL=~/opt/bin bash -ci "$(curl –fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
 ```
 
 * For OSX or Windows see: [Julia Download](https://julialang.org/downloads/)

### Install DFTforge
```bash
git clone https://github.com/ElectronicStructureTheory-KAIST/DFTforge.jl/
julia install.jl
```

### Run example

#### G-type AFM (anti ferromagnetic) NiO example
```bash
./example_NiO_OpenMX.sh
```
---

# DFTforge the DFT postprocessing environment
Simplify obtaining Hamiltonian, Eigenvalue, Eigenvector from DFT results and Pre-caching them for parallelized calculations.
Pre-caching results are easy to use when calculating in k and q space.
 * It consists of two parts DFTforge & DFTrefinery.


## DFTforge
The wrapper for calculating Hamiltonian, Eigenvalue, Eigenvector from DFT results.

 * read DFT results from OpenMX, Wannier90(?), EcalJ(?)
 * Calculate Phi, Enk, Hks,

### Etc. functionalities
 * K-point representation in INT (for unique K,Q).
 * Gennerate K-points.
 * Generalised argument parser (support TOML sytle input).



## DFTrefinery
use DFTforge for calcuating various properties
Wrapper module for easy access of DFT-forge, especially easy use of HDF5 cached eigenstate information.

 * Store Eigenvalues & Eigenvectors in HDF5 format (K,Q)
 * Read Stored Eigenvalues & Eigenvectors
 * Simple interface for K space, K,Q space calucations
 ** Generalised parallelized K,Q space calculation function wrapper.



Text can be **bold**, _italic_, or ~~strikethrough~~.

[Link to another page](./another-page.html).

There should be whitespace between paragraphs.

There should be whitespace between paragraphs. We recommend including a README, or a file with information about your project.

# Header 1

This is a normal paragraph following a header. GitHub is a code hosting platform for version control and collaboration. It lets you and others work together on projects from anywhere.

## Header 2

> This is a blockquote following a header.
>
> When something is important enough, you do it even if the odds are not in your favor.

### Header 3

```js
// Javascript code with syntax highlighting.
var fun = function lang(l) {
  dateformat.i18n = require('./lang/' + l)
  return true;
}
```

```ruby
# Ruby code with syntax highlighting
GitHubPages::Dependencies.gems.each do |gem, version|
  s.add_dependency(gem, "= #{version}")
end
```

#### Header 4

*   This is an unordered list following a header.
*   This is an unordered list following a header.
*   This is an unordered list following a header.

##### Header 5

1.  This is an ordered list following a header.
2.  This is an ordered list following a header.
3.  This is an ordered list following a header.

###### Header 6

| head1        | head two          | three |
|:-------------|:------------------|:------|
| ok           | good swedish fish | nice  |
| out of stock | good and plenty   | nice  |
| ok           | good `oreos`      | hmm   |
| ok           | good `zoute` drop | yumm  |

### There's a horizontal rule below this.

* * *

### Here is an unordered list:

*   Item foo
*   Item bar
*   Item baz
*   Item zip

### And an ordered list:

1.  Item one
1.  Item two
1.  Item three
1.  Item four

### And a nested list:

- level 1 item
  - level 2 item
  - level 2 item
    - level 3 item
    - level 3 item
- level 1 item
  - level 2 item
  - level 2 item
  - level 2 item
- level 1 item
  - level 2 item
  - level 2 item
- level 1 item

### Small image

![Octocat](https://assets-cdn.github.com/images/icons/emoji/octocat.png)

### Large image

![Branching](https://guides.github.com/activities/hello-world/branching.png)


### Definition lists can be used with HTML syntax.

<dl>
<dt>Name</dt>
<dd>Godzilla</dd>
<dt>Born</dt>
<dd>1952</dd>
<dt>Birthplace</dt>
<dd>Japan</dd>
<dt>Color</dt>
<dd>Green</dd>
</dl>

```
Long, single-line code blocks should not wrap. They should horizontally scroll if they are too long. This line should be long enough to demonstrate this.
```

```
The final element.
```
