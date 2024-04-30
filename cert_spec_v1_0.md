# MILP certificate file format version 1.0

## Preliminaries
Before we describe the mixed-integer linear programming (MILP) certificate [file format](#file-format), we introduce some jargon and conventions used in this document.

### Variable indices
Variables of the MILP problem are indexed from $0$ to $n-1$ where $n$ is the total number of variables. Throughout this document, $X_i$ refers to the variable indexed by $i$ regardless of its actual variable name.

### Numbers
We refer to two kinds of numbers: indices and values. Indices are integers. Values are rational numbers expressed as finite decimals or fractions.

### Constraints
The certificate file includes two kinds of constraints: the constraints of the MILP problem and derived constraints. Constraints have unique constraint names and are assigned indices starting from $0$ according to the order in which they appear in the certificate file.

Each constraint is of the form $a_{i_1}X_{i_1}+\dots+a_{i_p}X_{i_p} \text{sense}\ \beta,$ where `sense` is $=,\leq$ or $\geq,\beta$ is a rational value, $p\in\{1,\dots,n\},$ and for $j=1,\dots p,i_j$ is a variable index, $X_{i_j}$ is a variable with index $i_j,$ and $a_{i_j}$ is a rational value.

Bound constraints such as $0\leq\beta\leq 1$ are not treated specially.

### Rounding
Suppose that the constraint $a_{i_1}X_{i_1}+\dots+a_{i_p}X_{i_p} \text{sense}\ \beta,$ where `sense` is $\leq$ or $\geq,$ is such that fpr $j=1,\dots, p,a_{i_j}$
is and integer and is nonzero only if $i_j$ is an integer variable index, then the constraint $a_{i_1}X_{i_1}+\dots+a_{i_p}X_{i_p} \text{sense}\ \beta'$ with

$$
\beta' =
\left \{
\begin{array}{ll}
\lfloor\beta\rfloor & \text{if sense is } \leq \\
\lceil\beta\rceil & \text{if sense is } \geq \\
\end{array}
\right.
$$

is said to be obtained from rounding.

### Domination of constraints
We use a restricted form of constraint domination. A constraint that is $0\geq \beta$ with $\beta >0$ , $0\leq\beta$ with $\beta<0$, or $0=\beta$ with $\beta\neq 0$ is called an absurdity or falsehood. Such a constraint dominates any other constraint.

The constraint $a^Tx\geq\beta$ or $a^Tx=\beta$ dominates $a'^Tx\geq\beta′$ if $a=a′$ and $\beta\geq\beta'.$

Similarly, the constraint $a^Tx\leq\beta$ or $a^Tx=\beta$ dominates $a′^Tx\leq\beta′$ if $a=a′$ and $\beta\leq\beta'.$

Finally, the constraint $a^Tx=\beta$ dominates $a′^Tx=\beta′$ if $a=a′$ and $\beta=\beta′.$

### Suitable linear combinations of constraints
For a constraint named $C$, define $l(C)$ to be the left-hand side of $C$ and define $r(C)$ to be the right-hand side of $C$, and

$$
s(C)=\begin{cases}
1 & \text{if } C \text{ is a } \geq\text{-constraint}\\
0 & \text{if } C\text{ is a } =\text{-constraint} \\
-1 & \text{if } C\text{ is a } \leq\text{-constraint}
\end{cases}
$$

For example, for the constraint $C_1:  2X_1+3X_2\geq 1$, we have $l(C_1)=2X_1+3X_2,r(C_1)=1,s(C_1)=1.$
Let $C_1,\dots,C_k$ be constraint names. Let $\lambda_1,\dots,\lambda_k$ be values such that $\lambda_js(C_j)\geq 0$ for all $j=1,…,k,$ or $\lambda_js(C_j)\leq 0$ for all $j=1,…,k.$
Then $\lambda_1\cdot C_1+\dots +\lambda_k \cdot C_k$ is called a suitable linear combination and denotes the constraint

$$
\sum_i\lambda_i l(C_i)=\sum_i\lambda_i r(C_i) \text{ if }\lambda_j s(C_j) = 0 \text{ for all } j=1,\dots,k
$$

$$
\sum_i\lambda_i l(C_i)\geq\sum_i\lambda_i r(C_i) \text{ if }\lambda_j s(C_j) \geq 0 \text{ for all } j=1,\dots,k \text{ and } \lambda_q s(C_q)\neq 0 \text{ for some } q\in\{1,\dots,k\}
$$

$$
\sum_i\lambda_i l(C_i)\leq\sum_i\lambda_i r(C_i) \text{ if }\lambda_j s(C_j) \leq 0 \text{ for all } j=1,\dots,k \text{ and } \lambda_q s(C_q)\neq 0 \text{ for some } q\in\{1,\dots,k\}
$$

## Constraint format
Specifying a constraint $a_{i_1}X_{i_1}+\cdots+a_{i_p}X_{i_p}\text{ sense }\beta$ in constraint format means listing the details of the constraint in the following order: the constraint name, then the character representing 𝚜𝚎𝚗𝚜𝚎 (`E` when `sense` is $=,$ `L` when `sense` is $\leq,$ `G` when `sense` is $\geq$), then $\beta$, then
- $p$ $i_1$ $a_{i_1}$ $\cdots$ $i_p$ $a_{i_p}$;
- or the keyword `OBJ` for indicating that $a_{i_1}X_{i_1}+⋯+a_{i_p}X_{i_p}$ is the same as the objective function.

For example, the constraint $−1X_1+6X_2\leq3$ with constraint name `C3` specified in constraint format is

```
C3 L 3  2  1 -1  2 6
```

If the objective function specified in the certificate file is equal to $−1X_1+6X_2,$ no matter whether to be minimized or maximized, then we may also write the shorter form

```
C3 L 3 OBJ
```

### Reason format
Every constraint in the [DER section](#der-section) specified below is associated with a `reason`.
A `reason` associated with a constraint with constraint index `idx` must have one of the following forms:

- { asm }
- {  lin  $p$  $i_1$  $\lambda_1$ $\cdots$  $i_p$  $\lambda_p$  } where $p$ is a nonnegative integer, $i_1,\dots,i_p$ are distinct indices at least $0$ and less than `idx`, and $\lambda_1,\dots,\lambda_p$ are rational values such that $\lambda_1\cdot C_{i_1}+\dots+\lambda_p\cdot C_{i_p}$ is a suitable linear combination that dominates the associated constraint; here $C_{i_j}$ denotes the constraint with index $i_j$ for $j=1,\dots,p.$
- {  rnd  $p$  $i_1$  $\lambda_1$  $\cdots$  $i_p$  $\lambda_p$  } where $p$ is a nonnegative integer, $i_1,\dots,i_p$ are distinct indices at least $0$ and less than `idx`, and $\lambda_1,\dots,\lambda_p$ are rational values such that $\lambda_1\cdot C_{i_1}+⋯+\lambda_p\cdot C_{i_p}$, is a suitable linear combination that can be rounded and the rounded constraint dominates the associated constraint; again $C_{i_j}$ denotes the constraint with index $i_j$ for $j=1,…,p.$
- {  uns  $i_1$  $l_1$  $i_2$  $l_2$  } where $i_1, l_1, i_2, l_2$ are indices at least $0$ and less than `idx` such that the constraints indexed by $i_1$ and $i_2$ both dominate the associated constraint and the constraints indexed by $l_1$ and $l_2$ are, perhaps after reordering, $a_{i_1}X_{i_1}+\cdots a_{i_p}X_{i_p}\leq\beta$ and $a_{i_1}X_{i_1}+\cdots a_{i_p}X_{i_p}\geq\beta +1$ such that $\beta$ is an integer and for $j=1,…,p, a_{i_j}$ is an integer and $X_{i_j}$ is an integer variable.
- { sol }

## File format
The initial lines of the file may be comment lines. Each comment line must begin with the character "%".
The first line, after the comment lines, must be
```
VER 1.0
```
indicating that the certificate file conforms to version 1.0 of the certificate file format.

The remaining content is divided into seven sections and must appear in the following order:
-   `VAR` for the variables,
-   `INT` for integer variables,
-   `OBJ` for the objective function,
-   `CON` for the constraints,
-   `RTP` for the relation to prove,
-   `SOL` for the solutions to check,
-   `DER` for the derivations.

These sections have to be formatted as shown below.

### VAR section
The section begins with
```
VAR n
```
where $n$ is the number of variables, then followed by the n variable names separated by spaces or line breaks.
The variables are assigned indices from $0$ to $n−1$ according to the order in which the variable names are listed.
For example, the following specifies two variables `x` and `y`:
```
VAR 2
x y
```

### INT section

The section begins with
`INT i`
where $i$ is the number of integer variables, then followed by the $i$ integer variable indices separated by spaces or line breaks.
For example, the following specifies our variables `x` and `y` as integer variables:
```
INT 2
0 1
```

### `OBJ` section

The section begins with
`OBJ objsense`
where `objsense` is the keyword `min` for minimization or `max` for maximization, then followed by an integer $k$ and $2k$ numbers $i_1\ c_1\ i_2\ c_2\ \dots\ i_k\ c_k$
separated by spaces or line breaks where for $j=1,…,k, i_j$ is a variable index and $c_k$ is the objective function coefficient for the variable with index $i_j.$

For example, the `OBJ` section for the problem

$$
\begin{array}{ll}
\text{min} & x+y \\
\text{s.t.} & C_1:4x+y\geq 1 \\
 & C_2:4x-y\leq 2\\
\end{array}
$$

could look like the following
```
OBJ min
2  0 1  1 1
```

### CON section
The section begins with

`CON m  b`
where $m$ is the total number of constraints (including bound constraints) and $b$ is the total number of bound constraints, and then followed by $m$ constraints listed in [constraint format](#constraint-format).

The constraints in this section are assigned indices from $0$ to $m−1$ according to the order in which they are listed.
Bound constraints must appear at the beginning of the section **before** the nonbound constraints.

For example, the `CON` section for the problem

$$
\begin{array}{ll}
\text{min} & x+y \\
\text{s.t.} & C_1:4x+y\geq 1 \\
 & C_2:4x-y\leq 2\\
\end{array}
$$

could look like the following
```
CON 2 0
C1 G 1  2  0 4  1 1
C2 L 2  2  0 4  1 -1
```
Here, there are no bound constraints but any program that outputs a certificate is free to choose arbitrary names for the bound constraints.

### RTP section
The section is either

`RTP infeas`
for indicating infeasibity, or
`RTP range lb  ub`
where `lb` specifies the lower bound on the optimal value to be verified and `ub` specifies the upper bound on the optimal value to be verified. The specified lower bound can be `-inf` if no actual lower bound is to be verified and the specified upper bound can be `inf` if no actual upper bound is to be verified.
For example, if the `RTP` section is
```
RTP range 1 inf
```
only a lower bound of $1$ is to be verified.
In our example, the `RTP` section looks like this:
```
RTP range 1 1
```
### SOL section
The section begins with
`SOL s`
where s is the number of solutions to be verified for feasibility, and then followed by
$S_1 \\ \vdots \\ S_s$
such that for $j=1,…,s, S_j$ is of the form $\text{name } p\ i_1\ v_1\ \cdots\ i_p\ v_p,$ where `name` is the solution name, $p$ is a nonnegative integer, and for $r=1,…,p, i_r$ is a variable index and $v_r$ is the solution value for the variable with index $i_r$. All other variables are assumed to have the value zero.
At least one of the solutions specified should have objective function value at least the lower bound given in the [RTP section](#rtp-section) in the case of maximization and at most the upper bound given in the [RTP section](rtp-section) in the case of minimization.
For example, to specify the two solutions $(x,y)=(1,2)$ (which is feasible) and $(x,y)=(0,1)$ (which is optimal) where $x$ has variable index `0` and $y$ has variable index `1`, the `SOL` section could look like this:
```
SOL 2
feas 2  0 1  1 2
opt 1  1 1
```

### DER section
The section begins with
`DER d`
where `d` is the number of derived constraints to be verified, and then followed by
$D_1\\\vdots\\ D_d$
such that $D_1,\dots,D_d$ are assigned constraint indices from $m$ to $m+d−1$ according to the given order and for $j=1,…,d, D_j$ is of the form
```
constraint reason index
```
where `constraint` is the constraint to be verified in [constraint format](#constraint-format), `reason` is the reason associated with the constraint in [reason format](#reason-format), and `index` is $−1$ or is an index such that no constraint with a larger index refers to this constraint.
Each constraint in this section carries implicitly a set of assumptions deduced as follows:

- If `reason` is {  asm  }, then the set of assumptions contains simply the index of the constraint.
- If `reason` is {  lin  $p\  i_1\  \lambda_1\  \cdots \ i_p\   \lambda_p$  } or {  rnd  $p\  i_1\  \lambda_1\  \cdots\  i_p\  \lambda_p$  }, then the set of assumptions is the union of the sets of assumptions of the constraints indexed by $i_1,…,i_p.$
- If `reason` is {  uns  $i_1\  l_1\  i_2\  l_2$  }, then the set of assumptions is the union of the sets of assumptions of the constraint indexed by $i_1$ without $l_1$ and of the constraint indexed by $i_2$ without $l_2.$

If the `RTP` section is `RTP infeas`, then the last constraint in this section should be an absurdity with an empty set of assumptions.
If the `RTP` section is `RTP range lb ub`, then the last constraint in this section should be a constraint with an empty set of assumptions that dominates $\text{OBJ}\geq  lb$ in the case of minimization or $\text{OBJ}\leq ub$ in the case of maximization where `OBJ` denotes the objective function.

**Remark.** The reason for specifying `index` is to allow a verifier to discard the corresponding constraint from memory once the verifier has completed verifying constraint with index `index`.

For example, the `DER` section for the problem above could look like this:
```
C3 G -1/2  1  1 1   { lin 2  0 1/2  1 -1/2 } 3
C4 G 0     1  1 1   { rnd 1  2 1 } 4
C5 G 1/4   OBJ     { lin 2  0 1/4  3 3/4 } 5
C6 G 1     OBJ     { rnd 1  4 1 } 0
```

For a small example, we refer to the file [paper_eg3.vipr](code/paper_eg3.vipr).
