# NURBS-BOOK-ALGORITHMS-JS

Implementations of NURBS book algorithms from pseudo code to JS.

## Overview

This library provides JavaScript implementations of various algorithms from the NURBS book. The algorithms cover a wide range of topics including curve and surface basics, B-spline basis functions, B-spline curves and surfaces, rational B-spline curves and surfaces, fundamental geometric algorithms, advanced geometric algorithms, conics and circles, construction of common surfaces, curve and surface fitting, and advanced surface construction techniques.

## Table of Contents

1. [Curve and Surface Basics](#curve-and-surface-basics)
2. [B-Spline Basis Functions](#b-spline-basis-functions)
3. [B-Spline Curves and Surfaces](#b-spline-curves-and-surfaces)
4. [Rational B-Spline Curves and Surfaces](#rational-b-spline-curves-and-surfaces)
5. [Fundamental Geometric Algorithms](#fundamental-geometric-algorithms)
6. [Advanced Geometric Algorithms](#advanced-geometric-algorithms)
7. [Conics and Circles](#conics-and-circles)
8. [Construction of Common Surfaces](#construction-of-common-surfaces)
9. [Curve and Surface Fitting](#curve-and-surface-fitting)
10. [Advanced Surface Construction Techniques](#advanced-surface-construction-techniques)

## Curve and Surface Basics

### Horner's Algorithm for Power Basis Curve

**Function:** `Horner1(a, n, u0)`

**Description:** Computes a point on a power basis curve using Horner's algorithm.

**Parameters:**
- `a`: Array of coefficients.
- `n`: Degree of the curve.
- `u0`: Parameter value.

**Returns:** Computed point on the curve.

**Example:**
```javascript
let a = [/* coefficients array */];
let n = a.length - 1;
let u0 = /* value of u */;
let point = Horner1(a, n, u0);
console.log(point);
```

### Bernstein's Algorithm for Bernstein Polynomial

**Function:** `Bernstein(i, n, u)`

**Description:** Computes the value of a Bernstein polynomial.

**Parameters:**
- `i`: Index of the polynomial.
- `n`: Degree of the polynomial.
- `u`: Parameter value.

**Returns:** Value of the Bernstein polynomial.

**Example:**
```javascript
let i = /* index */;
let n = /* degree */;
let u = /* value of u */;
let B = Bernstein(i, n, u);
console.log(B);
```

### Algorithm for Computing All nth-Degree Bernstein Polynomials

**Function:** `AllBernstein(n, u)`

**Description:** Computes all nth-degree Bernstein polynomials.

**Parameters:**
- `n`: Degree of the polynomials.
- `u`: Parameter value.

**Returns:** Array of Bernstein polynomial values.

**Example:**
```javascript
let n = /* degree */;
let u = /* value of u */;
let Bs = AllBernstein(n, u);
console.log(Bs);
```

### Point on Bezier Curve Using Bernstein Polynomials

**Function:** `PointOnBezierCurve(P, n, u)`

**Description:** Computes a point on a Bezier curve using Bernstein polynomials.

**Parameters:**
- `P`: Array of control points.
- `n`: Degree of the curve.
- `u`: Parameter value.

**Returns:** Computed point on the Bezier curve.

**Example:**
```javascript
let P = [/* array of control points */];
let n = P.length - 1;
let u = /* value of u */;
let pointOnCurve = PointOnBezierCurve(P, n, u);
console.log(pointOnCurve);
```

### Point on Bezier Curve Using de Casteljau's Algorithm

**Function:** `deCasteljau1(P, n, u)`

**Description:** Computes a point on a Bezier curve using de Casteljau's algorithm.

**Parameters:**
- `P`: Array of control points.
- `n`: Degree of the curve.
- `u`: Parameter value.

**Returns:** Computed point on the Bezier curve.

**Example:**
```javascript
let P = [/* array of control points */];
let n = P.length - 1;
let u = /* value of u */;
let pointOnCurve = deCasteljau1(P, n, u);
console.log(pointOnCurve);
```

### Point on Power Basis Surface

**Function:** `Horner2(a, n, m, u0, v0)`

**Description:** Computes a point on a power basis surface.

**Parameters:**
- `a`: Array of arrays representing the surface control points.
- `n`: Degree in the u direction.
- `m`: Degree in the v direction.
- `u0`: Parameter value in the u direction.
- `v0`: Parameter value in the v direction.

**Returns:** Computed point on the power basis surface.

**Example:**
```javascript
let a = [/* array of arrays representing the surface control points */];
let n = /* the degree in the u direction */;
let m = /* the degree in the v direction */;
let u0 = /* value of u */;
let v0 = /* value of v */;
let surfacePoint = Horner2(a, n, m, u0, v0);
console.log(surfacePoint);
```

### Point on Bezier Surface Using de Casteljau's Algorithm

**Function:** `deCasteljau2(P, n, m, u0, v0)`

**Description:** Computes a point on a Bezier surface using de Casteljau's algorithm.

**Parameters:**
- `P`: Array of arrays of control points.
- `n`: Degree in the u direction.
- `m`: Degree in the v direction.
- `u0`: Parameter value in the u direction.
- `v0`: Parameter value in the v direction.

**Returns:** Computed point on the Bezier surface.

**Example:**
```javascript
let P = [/* array of arrays representing the surface control points */];
let n = /* the degree in the u direction */;
let m = /* the degree in the v direction */;
let u0 = /* value of u */;
let v0 = /* value of v */;
let surfacePoint = deCasteljau2(P, n, m, u0, v0);
console.log(surfacePoint);
```

## B-Spline Basis Functions

### Knot Span Index in NURBS Calculations

**Function:** `FindSpan(n, p, u, U)`

**Description:** Determines the knot span index in NURBS calculations.

**Parameters:**
- `n`: Number of knots minus one.
- `p`: Degree of the B-spline.
- `u`: Parameter value.
- `U`: Knot vector.

**Returns:** Knot span index.

**Example:**
```javascript
let n = /* the number of knots minus one */;
let p = /* the degree of the B-spline */;
let u = /* the parameter value */;
let U = [/* the knot vector */];
let span = FindSpan(n, p, u, U);
console.log(span);
```

### Nonvanishing Basis Functions

**Function:** `BasisFuns(i, p, u, U)`

**Description:** Computes the nonvanishing basis functions.

**Parameters:**
- `i`: Knot span.
- `p`: Degree of the B-spline.
- `u`: Parameter value.
- `U`: Knot vector.

**Returns:** Array of basis functions.

**Example:**
```javascript
let i = /* the knot span */;
let p = /* the degree of the B-spline */;
let u = /* the parameter value */;
let U = [/* the knot vector */];
let basisFuns = BasisFuns(i, p, u, U);
console.log(basisFuns);
```

### Span and Multiplicity of a Knot

**Function:** `findSpanMult(n, p, u, UP)`

**Description:** Finds the span and multiplicity of a knot in a knot vector.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `u`: Parameter value.
- `UP`: Knot vector.

**Returns:** Array containing the span and multiplicity of the knot.

**Example:**
```javascript
let n = /* number of control points minus 1 */;
let p = /* degree of curve */;
let u = /* parameter value */;
let UP = [/* knot vector */];
let [span, mult] = findSpanMult(n, p, u, UP);
console.log('Span:', span, 'Multiplicity:', mult);
```

### Nonvanishing Basis Functions and Their Derivatives

**Function:** `DersBasisFuns(i, p, u, U, n)`

**Description:** Computes the nonvanishing basis functions and their derivatives.

**Parameters:**
- `i`: Knot span.
- `p`: Degree of the B-spline.
- `u`: Parameter value.
- `U`: Knot vector.
- `n`: Order of the derivative.

**Returns:** Array of derivatives of the basis functions.

**Example:**
```javascript
let i = /* the knot span */;
let p = /* the degree of the B-spline */;
let u = /* the parameter value */;
let U = [/* the knot vector */];
let n = /* the order of the derivative */;
let ders = DersBasisFuns(i, p, u, U, n);
console.log(ders);
```

### Basis Function Nip

**Function:** `OneBasisFun(p, m, U, i, u)`

**Description:** Computes the basis function Nip.

**Parameters:**
- `p`: Degree of the B-spline.
- `m`: Upper index of the knot vector.
- `U`: Knot vector.
- `i`: Knot span.
- `u`: Parameter value.

**Returns:** Value of the basis function.

**Example:**
```javascript
let p = /* the degree of the B-spline */;
let m = /* the upper index of U */;
let U = [/* the knot vector */];
let i = /* the knot span */;
let u = /* the parameter value */;
let Nip = OneBasisFun(p, m, U, i, u);
console.log(Nip);
```

### Derivatives of Basis Function Nip

**Function:** `DersOneBasisFun(p, m, U, i, u, n)`

**Description:** Computes the derivatives of the basis function Nip.

**Parameters:**
- `p`: Degree of the B-spline.
- `m`: Upper index of the knot vector.
- `U`: Knot vector.
- `i`: Knot span.
- `u`: Parameter value.
- `n`: Derivative order.

**Returns:** Array of derivatives of the basis function.

**Example:**
```javascript
let p = /* the degree of the B-spline */;
let m = /* the upper index of U */;
let U = [/* the knot vector */];
let i = /* the knot span */;
let u = /* the parameter value */;
let n = /* the derivative order */;
let derivatives = DersOneBasisFun(p, m, U, i, u, n);
console.log(derivatives);
```

## B-Spline Curves and Surfaces

### Curve Point

**Function:** `CurvePoint(n, p, U, P, u)`

**Description:** Computes a point on a B-spline curve.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `P`: Array of control points.
- `u`: Parameter value.

**Returns:** Computed point on the curve.

**Example:**
```javascript
let n = /* the number of control points minus 1 */;
let p = /* the degree of the curve */;
let U = [/* the knot vector */];
let P = [/* the array of control points */];
let u = /* the parameter value */;
let curvePoint = CurvePoint(n, p, U, P, u);
console.log(curvePoint);
```

### Curve Derivatives

**Function:** `CurveDerivsAlg1(n, p, U, P, u, d)`

**Description:** Computes the derivatives of a B-spline curve.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `P`: Array of control points.
- `u`: Parameter value.
- `d`: Derivative order.

**Returns:** Array of derivatives of the curve.

**Example:**
```javascript
let n = /* the number of control points minus 1 */;
let p = /* the degree of the curve */;
let U = [/* the knot vector */];
let P = [/* the array of control points */];
let u = /* the parameter value */;
let d = /* the derivative order */;
let curveDerivatives = CurveDerivsAlg1(n, p, U, P, u, d);
console.log(curveDerivatives);
```

### Control Points of Curve Derivatives

**Function:** `CurveDerivCpts(n, p, U, P, d, r1, r2)`

**Description:** Computes the control points of the derivatives of a B-spline curve.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `P`: Array of control points.
- `d`: Derivative order.
- `r1`: Lower index of the range of control points.
- `r2`: Upper index of the range of control points.

**Returns:** Array of control points for the curve derivatives.

**Example:**
```javascript
let n = /* the number of control points minus 1 */;
let p = /* the degree of the curve */;
let U = [/* the knot vector */];
let P = [/* the control points as an array of [x, y, z] points */];
let d = /* the derivative order */;
let r1 = /* the lower index of the range of control points */;
let r2 = /* the upper index of the range of control points */;
let derivativeControlPoints = CurveDerivCpts(n, p, U, P, d, r1, r2);
console.log(derivativeControlPoints);
```

### Curve Derivatives (Alternative Algorithm)

**Function:** `CurveDerivsAlg2(n, p, U, P, u, d)`

**Description:** Computes the derivatives of a B-spline curve using an alternative algorithm.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `P`: Array of control points.
- `u`: Parameter value.
- `d`: Derivative order.

**Returns:** Array of derivatives of the curve.

**Example:**
```javascript
let n = /* the number of control points minus 1 */;
let p = /* the degree of the curve */;
let U = [/* the knot vector */];
let P = [/* the array of control points as [x, y, z] */];
let u = /* the parameter value */;
let d = /* the derivative order */;
let curveDerivatives = CurveDerivsAlg2(n, p, U, P, u, d);
console.log(curveDerivatives);
```

### Surface Point

**Function:** `SurfacePoint(n, p, U, m, q, V, P, u, v)`

**Description:** Computes a point on a B-spline surface.

**Parameters:**
- `n`: Number of control points in the u direction minus one.
- `p`: Degree of the surface in the u direction.
- `U`: Knot vector in the u direction.
- `m`: Number of control points in the v direction minus one.
- `q`: Degree of the surface in the v direction.
- `V`: Knot vector in the v direction.
- `P`: Array of arrays of control points.
- `u`: Parameter value in the u direction.
- `v`: Parameter value in the v direction.

**Returns:** Computed point on the surface.

**Example:**
```javascript
let n = /* number of control points in u direction minus 1 */;
let p = /* degree of the surface in u direction */;
let U = [/* knot vector in u direction */];
let m = /* number of control points in v direction minus 1 */;
let q = /* degree of the surface in v direction */;
let V = [/* knot vector in v direction */];
let P = [/* control point grid as an array of arrays of [x, y, z] points */];
let u = /* parameter value in u direction */;
let v = /* parameter value in v direction */;
let surfacePoint = SurfacePoint(n, p, U, m, q, V, P, u, v);
console.log(surfacePoint);
```

### Surface Derivatives

**Function:** `SurfaceDerivsAlg1(n, p, U, m, q, V, P, u, v, d)`

**Description:** Computes the derivatives of a B-spline surface.

**Parameters:**
- `n`: Number of control points in the u direction minus one.
- `p`: Degree of the surface in the u direction.
- `U`: Knot vector in the u direction.
- `m`: Number of control points in the v direction minus one.
- `q`: Degree of the surface in the v direction.
- `V`: Knot vector in the v direction.
- `P`: Array of arrays of control points.
- `u`: Parameter value in the u direction.
- `v`: Parameter value in the v direction.
- `d`: Order of derivatives.

**Returns:** Array of derivatives of the surface.

**Example:**
```javascript
let n = /* number of control points in u direction minus 1 */;
let p = /* degree of the surface in u direction */;
let U = [/* knot vector in u direction */];
let m = /* number of control points in v direction minus 1 */;
let q = /* degree of the surface in v direction */;
let V = [/* knot vector in v direction */];
let P = [/* control point grid as an array of arrays of [x, y, z] points */];
let u = /* parameter value in u direction */;
let v = /* parameter value in v direction */;
let d = /* order of derivatives */;
let surfaceDerivatives = SurfaceDerivsAlg1(n, p, U, m, q, V, P, u, v, d);
console.log(surfaceDerivatives);
```

### Control Points of Derivative Surfaces

**Function:** `SurfaceDerivCpts(n, p, U, m, q, V, P, d, r1, r2, s1, s2)`

**Description:** Computes the control points of the derivatives of a B-spline surface.

**Parameters:**
- `n`: Number of control points in the u direction minus one.
- `p`: Degree of the surface in the u direction.
- `U`: Knot vector in the u direction.
- `m`: Number of control points in the v direction minus one.
- `q`: Degree of the surface in the v direction.
- `V`: Knot vector in the v direction.
- `P`: Array of arrays of control points.
- `d`: Order of derivatives.
- `r1`: Lower index of the range of control points in the u direction.
- `r2`: Upper index of the range of control points in the u direction.
- `s1`: Lower index of the range of control points in the v direction.
- `s2`: Upper index of the range of control points in the v direction.

**Returns:** Array of control points for the surface derivatives.

**Example:**
```javascript
let n = /* the number of control points in u direction minus 1 */;
let p = /* the degree in u direction */;
let U = [/* the knot vector in u direction */];
let m = /* the number of control points in v direction minus 1 */;
let q = /* the degree in v direction */;
let V = [/* the knot vector in v direction */];
let P = [/* the control point grid */];
let d = /* the derivative order */;
let r1 = /* the start index in u direction */;
let r2 = /* the end index in u direction */;
let s1 = /* the start index in v direction */;
let s2 = /* the end index in v direction */;
let derivativeControlPoints = SurfaceDerivCpts(n, p, U, m, q, V, P, d, r1, r2, s1, s2);
console.log(derivativeControlPoints);
```

### Surface Derivatives (Alternative Algorithm)

**Function:** `SurfaceDerivsAlg2(n, p, U, m, q, V, P, u, v, d)`

**Description:** Computes the derivatives of a B-spline surface using an alternative algorithm.

**Parameters:**
- `n`: Number of control points in the u direction minus one.
- `p`: Degree of the surface in the u direction.
- `U`: Knot vector in the u direction.
- `m`: Number of control points in the v direction minus one.
- `q`: Degree of the surface in the v direction.
- `V`: Knot vector in the v direction.
- `P`: Array of arrays of control points.
- `u`: Parameter value in the u direction.
- `v`: Parameter value in the v direction.
- `d`: Order of derivatives.

**Returns:** Array of derivatives of the surface.

**Example:**
```javascript
let n = /* the number of control points in u direction minus 1 */;
let p = /* the degree in u direction */;
let U = [/* the knot vector in u direction */];
let m = /* the number of control points in v direction minus 1 */;
let q = /* the degree in v direction */;
let V = [/* the knot vector in v direction */];
let P = [/* the control point grid */];
let u = /* parameter value in u direction */;
let v = /* parameter value in v direction */;
let d = /* the order of derivatives */;
let surfaceDerivatives = SurfaceDerivsAlg2(n, p, U, m, q, V, P, u, v, d);
console.log(surfaceDerivatives);
```

## Rational B-Spline Curves and Surfaces

### Point on Rational B-Spline Curve

**Function:** `CurvePoint(n, p, U, Pw, u)`

**Description:** Computes a point on a rational B-spline curve.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `Pw`: Array of control points with weights.
- `u`: Parameter value.

**Returns:** Computed point on the curve.

**Example:**
```javascript
let n = /* the number of control points minus 1 */;
let p = /* the degree of the curve */;
let U = [/* the knot vector */];
let Pw = [/* the control points with weights as an array of [x, y, z, w] */];
let u = /* the parameter value */;
let curvePoint = CurvePoint(n, p, U, Pw, u);
console.log(curvePoint);
```

### Derivatives of Rational B-Spline Curve

**Function:** `RatCurveDerivs(Aders, wders, d)`

**Description:** Computes the derivatives of a rational B-spline curve from the derivatives of the weighted points.

**Parameters:**
- `Aders`: Array of derivatives of the weighted points.
- `wders`: Array of derivatives of the weights.
- `d`: Derivative order.

**Returns:** Array of derivatives of the curve points in Cartesian coordinates.

**Example:**
```javascript
let Aders = [/* array of derivatives of the weighted points */];
let wders = [/* array of derivatives of the weights */];
let d = /* the derivative order */;
let curveDerivatives = RatCurveDerivs(Aders, wders, d);
console.log(curveDerivatives);
```

### Point on Rational B-Spline Surface

**Function:** `SurfacePoint(n, p, U, m, q, V, Pw, u, v)`

**Description:** Computes a point on a rational B-spline surface.

**Parameters:**
- `n`: Number of control points in the u direction minus one.
- `p`: Degree of the surface in the u direction.
- `U`: Knot vector in the u direction.
- `m`: Number of control points in the v direction minus one.
- `q`: Degree of the surface in the v direction.
- `V`: Knot vector in the v direction.
- `Pw`: Array of arrays of control points with weights.
- `u`: Parameter value in the u direction.
- `v`: Parameter value in the v direction.

**Returns:** Computed point on the surface.

**Example:**
```javascript
let n = /* the number of control points in u direction minus 1 */;
let p = /* the degree of the surface in u direction */;
let U = [/* the knot vector in u direction */];
let m = /* the number of control points in v direction minus 1 */;
let q = /* the degree of the surface in v direction */;
let V = [/* the knot vector in v direction */];
let Pw = [/* the weighted control points as an array of [x, y, z, w] */];
let u = /* the parameter value in u direction */;
let v = /* the parameter value in v direction */;
let surfacePoint = SurfacePoint(n, p, U, m, q, V, Pw, u, v);
console.log(surfacePoint);
```

### Derivatives of Rational B-Spline Surface

**Function:** `RatSurfaceDerivs(Aders, wders, d)`

**Description:** Computes the derivatives of a rational B-spline surface from the derivatives of the weighted points.

**Parameters:**
- `Aders`: Array of derivatives of the weighted points.
- `wders`: Array of derivatives of the weights.
- `d`: Derivative order.

**Returns:** Array of derivatives of the surface points in Cartesian coordinates.

**Example:**
```javascript
let Aders = [/* array of derivatives of the weighted points */];
let wders = [/* array of derivatives of the weights */];
let d = /* the derivative order */;
let surfaceDerivatives = RatSurfaceDerivs(Aders, wders, d);
console.log(surfaceDerivatives);
```

## Fundamental Geometric Algorithms

### Knot Insertion on a NURBS Curve

**Function:** `CurveKnotIns(n, p, UP, Pw, u, k, s, r)`

**Description:** Performs knot insertion on a NURBS curve.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `UP`: Original knot vector.
- `Pw`: Original control points.
- `u`: Knot to insert.
- `k`: Knot span.
- `s`: Multiplicity of the knot.
- `r`: Number of times to insert the knot.

**Returns:** Object containing the new number of control points, new knot vector, and new control points.

**Example:**
```javascript
let n = /* number of control points minus 1 */;
let p = /* degree of curve */;
let UP = [/* original knot vector */];
let Pw = [/* original weighted control points as an array of {x, y, z, w} */];
let u = /* knot to insert */;
let k = /* knot span */;
let s = /* multiplicity of knot u */;
let r = /* number of times to insert knot u */;
let { nq, UQ, Qw } = CurveKnotIns(n, p, UP, Pw, u, k, s, r);
console.log(UQ); // New knot vector
console.log(Qw); // New control points
```

### Point on Rational B-Spline Curve by Corner Cutting

**Function:** `curvePntByCornerCut(np, UP, w, u)`

**Description:** Computes a point on a rational B-spline curve using corner cutting.

**Parameters:**
- `np`: Array of control points.
- `UP`: Knot vector.
- `w`: Weight.
- `u`: Parameter value.

**Returns:** Computed point on the curve.

**Example:**
```javascript
let np = [/* control points as an array of {x, y, z, w} */];
let UP = [/* knot vector */];
let w = /* weight */;
let u = /* parameter value */;
let C = curvePntByCornerCut(np, UP, w, u);
console.log(C); // Computed point on the curve
```

### Surface Knot Insertion in NURBS

**Function:** `SurfaceKnotIns(np, p, UP, mp, q, VP, Pw, dir, uv, k, s, r, nq, UQ, mq, VQ, Qw)`

**Description:** Inserts a knot into a NURBS surface in either the U or V direction.

**Parameters:**
- `np`: Number of control points in the U direction minus one.
- `p`: Degree in the U direction.
- `UP`: Knot vector in the U direction.
- `mp`: Number of control points in the V direction minus one.
- `q`: Degree in the V direction.
- `VP`: Knot vector in the V direction.
- `Pw`: Array of control points.
- `dir`: Direction for knot insertion (0 for U, 1 for V).
- `uv`: Knot value to insert.
- `k`: Span where the knot is to be inserted.
- `s`: Multiplicity of the knot.
- `r`: Number of times the knot is to be inserted.
- `nq`: New number of control points in the U (or V) direction minus one.
- `UQ`: New knot vector in the U (or V) direction.
- `mq`: New number of control points in the V (or U) direction minus one.
- `VQ`: New knot vector in the V (or U) direction.
- `Qw`: New control points array after insertion.

**Returns:** Object containing the new control points and knot vectors.

**Example:**
```javascript
let np = /* number of control points in U minus one */;
let p = /* degree in U direction */;
let UP = /* knot vector in U */;
let mp = /* number of control points in V minus one */;
let q = /* degree in V direction */;
let VP = /* knot vector in V */;
let Pw = /* control points array */;
let dir = /* direction for knot insertion (0 for U, 1 for V) */;
let uv = /* knot value to insert */;
let k = /* span where knot is to be inserted */;
let s = /* multiplicity of knot */;
let r = /* number of times knot is to be inserted */;
let nq = /* new number of control points in U (or V) minus one */;
let UQ = /* new knot vector in U (or V) */;
let mq = /* new number of control points in V (or U) minus one */;
let VQ = /* new knot vector in V (or U) */;
let Qw = /* new control points array after insertion */;
let result = SurfaceKnotIns(np, p, UP, mp, q, VP, Pw, dir, uv, k, s, r, nq, UQ, mq, VQ, Qw);
console.log('New Knot Vector:', result.UQ || result.VQ);
console.log('New Control Points:', result.Qw);
```

### Knot Refinement in NURBS Curve

**Function:** `RefineKnotVectCurve(n, p, U, Pw, X, r, Ubar, Qw)`

**Description:** Refines a NURBS curve by inserting a given set of knots.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Original knot vector.
- `Pw`: Array of control points.
- `X`: Array of new knots to insert.
- `r`: Number of new knots minus one.
- `Ubar`: New knot vector.
- `Qw`: New control points.

**Example:**
```javascript
let n = /* number of control points minus one */;
let p = /* degree of the curve */;
let U = /* original knot vector */;
let Pw = /* control points */;
let X = /* new knots to insert */;
let r = X.length - 1;
let Ubar = new Array(U.length + r + 1); // New knot vector
let Qw = new Array(Pw.length + r + 1); // New control points
RefineKnotVectCurve(n, p, U, Pw, X, r, Ubar, Qw);
console.log('New Knot Vector:', Ubar);
console.log('New Control Points:', Qw);
```

### Knot Refinement in NURBS Curve (Alternative Algorithm)

**Function:** `RefineKnotVectCurve(n, p, U, Pw, X, r)`

**Description:** Refines a NURBS curve by inserting a set of new knots into the original knot vector.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Original knot vector.
- `Pw`: Array of control points.
- `X`: Array of new knots to insert.
- `r`: Number of new knots minus one.

**Returns:** Object containing the new knot vector and new control points.

**Example:**
```javascript
let n = 3; // Example value for number of control points minus one
let p = 2; // Example value for degree of the curve
let U = [0, 0, 0, 1, 1, 1]; // Example original knot vector
let Pw = [/* Control points array */];
let X = [0.5, 0.75]; // Example array of new knots to be inserted
let r = X.length - 1;
let result = RefineKnotVectCurve(n, p, U, Pw, X, r);
console.log('New Knot Vector:', result.Ubar);
console.log('New Control Points:', result.Qw);
```

### Decomposing a NURBS Curve into Bézier Segments

**Function:** `DecomposeCurve(n, p, U, Pw)`

**Description:** Decomposes a NURBS curve into its constituent Bézier segments.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `Pw`: Array of control points.

**Returns:** Array of Bézier segments.

**Example:**
```javascript
let n = /* number of control points minus one */;
let p = /* degree of the curve */;
let U = /* knot vector */;
let Pw = /* control points array */;
let bezierSegments = DecomposeCurve(n, p, U, Pw);
console.log('Bézier Segments:', bezierSegments);
```

### Decomposing a NURBS Surface into Bézier Patches

**Function:** `DecomposeSurface(n, p, U, m, q, V, Pw, dir)`

**Description:** Decomposes a NURBS surface into its constituent Bézier patches in either the U or V direction.

**Parameters:**
- `n`: Number of control points in the U direction minus one.
- `p`: Degree in the U direction.
- `U`: Knot vector in the U direction.
- `m`: Number of control points in the V direction minus one.
- `q`: Degree in the V direction.
- `V`: Knot vector in the V direction.
- `Pw`: Array of control points.
- `dir`: Direction for decomposition (0 for U, 1 for V).

**Returns:** Array of Bézier patches.

**Example:**
```javascript
let n = /* number of control points in U minus one */;
let p = /* degree in U direction */;
let U = /* knot vector in U */;
let m = /* number of control points in V minus one */;
let q = /* degree in V direction */;
let V = /* knot vector in V */;
let Pw = /* control points array */;
let dir = /* direction for decomposition (0 for U, 1 for V) */;
let bezierPatches = DecomposeSurface(n, p, U, m, q, V, Pw, dir);
console.log('Bézier Patches:', bezierPatches);
```

### Removing a Knot from a NURBS Curve

**Function:** `RemoveCurveKnot(n, p, U, Pw, u, r, s, num)`

**Description:** Attempts to remove a knot from a NURBS curve a specified number of times.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `Pw`: Array of control points.
- `u`: Knot to be removed.
- `r`: Knot span index.
- `s`: Multiplicity of the knot.
- `num`: Number of times to remove the knot.

**Example:**
```javascript
let n, p, U, Pw, u, r, s, num;
// Initialize these variables as per your curve specifications
RemoveCurveKnot(n, p, U, Pw, u, r, s, num);
// The control points Pw and knot vector U are modified in place
```

### Degree Elevation of a NURBS Curve

**Function:** `DegreeElevateCurve(n, p, U, Pw, t)`

**Description:** Increases the degree of a NURBS curve by a specified factor.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Knot vector.
- `Pw`: Array of control points.
- `t`: Degree elevation factor.

**Returns:** Object containing the new degree, new knot vector, and new control points.

**Example:**
```javascript
let n = /* number of control points minus one */;
let p = /* degree of the curve */;
let U = /* original knot vector */;
let Pw = /* control points */;
let t = /* number of times to elevate degree */;
let { nh, Uh, Qw } = DegreeElevateCurve(n, p, U, Pw, t);
console.log('New degree:', nh);
console.log('New knot vector:', Uh);
console.log('New control points:', Qw);
```

### Degree Elevation of a NURBS Surface

**Function:** `DegreeElevateSurface(n, p, U, m, q, V, Pw, dir, t)`

**Description:** Increases the degree of a NURBS surface in either the U or V direction by a specified factor.

**Parameters:**
- `n`: Number of control points in the U direction minus one.
- `p`: Degree in the U direction.
- `U`: Knot vector in the U direction.
- `m`: Number of control points in the V direction minus one.
- `q`: Degree in the V direction.
- `V`: Knot vector in the V direction.
- `Pw`: Array of control points.
- `dir`: Direction for degree elevation (0 for U, 1 for V).
- `t`: Degree elevation factor.

**Returns:** Object containing the new degree, new knot vectors, and new control points.

**Example:**
```javascript
let n, p, U, m, q, V, Pw, dir, t;
// Initialize these variables as per your surface specifications
let result = DegreeElevateSurface(n, p, U, m, q, V, Pw, dir, t);
console.log('New degree in U direction:', result.nh);
console.log('New knot vector in U direction:', result.Uh);
console.log('New degree in V direction:', result.mh);
console.log('New knot vector in V direction:', result.Vh);
console.log('New control points:', result.Qw);
```

### Degree Reduction of a NURBS Curve

**Function:** `DegreeReduceCurve(n, p, U, Qw)`

**Description:** Attempts to reduce the degree of a NURBS curve from `p` to `p-1`.

**Parameters:**
- `n`: Number of control points minus one.
- `p`: Degree of the curve.
- `U`: Original knot vector.
- `Qw`: Array of control points.

**Returns:** Object containing the new degree, new knot vector, and new control points.

**Example:**
```javascript
let n = /* number of control points minus one */;
let p = /* degree of the curve */;
let U = /* original knot vector */;
let Qw = /* control points */;
let { nh, Uh, Pw } = DegreeReduceCurve(n, p, U, Qw);
console.log('New degree:', nh);
console.log('New knot vector:', Uh);
console.log('New control points:', Pw);
```

## Advanced Geometric Algorithms

### Matrix to Convert Bézier Form to Power Form

**Function:** `BezierToPowerMatrix(p)`

**Description:** Computes the matrix to convert a Bézier curve to its power form.

**Parameters:**
- `p`: Degree of the Bézier curve.

**Returns:** Matrix for converting Bézier form to power form.

**Example:**
```javascript
let p = /* degree of the Bézier curve */;
let matrix = BezierToPowerMatrix(p);
console.log(matrix);
```

### Matrix to Convert Power Form to Bézier Form

**Function:** `PowerToBezierMatrix(p, M)`

**Description:** Computes the matrix to convert a curve from its power form to the Bézier form.

**Parameters:**
- `p`: Degree of the curve.
- `M`: Matrix to be inverted (from Bézier to power form).

**Returns:** Inverse matrix for converting power form to Bézier form.

**Example:**
```javascript
let p = /* degree of the curve */;
let M = /* matrix from Bézier to power form (from A6.1) */;
let MI = PowerToBezierMatrix(p, M);
console.log(MI);
```

## Conics and Circles

### NURBS Circular Arc

**Function:** `MakeNurbsCircle(O, X, Y, r, ths, the)`

**Description:** Creates an arbitrary NURBS circular arc.

**Parameters:**
- `O`: Center point.
- `X`: Vector defining the x
