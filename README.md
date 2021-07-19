# Directions

This is a programming task as part of the interview process for a Software Engineer position at Metalenz. 

We have provided a skeleton module, unit tests, and a prompt using the Julia language below. The project is to complete the problems, filling in the necessary functions in `/src/DielectricRaster.jl` as needed. Note that **Question 6 is entirely optional**, and you are encouraged to refer to any documentation that you find useful.

You can check your work with the provided unit tests using the built-in Julia unit testing tool (`julia> ]test DielectricRaster`). You may add additional tests to the testing suite, and you should make use of any additional practices that you find useful.

You are allowed to use dependencies in addition to ones included, but no external packages are required. If you choose to include further dependencies, please document the packages and take reasonable steps to ensure your work is reproducible (i.e. update the `Project.toml`). You may use any version of Julia you prefer from 1.4 onwards.

Your submission should be in the form of a git repository. You can submit your project via email either by a repository URL, such as GitHub link, or a tarball containing the repository. 

This assignment is intended to be completable in four hours cumulative. 

If you encounter any issues, don't hesitate to reach out to Sam Buercklin (sam@metalenz.com) and Pawel Latawiec (pawel@metalenz.com) to clear up any problems.

As a final note, please take your time and don't worry if a problem is giving you trouble. We are interested in your skills as both a software engineer and as a computational developer; partial solutions count, as does your process towards arriving at an answer. 

# Rastering Dielectrics

For our simulations, we typically build up a system to solve from a set of primitives. The primitives give basic rules which govern how they are rastered onto an underlying array which represents pixels of a simulation.

We have the following:

```julia
module DielectricRaster

struct Square{T}
    x::Float64
    y::Float64
    w::Float64
    dielectric::T
end

struct GriddedArray{T}
    xlims::Tuple{Float64, Float64}
    ylims::Tuple{Float64, Float64}
    pixels::Matrix{T}
end

const TensorGriddedArray = GriddedArray{SMatrix{2, 2, Float64}}

end
```

A `Square` represents a `w`idth, a center `(x, y)` and some dielectric constant (related to the refractive index of a material squared). A `GriddedArray` represents a set of `pixels` from the limits represented by `xlims` and `ylims`. 

For instance, if we have `A = GriddedArray((0.0, 1.0), (0.0, 2.0), fill(1.0, 2, 2))`, then `A[1, 1]` represents the pixel with *center* at `(0.0, 0.0)`, lower-left corner `(-0.5, -1.0)`, upper-right corner `(+0.5, +1.0)`. `A[end, end]` represents the pixel with `center` at `(1.0, 2.0)`, lower-left corner `(+0.5, +1.0)`, upper-right corner `(+1.5, +3.0)`.
## Q1

For solving PDEs, we frequently wish to represent a set of rule-based primitives like `Square` on an underlying matrix. Write the function `raster_area!(A::GriddedArray, s::Square)` which populates each pixel in `A` with the fraction of area that pixel intersects with `s`. For instance, pixels entirely within the boundary of `Square` are set to `1.0`, pixels entirely outside of the boundary of `Square` are set to `0.0`, and pixels on the boundary are in-between.

## Q2

Now suppose the current value of `GriddedArray` represents an underlying dielectric. For each pixel find the average of the dielectric between the `GriddedArray` and the `Square` on each pixel: `\langle \varepsilon \rangle = E[epsilon]` where `E[...]` denotes the average or expectation of its argument.

## Q3

For rastering a dielectric, there is a more precise formula we can use. One part of this formula is to calculate the harmonic mean of the dielectric over the pixel. That is, we take the inverse of the average of the inverse of the dielectric. Write `raster_harmonic!(A::GriddedArray, s::Square)` which does that. That is to say, calculate `\langle \varepsilon^{-1} \rangle^{-1} = E[epsilon^-1]^-1`.

NOTE: Our `GriddedArray` must have non-zero entries, since a `0.0` dielectric constant is not well defined.

## Q4

Surfaces are important in E&M simulations. Write a function `surface_normal(s::Square, coord)` which returns the surface normal vector of the `s::Square` nearest to the `coord`inate. In the case where the nearest element of `s` is a corner (or vertex), return the vector which is the mean of the surface normals of the adjacent faces.

## Q5
`surface_normal` should return a `Vector` of some type. Write a function which transforms from the global coordinate frame used in `GriddedArray` to one where a vector returned by `surface_normal` is a basis vector, specifically the last basis vector. If passed a vector or array `A`, then the result of `orient_along_normal(A, normal)` should have the last axis parallel to `normal`.

## Q6 (Optional)

Now we can finally raster our dielectric. First of all, our `dielectric` must be made into a tensor- each pixel is a 2x2 matrix which represents the dielectric tensor. Specifically, in the local reference from of the surface normal (where the surface normal direction is the second component), the dielectric tensor for each pixel is given as

```
\begin{bmatrix}
\langle \varepsilon \rangle & 0 \\

0 &  \langle \varepsilon^{-1} \rangle^{-1}
\end{bmatrix}
```
or equivalently in pseudocode,
```
[ E[epsilon]        0; 
      0       E[epsilon^-1]^-1 ]

```
Write a function `raster_dielectric(A::TensorGriddedArray, B::GriddedArray, s::Square)` which correctly computes this tensor within the global reference frame. `B` represents the dielectric background media and `A` is the destination `TensorGriddedArray`.

First, for each pixel compute the values for `\langle \varepsilon \rangle` and `\langle \varepsilon^-1 \rangle^-1`. Next, construct the dielectric tensor in the **local** reference frame above using these values. Finally, transform the tensor at each pixel to the global reference frame using the tools you have developed. 


# Closing

Good luck! We're looking forward to going over your solution, and please reach out if you have any questions or concerns. 